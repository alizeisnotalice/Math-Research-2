import Foundation
import Metal
import Darwin

struct Input: Codable {
    let logEta, r, hR, bR2, a, logLambda, logK: Float
    let n: UInt32
}
struct Answer: Codable {
    let metric, t, realGap, firstGap, analyticGap, sourceGap, remainder, logAmplitude: Float
}
struct Row {
    let input: Input
    let answer: Answer
    let tested, realAbove, correctedAbove, firstDisagrees, nonfinite, index, salt, padding: UInt32
}
struct Hit: Codable {
    let input: Input
    let answer: Answer
    let index, salt: UInt32
}
struct Ledger: Codable {
    var tested = 0
    var realAbove = 0
    var correctedAbove = 0
    var firstDisagrees = 0
    var nonfinite = 0
    var best: Hit?

    mutating func update(_ row: Row) {
        tested += Int(row.tested)
        realAbove += Int(row.realAbove)
        correctedAbove += Int(row.correctedAbove)
        firstDisagrees += Int(row.firstDisagrees)
        nonfinite += Int(row.nonfinite)
        if row.answer.metric.isFinite && row.answer.metric > (best?.answer.metric ?? -.infinity) {
            best = Hit(input: row.input, answer: row.answer, index: row.index, salt: row.salt)
        }
    }
}
struct Report: Codable {
    let device: String
    let os: String
    let gridSeconds, commandSeconds, commandWaitSeconds: Double
    let batches, batchSize: Int
    let startUTC, endUTC: String
    let summaries: [String: Ledger]
}
let args = CommandLine.arguments
let requestedSeconds = args.count > 1 ? Double(args[1])! : 65
let suffix = args.count > 2 ? args[2] : "main"
let batchSize = 1 << 19
precondition(MemoryLayout<Row>.stride == 96)
guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
    fatalError("Missing Metal device")
}
let source = try String(contentsOfFile:"action.metal", encoding: .utf8)
let options = MTLCompileOptions()
options.mathMode = .safe
let library = try device.makeLibrary(source: source, options: options)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "action")!)
precondition(pipeline.threadExecutionWidth == 32)
let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
precondition(width % 32 == 0)
let keys = [4,6,8].flatMap { n in
    ["critical_asymmetric_small", "critical_symmetric_small",
     "critical_finite", "near_critical"].map { "n\(n)_\($0)" }
}
let rowCount = batchSize / 32 * keys.count
let output = device.makeBuffer(length: rowCount * MemoryLayout<Row>.stride, options: .storageModeShared)!
let rows = output.contents().bindMemory(to: Row.self, capacity: rowCount)
let iso = ISO8601DateFormatter()
iso.formatOptions = [.withInternetDateTime, .withFractionalSeconds]
let startUTC = iso.string(from: Date())
print("start=\(startUTC) device=\(device.name) group=\(width) pid=\(getpid())")
fflush(stdout)

let start = ProcessInfo.processInfo.systemUptime
var waitSeconds = 0.0
var gpuSeconds = 0.0
var batches = 0
var summaries = Array(repeating: Ledger(), count: keys.count)
repeat {
    let base = UInt64(batches) * UInt64(batchSize)
    let salt = UInt32(base >> 32)
    var parameters = SIMD2(UInt32(truncatingIfNeeded: base), salt)
    let command = queue.makeCommandBuffer()!
    let encoder = command.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(output, offset: 0, index: 0)
    encoder.setBytes(&parameters, length: MemoryLayout<SIMD2<UInt32>>.stride, index: 1)
    encoder.dispatchThreads(MTLSize(width: batchSize, height: 1, depth: 1),
                            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
    encoder.endEncoding()
    let beginWait = ProcessInfo.processInfo.systemUptime
    command.commit()
    command.waitUntilCompleted()
    waitSeconds += ProcessInfo.processInfo.systemUptime - beginWait
    guard command.status == .completed else {
        fatalError("Metal command failed: \(String(describing: command.error))")
    }
    gpuSeconds += command.gpuEndTime - command.gpuStartTime
    for i in 0..<rowCount {
        summaries[i % keys.count].update(rows[i])
    }
    batches += 1
    if batches == 1 || batches % 25 == 0 {
        print("utc=\(iso.string(from: Date())) batches=\(batches) evaluated=\(batches*batchSize) gpu_seconds=\(gpuSeconds)")
        fflush(stdout)
    }
} while ProcessInfo.processInfo.systemUptime - start < requestedSeconds

let total = ProcessInfo.processInfo.systemUptime - start
let endUTC = iso.string(from: Date())
let report = Report(device: device.name, os: ProcessInfo.processInfo.operatingSystemVersionString,
                    gridSeconds: total, commandSeconds: gpuSeconds,
                    commandWaitSeconds: waitSeconds, batches: batches,
                    batchSize: batchSize, startUTC: startUTC, endUTC: endUTC,
                    summaries: Dictionary(uniqueKeysWithValues: zip(keys, summaries)))
let json = JSONEncoder()
json.outputFormatting = [.prettyPrinted, .sortedKeys]
try json.encode(report).write(to: URL(fileURLWithPath: "action-\(suffix).json"), options: .atomic)
print("end=\(endUTC) batches=\(batches) evaluated=\(batches*batchSize) gpu_seconds=\(gpuSeconds) wall_seconds=\(total)")
