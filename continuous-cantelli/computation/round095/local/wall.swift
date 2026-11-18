import Foundation
import Metal
import Darwin
struct Point: Codable {
    let k,v: Float
    let family: UInt32
}
struct Row {
    let p: Point
    let logR,logLinear,oldDelta,upperDelta,lowerDelta,score: Float
    let tested,oldAbove,upperAbove,nonfinite,index,salt: UInt32
}
struct Hit: Codable {
    let p: Point
    let logR,logLinear,oldDelta,upperDelta,lowerDelta,score: Float
    let index,salt: UInt32
    init(_ row: Row) {
        p=row.p
        logR=row.logR
        logLinear=row.logLinear
        oldDelta=row.oldDelta
        upperDelta=row.upperDelta
        lowerDelta=row.lowerDelta
        score=row.score
        index=row.index
        salt=row.salt
    }
}
struct Ledger: Codable {
    var tested=0
    var nonfinite=0
    var oldAbove=0
    var upperAbove=0
    var best: Hit?
    mutating func update(_ row: Row) {
        tested += Int(row.tested)
        nonfinite += Int(row.nonfinite)
        oldAbove += Int(row.oldAbove)
        upperAbove += Int(row.upperAbove)
        if row.score.isFinite && row.score > (best?.score ?? -.infinity) { best=Hit(row) }
    }
}
struct Report: Codable {
    let device: String
    let os: String
    let order: Int
    let gridSeconds: Double
    let commandSeconds: Double
    let commandWaitSeconds: Double
    let batches: Int
    let batchSize: Int
    let startUTC: String
    let endUTC: String
    let startIndex, endIndexExclusive: UInt64
    let summaries: [String: Ledger]
}

func gaussLegendre(_ count: Int) -> [SIMD2<Float>] {
    var nodes = [SIMD2<Float>](repeating: SIMD2.zero, count: count)
    for i in 0..<(count/2) {
        var z = cos(Double.pi * (Double(i) + 0.75) / (Double(count) + 0.5))
        var derivative = 0.0
        for _ in 0..<12 {
            var previous = 0.0
            var value = 1.0
            for j in 1...count {
                let secondPrevious = previous
                previous = value
                value = ((2 * Double(j) - 1) * z * previous
                         - (Double(j) - 1) * secondPrevious) / Double(j)
            }
            derivative = Double(count) * (z * value - previous) / (z * z - 1)
            let adjustment = value / derivative
            z -= adjustment
            if abs(adjustment) < 2e-16 { break }
        }
        let weight = Float(1 / ((1-z*z) * derivative*derivative))
        nodes[i] = SIMD2(Float((1-z)/2), weight)
        nodes[count-1-i] = SIMD2(Float((1+z)/2), weight)
    }
    return nodes
}

let args = CommandLine.arguments
let requestedSeconds = args.count > 1 ? Double(args[1])! : 65
let order = args.count > 2 ? Int(args[2])! : 96
let suffix = args.count > 3 ? args[3] : "main"
let batchSize = 1 << 17
let startIndex = args.count > 5 ? UInt64(args[5])! : 0
precondition(startIndex % UInt64(batchSize) == 0)
precondition(MemoryLayout<Row>.stride == 60)
guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
    fatalError("Missing Metal device")
}
let source = try String(contentsOfFile:"wall.metal", encoding: .utf8)
let options = MTLCompileOptions()
options.mathMode = .safe
let library = try device.makeLibrary(source: source, options: options)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "wall")!)
precondition(pipeline.threadExecutionWidth == 32)
let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
precondition(width % 32 == 0)
let nodes = gaussLegendre(order)
let nodeBuffer = nodes.withUnsafeBytes {
    device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
}
let keys = ["all","far_weak_wall","rational_upper","interior","transition",
            "small_wall","naive_bound_failures","strong_subcritical","near_observation"]
let rowCount = batchSize / 32 * keys.count
let output = device.makeBuffer(length: rowCount * MemoryLayout<Row>.stride, options: .storageModeShared)!
let rows = output.contents().bindMemory(to: Row.self, capacity: rowCount)
let iso = ISO8601DateFormatter()
iso.formatOptions = [.withInternetDateTime, .withFractionalSeconds]
if args.count > 4, let timestamp = Double(args[4]), timestamp > 0 {
    let wait = timestamp - Date().timeIntervalSince1970
    if wait > 0 { Thread.sleep(forTimeInterval: wait) }
}
let startUTC = iso.string(from: Date())
print("start=\(startUTC) device=\(device.name) order=\(order) first_index=\(startIndex) group=\(width) pid=\(getpid())")
fflush(stdout)

let start = ProcessInfo.processInfo.systemUptime
var waitSeconds = 0.0
var gpuSeconds = 0.0
var batches = 0
var summaries = Array(repeating: Ledger(), count: keys.count)
repeat {
    let base = startIndex + UInt64(batches) * UInt64(batchSize)
    let salt = UInt32(base >> 32)
    var parameters = SIMD4(UInt32(batchSize), UInt32(order), UInt32(truncatingIfNeeded: base), salt)
    let command = queue.makeCommandBuffer()!
    let encoder = command.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(nodeBuffer, offset: 0, index: 0)
    encoder.setBuffer(output, offset: 0, index: 1)
    encoder.setBytes(&parameters, length: MemoryLayout<SIMD4<UInt32>>.stride, index: 2)
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
        print("utc=\(iso.string(from: Date())) batches=\(batches) tested=\(summaries[0].tested) gpu_seconds=\(gpuSeconds) max_fixed_wall_deviation=\(summaries[1].best?.score ?? .nan)")
        fflush(stdout)
    }
} while ProcessInfo.processInfo.systemUptime - start < requestedSeconds

let total = ProcessInfo.processInfo.systemUptime - start
precondition(summaries[0].tested == batches*batchSize)
let endUTC = iso.string(from: Date())
let report = Report(device: device.name, os: ProcessInfo.processInfo.operatingSystemVersionString,
                    order: order, gridSeconds: total, commandSeconds: gpuSeconds,
                    commandWaitSeconds: waitSeconds, batches: batches,
                    batchSize: batchSize, startUTC: startUTC, endUTC: endUTC,
                    startIndex: startIndex, endIndexExclusive: startIndex + UInt64(batches)*UInt64(batchSize),
                    summaries: Dictionary(uniqueKeysWithValues: zip(keys, summaries)))
let json = JSONEncoder()
json.outputFormatting = [.prettyPrinted, .sortedKeys]
try json.encode(report).write(to: URL(fileURLWithPath: "wall-\(suffix).json"), options: .atomic)
print("end=\(endUTC) batches=\(batches) tested=\(summaries[0].tested) gpu_seconds=\(gpuSeconds) wall_seconds=\(total)")
