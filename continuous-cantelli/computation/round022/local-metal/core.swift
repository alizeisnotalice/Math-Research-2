import Foundation
import Metal
import Darwin

struct Candidate {
    let m: Float
    let r: Float
    let epsilon: Float
    let delta: Float
}

struct Hit: Codable {
    let ratio: Float
    let M: Float
    let r: Float
    let epsilon: Float
    let delta: Float
}

struct Ledger: Codable {
    var tested = 0
    var aboveOne = 0
    var nonFinite = 0
    var hit: Hit?

    mutating func update(_ candidate: Candidate, value: Float) {
        tested += 1
        guard value.isFinite else {
            nonFinite += 1
            return
        }
        aboveOne += value > 1 ? 1 : 0
        if value > (hit?.ratio ?? -.infinity) {
            hit = Hit(ratio: value, M: candidate.m, r: candidate.r,
                      epsilon: candidate.epsilon, delta: candidate.delta)
        }
    }
}

struct Report: Codable {
    let device: String
    let os: String
    let arguments: [String]
    let tested: Int
    let metalSeconds: Double
    let commandWaitSeconds: Double
    let totalSeconds: Double
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

func geometric(_ minimumPower: Double, _ maximumPower: Double, _ count: Int) -> [Double] {
    (0..<count).map {
        pow(10, minimumPower + Double($0) / Double(count-1) * (maximumPower-minimumPower))
    }
}

let args = CommandLine.arguments
let mode = args.count > 1 ? args[1] : "full"
let directory = args.count > 2 ? args[2] : "."
let order = args.count > 3 ? Int(args[3])! : 64
let pwd = URL(fileURLWithPath: directory, isDirectory: true)
let outputURL = pwd.appendingPathComponent("metal-\(mode)-\(order).json")
guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
    fatalError("Missing Metal device")
}
let source = try String(contentsOf: pwd.appendingPathComponent("core.metal"), encoding: .utf8)
let options = MTLCompileOptions()
options.mathMode = .safe
let library = try device.makeLibrary(source: source, options: options)
guard let kernel = library.makeFunction(name: "core_ratio") else {
    fatalError("Missing Metal kernel")
}
let pipeline = try device.makeComputePipelineState(function: kernel)
let nodes = gaussLegendre(order)
let nodeBuffer = nodes.withUnsafeBytes {
    device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
}
let batchSize = 1 << 19
let input = device.makeBuffer(length: batchSize * MemoryLayout<Candidate>.stride, options: .storageModeShared)!
let output = device.makeBuffer(length: batchSize * MemoryLayout<Float>.stride, options: .storageModeShared)!
let candidates = input.contents().bindMemory(to: Candidate.self, capacity: batchSize)
let values = output.contents().bindMemory(to: Float.self, capacity: batchSize)
let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
print("device=\(device.name) mode=\(mode) order=\(order) group=\(width) pid=\(getpid())")
fflush(stdout)

let start = ProcessInfo.processInfo.systemUptime
var waitSeconds = 0.0
var gpuSeconds = 0.0
var tested = 0
var count = 0
var summaries: [String: Ledger] = [:]

func runBatch() {
    if count == 0 { return }
    let command = queue.makeCommandBuffer()!
    let encoder = command.makeComputeCommandEncoder()!
    var counts = SIMD2(UInt32(count), UInt32(order))
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(input, offset: 0, index: 0)
    encoder.setBuffer(nodeBuffer, offset: 0, index: 1)
    encoder.setBuffer(output, offset: 0, index: 2)
    encoder.setBytes(&counts, length: MemoryLayout<SIMD2<UInt32>>.stride, index: 3)
    encoder.dispatchThreads(MTLSize(width: count, height: 1, depth: 1),
                            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
    encoder.endEncoding()
    let waitStart = ProcessInfo.processInfo.systemUptime
    command.commit()
    command.waitUntilCompleted()
    waitSeconds += ProcessInfo.processInfo.systemUptime - waitStart
    if command.status != .completed {
        fatalError("Metal command failed: \(String(describing: command.error))")
    }
    gpuSeconds += command.gpuEndTime - command.gpuStartTime
    for i in 0..<count {
        let c = candidates[i]
        let ratio = values[i]
        summaries["all", default: Ledger()].update(c, value: ratio)
        if abs(c.delta) <= 1 {
            summaries["abs_y_minus_1_le_1", default: Ledger()].update(c, value: ratio)
        }
        if c.m <= 1.01 && abs(c.delta) <= 1 {
            summaries["near_M1_and_contact", default: Ledger()].update(c, value: ratio)
        }
        if c.m <= 1.01 && abs(c.delta) <= 1 && c.r <= 0.25 {
            summaries["near_M1_contact_r_le_0.25", default: Ledger()].update(c, value: ratio)
        }
        if abs(c.delta) <= 1e-3 && c.r <= 0.25 {
            summaries["near_contact_r_le_0.25", default: Ledger()].update(c, value: ratio)
        }
        if c.m <= 1.01 && abs(c.delta) <= 1e-3 && c.r <= 0.25 && c.epsilon <= 1e-8 {
            summaries["near_M1_near_contact_small_core_small_epsilon", default: Ledger()].update(c, value: ratio)
        }
        if c.m <= 1.01 && abs(c.delta) <= 1e-3 && c.r <= 0.25 && c.epsilon <= 1e-8 && c.delta == 0 {
            summaries["near_M1_at_contact_small_core_small_epsilon", default: Ledger()].update(c, value: ratio)
        }
    }
    tested += count
    count = 0
    print("tested=\(tested) cumulative_metal_seconds=\(gpuSeconds) best=\(summaries["all"]?.hit?.ratio ?? .nan)")
    fflush(stdout)
}

func submit(_ m: Double, _ r: Double, _ epsilon: Double, _ delta: Double) {
    let mf = Float(m)
    let radiusLimit = min(Double(mf), pow(Double(mf), 1.0 / 3.0))
    var rf = Float(r)
    if Double(rf) > radiusLimit { rf = Float(radiusLimit).nextDown }
    candidates[count] = Candidate(m: mf, r: rf, epsilon: Float(epsilon), delta: Float(delta))
    count += 1
    if count == batchSize { runBatch() }
}

let mOffsets = geometric(-6, log10(31), mode == "smoke" ? 5 : 72)
let fractions = mode == "smoke" ? [0.025, 0.25, 0.5, 0.75, 1.0] :
    geometric(-4, -1, 8) + (1...32).map { Double($0)/32 }
let epsilons = geometric(-12, 1, mode == "smoke" ? 8 : 64)
let contactDeltas = geometric(-8, -0.025, mode == "smoke" ? 9 : 56)
let standardized = mode == "smoke" ? [-12.0, -9, -6, -3, -1, 0, 1, 3, 6, 9, 12] :
    (-72...72).map { Double($0) / 6 }
for offset in mOffsets {
    let m = 1 + offset
    let maximumRadius = min(m, pow(m, 1.0 / 3.0))
    for epsilon in epsilons {
        let scale = sqrt(m + epsilon)
        for fraction in fractions {
            let r = maximumRadius * fraction
            submit(m, r, epsilon, 0)
            for delta in contactDeltas {
                submit(m, r, epsilon, delta)
                submit(m, r, epsilon, -delta)
            }
            for y in standardized {
                submit(m, r, epsilon, y * scale - 1)
            }
        }
    }
}
runBatch()
let total = ProcessInfo.processInfo.systemUptime - start
let report = Report(device: device.name, os: ProcessInfo.processInfo.operatingSystemVersionString,
                    arguments: Array(args.dropFirst()), tested: tested, metalSeconds: gpuSeconds,
                    commandWaitSeconds: waitSeconds, totalSeconds: total, summaries: summaries)
let json = JSONEncoder()
json.outputFormatting = [.prettyPrinted, .sortedKeys]
try json.encode(report).write(to: outputURL, options: .atomic)
