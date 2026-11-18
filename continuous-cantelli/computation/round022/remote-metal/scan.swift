import Foundation
import Metal
import Darwin

struct Core {
    var x: Float
    var c: Float
    var d: Float
    var radius: Float
    var variance: Float
}

struct Settings {
    var shiftCount: UInt32
    var epsilonCount: UInt32
    var panelCount: UInt32
    var epsilonMin: Float
    var epsilonMax: Float
}

struct Parameters {
    var t: Double
    var x: Double
    var c: Double
    var d: Double
    var radius: Double
    var variance: Double
}

struct Hit: Codable {
    var ratio: Float
    var t: Double
    var x: Double
    var c: Double
    var d: Double
    var radius: Double
    var variance: Double
    var shift: Double
    var epsilon: Double
}

func polynomial(_ t: Double) -> Double {
    ((5913 * t + 954) * t + 96) * t - 8
}

func logGrid(_ low: Double, _ high: Double, _ i: Int, _ count: Int) -> Double {
    exp(log(low) + (log(high) - log(low)) * Double(i) / Double(count - 1))
}

func checkMoments(_ p: Parameters) -> Bool {
    let n = 64
    let spacing = p.radius / Double(n)
    let normalization = sqrt(2 * Double.pi)
    var coreMass = 0.0
    var coreV = 0.0
    for i in 0...n {
        let u = Double(i) * spacing
        let left = exp(-0.5 * pow(p.x - u, 2)) / normalization
        let right = exp(-0.5 * pow(p.x + u, 2)) / normalization
        let weight = i == 0 || i == n ? 1.0 : (i.isMultiple(of: 2) ? 2.0 : 4.0)
        coreMass += weight * (left + right)
        coreV += weight * (left * p.c * u + right * p.d * pow(u, 3))
    }
    coreMass *= spacing / 3
    coreV *= spacing / 3
    return coreV <= p.variance - 1 + 1e-9 && p.variance * coreMass - coreV <= 1 + 1e-9
}

func makeGrid(tCount: Int, mCount: Int, radius: Double, mMin: Double, mMax: Double) -> [Parameters] {
    var a = 0.0
    var b = 0.1
    for _ in 0..<55 {
        let m = (a + b) / 2
        if polynomial(m) > 0 { b = m } else { a = m }
    }
    let root = (a + b) / 2
    var grid = [Parameters]()
    for ti in 0..<tCount {
        let t = root + logGrid(1e-4, 30, ti, tCount)
        let c = sqrt(polynomial(t) / (21 * t + 6))
        let d = t / c
        let x = 3 * d - 2 / c
        for mi in 0..<mCount {
            let m = 1 + logGrid(mMin - 1, mMax - 1, mi, mCount)
            let r = min(radius, m / c, pow(m / d, 1.0 / 3.0))
            let p = Parameters(t: t, x: x, c: c, d: d, radius: r, variance: m)
            if checkMoments(p) { grid.append(p) }
        }
    }
    return grid
}

func run() throws {
    let args = CommandLine.arguments
    guard args.count == 8 || args.count == 13 else {
        throw NSError(domain: "Usage: scan t-count m-count shift-count epsilon-count panels chunk-cores output [radius m-min m-max epsilon-min epsilon-max]", code: 1)
    }
    let n = args[1...6].compactMap(Int.init)
    guard n.count == 6, n[0...3].allSatisfy({ $0 >= 2 }), n[4] > 0, n[5] > 0 else {
        throw NSError(domain: "Invalid grid dimensions", code: 2)
    }
    let defaults = [0.25, 1.0001, 6.0, 1e-5, 1.0]
    let bounds = args.count == 13 ? args[8...12].compactMap(Double.init) : defaults
    guard bounds.count == 5, bounds[0] > 0, bounds[1] > 1, bounds[2] > bounds[1], bounds[3] > 0, bounds[4] > bounds[3] else {
        throw NSError(domain: "Invalid grid bounds", code: 7)
    }
    let grid = makeGrid(tCount: n[0], mCount: n[1], radius: bounds[0], mMin: bounds[1], mMax: bounds[2])
    guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
        throw NSError(domain: "Metal device or queue unavailable", code: 3)
    }
    let stable = ProcessInfo.processInfo.environment["CC20_STABLE"] == "1"
    let source = try String(contentsOfFile: stable ? "density_stable.metal" : "density.metal", encoding: .utf8)
    let options = MTLCompileOptions()
    options.mathMode = .safe
    let library = try device.makeLibrary(source: source, options: options)
    guard let kernel = library.makeFunction(name: "scan_density") else {
        throw NSError(domain: "Kernel unavailable", code: 4)
    }
    let pipeline = try device.makeComputePipelineState(function: kernel)
    var settings = Settings(shiftCount: UInt32(n[2]), epsilonCount: UInt32(n[3]), panelCount: UInt32(n[4]), epsilonMin: Float(bounds[3]), epsilonMax: Float(bounds[4]))
    let pointsPerCore = n[2] * n[3]
    let threadCount = min(pipeline.maxTotalThreadsPerThreadgroup, 4 * pipeline.threadExecutionWidth)
    let floatStride = MemoryLayout<Float>.stride
    var best = [Hit]()
    var bestByRegime = [String: Hit]()
    var timeSum = 0.0
    var overflowCount = 0
    var validCount = 0
    print("device=\(device.name) cores=\(grid.count) panels=\(n[4]) shifts=\(n[2]) epsilons=\(n[3]) threads_per_group=\(threadCount)")
    fflush(stdout)
    for start in stride(from: 0, to: grid.count, by: n[5]) {
        let end = min(start + n[5], grid.count)
        let tile = grid[start..<end].map { p in
            Core(x: Float(p.x), c: Float(p.c), d: Float(p.d), radius: Float(p.radius), variance: Float(stable ? p.variance - 1 : p.variance))
        }
        let pointCount = tile.count * pointsPerCore
        guard let input = device.makeBuffer(bytes: tile, length: tile.count * MemoryLayout<Core>.stride, options: .storageModeShared),
              let output = device.makeBuffer(length: pointCount * floatStride, options: .storageModeShared),
              let command = queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
            throw NSError(domain: "Metal allocation failed", code: 5)
        }
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(input, offset: 0, index: 0)
        encoder.setBytes(&settings, length: MemoryLayout<Settings>.stride, index: 1)
        encoder.setBuffer(output, offset: 0, index: 2)
        encoder.dispatchThreads(MTLSize(width: pointCount, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: threadCount, height: 1, depth: 1))
        encoder.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else { throw command.error ?? NSError(domain: "Metal dispatch failed", code: 6) }
        let gpuSeconds = command.gpuEndTime - command.gpuStartTime
        timeSum += gpuSeconds
        let values = output.contents().bindMemory(to: Float.self, capacity: pointCount)
        var tileBest: Float = 0
        for i in 0..<pointCount {
            let ratio = values[i]
            guard ratio.isFinite else { overflowCount += 1; continue }
            validCount += 1
            tileBest = max(tileBest, ratio)
            let core = grid[start + i / pointsPerCore]
            let regime = core.x < -2 ? "xBelowMinusTwo" : (core.x < 0 ? "xBetweenMinusTwoAndZero" : "xNonnegative")
            let beatsRegime = (bestByRegime[regime]?.ratio ?? 0) < ratio
            if !beatsRegime && best.count == 12 && best[11].ratio >= ratio { continue }
            let rem = i % pointsPerCore
            let shift = core.radius * (2 * Double(rem / n[3]) / Double(n[2] - 1) - 1)
            let epsilon = logGrid(Double(settings.epsilonMin), Double(settings.epsilonMax), rem % n[3], n[3])
            let hit = Hit(ratio: ratio, t: core.t, x: core.x, c: core.c, d: core.d, radius: core.radius, variance: core.variance, shift: shift, epsilon: epsilon)
            if beatsRegime { bestByRegime[regime] = hit }
            best.append(hit)
            best.sort { $0.ratio > $1.ratio }
            if best.count > 12 { best.removeLast() }
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(gpuSeconds) max_ratio=\(tileBest)")
        fflush(stdout)
    }
    let encoder = JSONEncoder()
    encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
    let results: [String: Any] = [
        "device": device.name, "validCount": validCount, "overflowCount": overflowCount,
        "summedGPUSeconds": timeSum, "coreCount": grid.count, "dimensions": Array(n[0...4]), "bounds": bounds, "stableDifferenceOfSquares": stable,
        "best": try JSONSerialization.jsonObject(with: encoder.encode(best)),
        "bestByRegime": try JSONSerialization.jsonObject(with: encoder.encode(bestByRegime))
    ]
    let data = try JSONSerialization.data(withJSONObject: results, options: [.prettyPrinted, .sortedKeys])
    try data.write(to: URL(fileURLWithPath: args[7]), options: .atomic)
    print("valid=\(validCount) overflow=\(overflowCount) gpu_sum_seconds=\(timeSum) best_ratio=\(best.first?.ratio ?? 0)")
}

try run()
