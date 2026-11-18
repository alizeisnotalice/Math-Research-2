import Foundation
import Metal
import Darwin

struct Core {
    var center: Float
    var coefficient: Float
    var radius: Float
    var excessVariance: Float
}

struct Settings {
    var shiftCount: UInt32
    var epsilonCount: UInt32
    var panels: UInt32
}

struct Parameters {
    var center: Double
    var coefficient: Double
    var radius: Double
    var variance: Double
    var positionGap: Double
    var endpointFraction: Double
    var coreMass: Double
    var coreVarianceMoment: Double
}

struct Hit: Codable {
    var ratio: Float
    var center: Double
    var coefficient: Double
    var radius: Double
    var variance: Double
    var positionGap: Double
    var endpointFraction: Double
    var coreMass: Double
    var coreVarianceMoment: Double
    var shift: Double
    var epsilon: Double
}

func logGrid(_ low: Double, _ high: Double, _ i: Int, _ count: Int) -> Double {
    exp(log(low) + (log(high) - log(low)) * Double(i) / Double(count - 1))
}

func coreMoments(center: Double, radius: Double) -> (Double, Double) {
    let panels = 64
    let normalization = sqrt(2 * Double.pi)
    var mass = 0.0
    var unitQuartic = 0.0
    for i in 0...panels {
        let z = Double(i) / Double(panels)
        let u = z * radius
        let left = center - u
        let right = center + u
        let densityPair = (exp(-0.5 * left * left) + exp(-0.5 * right * right)) / normalization
        let weight = i == 0 || i == panels ? 1.0 : (i.isMultiple(of: 2) ? 2.0 : 4.0)
        mass += weight * densityPair
        unitQuartic += weight * densityPair * z * z * z * z
    }
    let step = radius / (3 * Double(panels))
    return (mass * step, unitQuartic * step)
}

func makeGrid(_ n: [Int]) -> [Parameters] {
    var grid = [Parameters]()
    for mi in 0..<n[0] {
        let m = 1 + logGrid(1.0 / 1024, 11, mi, n[0])
        let threshold = sqrt(m * log(m) / (m - 1))
        for ai in 0..<n[1] {
            let gap = logGrid(1e-5, 4, ai, n[1])
            let a = threshold + gap
            for ri in 0..<n[2] {
                let r = logGrid(1.0 / 32, 1.5, ri, n[2])
                let r4 = r * r * r * r
                let (mass, unitQuartic) = coreMoments(center: a, radius: r)
                for bi in 0..<n[3] {
                    let fraction = logGrid(1.0 / 256, 1, bi, n[3])
                    let coreVariance = fraction * m * unitQuartic
                    if coreVariance <= m - 1 + 1e-10 && m * mass - coreVariance <= 1 + 1e-10 {
                        grid.append(Parameters(center: a, coefficient: fraction * m / r4, radius: r, variance: m,
                                               positionGap: gap, endpointFraction: fraction, coreMass: mass,
                                               coreVarianceMoment: coreVariance))
                    }
                }
            }
        }
    }
    return grid
}

func run() throws {
    let args = CommandLine.arguments
    guard args.count == 10 else {
        throw NSError(domain: "Usage: scan m-count gap-count r-count endpoint-fraction-count shift-count epsilon-count panels chunk-cores output", code: 1)
    }
    let n = args[1...8].compactMap(Int.init)
    guard n.count == 8, n[0...5].allSatisfy({ $0 >= 2 }), n[6] > 0, n[7] > 0 else {
        throw NSError(domain: "Invalid grid dimensions", code: 2)
    }
    let grid = makeGrid(n)
    guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
        throw NSError(domain: "Metal device or queue unavailable", code: 3)
    }
    let source = try String(contentsOfFile: "quartic.metal", encoding: .utf8)
    let options = MTLCompileOptions()
    options.mathMode = .safe
    let library = try device.makeLibrary(source: source, options: options)
    guard let kernel = library.makeFunction(name: "quartic_scan") else {
        throw NSError(domain: "Kernel unavailable", code: 4)
    }
    let pipeline = try device.makeComputePipelineState(function: kernel)
    var settings = Settings(shiftCount: UInt32(n[4]), epsilonCount: UInt32(n[5]), panels: UInt32(n[6]))
    let pointsPerCore = n[4] * n[5]
    let threadCount = min(pipeline.maxTotalThreadsPerThreadgroup, 4 * pipeline.threadExecutionWidth)
    let floatStride = MemoryLayout<Float>.stride
    var bestByRegime = [String: Hit]()
    var timeSum = 0.0
    var overflowCount = 0
    var validCount = 0
    print("device=\(device.name) cores=\(grid.count) panels=\(n[6]) shifts=\(n[4]) epsilons=\(n[5]) threads_per_group=\(threadCount)")
    fflush(stdout)
    for start in stride(from: 0, to: grid.count, by: n[7]) {
        let end = min(start + n[7], grid.count)
        let tile = grid[start..<end].map { p in
            Core(center: Float(p.center), coefficient: Float(p.coefficient), radius: Float(p.radius), excessVariance: Float(p.variance - 1))
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
        for ci in 0..<tile.count {
            let core = grid[start + ci]
            let coreStart = ci * pointsPerCore
            var maximum: Float = 0
            var bestIndex = 0
            for ri in 0..<pointsPerCore {
                let ratio = values[coreStart + ri]
                if !ratio.isFinite { overflowCount += 1; continue }
                validCount += 1
                if ratio > maximum { maximum = ratio; bestIndex = ri }
            }
            tileBest = max(maximum, tileBest)
            var regimes = ["overall"]
            if core.radius <= 0.25001 { regimes.append("radiusAtMostQuarter") }
            if core.radius <= 0.12501 { regimes.append("radiusAtMostEighth") }
            if core.positionGap <= 0.01 { regimes.append("closeToPositionalBoundary") }
            if core.variance >= 2 { regimes.append("varianceAtLeastTwo") }
            if core.center <= 2 { regimes.append("centerAtMostTwo") }
            for regime in regimes where maximum > (bestByRegime[regime]?.ratio ?? 0) {
                let shift = 1.25 * core.radius * (2 * Double(bestIndex / n[5]) / Double(n[4] - 1) - 1)
                let epsilon = logGrid(1e-7, 0.25, bestIndex % n[5], n[5])
                bestByRegime[regime] = Hit(ratio: maximum, center: core.center, coefficient: core.coefficient, radius: core.radius,
                    variance: core.variance, positionGap: core.positionGap, endpointFraction: core.endpointFraction,
                    coreMass: core.coreMass, coreVarianceMoment: core.coreVarianceMoment, shift: shift, epsilon: epsilon)
            }
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(gpuSeconds) max_ratio=\(tileBest)")
        fflush(stdout)
    }
    let encoder = JSONEncoder()
    encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
    let results: [String: Any] = ["device": device.name, "validCount": validCount, "overflowCount": overflowCount,
        "summedGPUSeconds": timeSum, "coreCount": grid.count, "dimensions": Array(n[0...6]),
        "bestByRegime": try JSONSerialization.jsonObject(with: encoder.encode(bestByRegime))]
    let data = try JSONSerialization.data(withJSONObject: results, options: [.prettyPrinted, .sortedKeys])
    try data.write(to: URL(fileURLWithPath: args[9]), options: .atomic)
    print("valid=\(validCount) overflow=\(overflowCount) gpu_sum_seconds=\(timeSum) best_ratio=\(bestByRegime["overall"]?.ratio ?? 0)")
}

try run()
