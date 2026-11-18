import CryptoKit
import Foundation
import Metal

let kernelSource = """
#include <metal_stdlib>
using namespace metal;

kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &stride [[buffer(2)]],
                       uint index [[thread_position_in_grid]]) {
    int own = input[index];
    int other = input[index ^ stride];
    output[index] = (index & stride) == 0 ? own + other : other - own;
}
"""

struct CheckFailure: Error {
    let reason: String
}

func check(_ condition: Bool, _ reason: String) throws {
    guard condition else { throw CheckFailure(reason: reason) }
}

func hash(_ data: Data) -> String {
    SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func parity(_ word: Int) -> Int {
    word.nonzeroBitCount & 1
}

func sign(_ bit: Int) -> Int32 {
    bit == 0 ? 1 : -1
}

func directWalsh(_ table: [Int32]) -> [Int32] {
    (0..<table.count).map { frequency in
        table.indices.reduce(into: Int32(0)) { sum, point in
            sum += table[point] * sign(parity(frequency & point))
        }
    }
}

final class MetalWalsh {
    let device: MTLDevice
    let queue: MTLCommandQueue
    let pipeline: MTLComputePipelineState
    var dispatches = 0
    var gpuSeconds: Double = 0

    init() throws {
        guard let device = MTLCreateSystemDefaultDevice(),
              let queue = device.makeCommandQueue()
        else { throw CheckFailure(reason: "Metal device or command queue unavailable") }
        let library = try device.makeLibrary(source: kernelSource, options: nil)
        guard let function = library.makeFunction(name: "fwht_stage") else {
            throw CheckFailure(reason: "Metal FWHT shader unavailable")
        }
        self.device = device
        self.queue = queue
        self.pipeline = try device.makeComputePipelineState(function: function)
    }

    func transform(_ input: [Int32], rowSize: Int) throws -> [Int32] {
        let bytes = input.count * MemoryLayout<Int32>.stride
        let first = input.withUnsafeBytes {
            device.makeBuffer(bytes: $0.baseAddress!, length: bytes, options: .storageModeShared)
        }
        guard var read = first,
              var write = device.makeBuffer(length: bytes, options: .storageModeShared),
              let command = queue.makeCommandBuffer()
        else { throw CheckFailure(reason: "Metal buffer allocation failed") }
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        var stride = 1
        while stride < rowSize {
            guard let encoder = command.makeComputeCommandEncoder() else {
                throw CheckFailure(reason: "Metal command encoder unavailable")
            }
            var step = UInt32(stride)
            encoder.setComputePipelineState(pipeline)
            encoder.setBuffer(read, offset: 0, index: 0)
            encoder.setBuffer(write, offset: 0, index: 1)
            encoder.setBytes(&step, length: MemoryLayout<UInt32>.stride, index: 2)
            encoder.dispatchThreads(
                MTLSize(width: input.count, height: 1, depth: 1),
                threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
            )
            encoder.endEncoding()
            swap(&read, &write)
            stride <<= 1
        }
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error { throw error }
        try check(command.status == .completed, "Metal command did not complete")
        dispatches += 1
        gpuSeconds += command.gpuEndTime - command.gpuStartTime
        let address = read.contents().bindMemory(to: Int32.self, capacity: input.count)
        return Array(UnsafeBufferPointer(start: address, count: input.count))
    }
}

struct Slice {
    let origin: Int
    let basis: [Int]
    let points: [Int]
    var codim: Int { 4 - basis.count }

    init(_ selected: [Int]) throws {
        let offset = selected[0]
        var vectors: [Int] = []
        var span = [0]
        for point in selected {
            let difference = point ^ offset
            if !span.contains(difference) {
                vectors.append(difference)
                span.append(contentsOf: span.map { $0 ^ difference })
            }
        }
        let ordered = span.map { $0 ^ offset }
        try check(ordered.sorted() == selected, "A slice is not affine")
        self.origin = offset
        self.basis = vectors
        self.points = ordered
    }
}

func partition(_ kind: Int) throws -> [Slice] {
    let row = (0..<16).map { point -> Int in
        let x0 = point & 1
        let x1 = (point >> 1) & 1
        let x2 = (point >> 2) & 1
        let x3 = (point >> 3) & 1
        if kind == 0 {
            if x0 == 0 { return 0 }
            if x1 == x2 { return 1 }
            if x3 == 0 { return 2 }
            return 3 + x2
        }
        if x0 == 0 {
            return (x1 == x3) ? 0 : 1 + x2
        }
        if x2 == x3 { return 3 }
        if x1 == 0 { return 4 }
        return 5 + x3
    }
    return try (0...row.max()!).map { index in try Slice(row.indices.filter { row[$0] == index }) }
}

func rank(_ vectors: [Int]) -> Int {
    var pivots: [Int: Int] = [:]
    for value in vectors {
        var vector = value
        while vector != 0 {
            let lead = Int.bitWidth - vector.leadingZeroBitCount - 1
            if let prior = pivots[lead] { vector ^= prior }
            else {
                pivots[lead] = vector
                break
            }
        }
    }
    return pivots.count
}

struct RNG {
    var state: UInt64
    mutating func next() -> UInt64 {
        state &+= 0x9e3779b97f4a7c15
        var output = state
        output = (output ^ (output >> 30)) &* 0xbf58476d1ce4e5b9
        output = (output ^ (output >> 27)) &* 0x94d049bb133111eb
        return output ^ (output >> 31)
    }
}

func upperSqrt(_ number: UInt64, _ scale: UInt64) -> UInt64 {
    let target = number * scale * scale
    var low: UInt64 = 0
    var high = 16 * scale
    while high - low > 1 {
        let mid = (high + low) / 2
        if mid * mid >= target { high = mid }
        else { low = mid }
    }
    return high
}

func feature(
    mode: Int, piece: Int, u: Int, y: Int, dimensions: Int, rng: inout RNG
) -> Int32 {
    let quadratic = ((u & 3) == 3 ? 1 : 0)
    switch mode {
    case 0:
        // Distinct single data characters, with nonconstant slice phases.
        return sign(parity(piece & y) ^ quadratic)
    case 1:
        // Identical nonaffine data function in every slice.
        return sign((y & 3) == 3 ? 1 : 0)
    case 2:
        // Arbitrary seeded truth tables in each intrinsic domain.
        return sign(Int(rng.next() & 1))
    default:
        // Vary local addressing and the data Fourier bands together.
        let first = parity((piece + u) & y)
        let perturb = (piece & 1 == 1 && y & 3 == 3) ? 1 : 0
        let third = (dimensions == 3 && u & 1 != 0 && y & 4 != 0) ? 1 : 0
        return sign(first ^ quadratic ^ perturb ^ third)
    }
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R41SliceLiftMetal output.json")
    let metal = try MetalWalsh()
    let seed: UInt64 = 0xA280410AB1234567
    var rng = RNG(state: seed)
    var exampleRows: [[String: Any]] = []
    var summaryRows: [[String: Any]] = []
    var count = 0
    var totalGPUCoefficients = 0
    var exactChecks = 0
    var comparedCollisionCases = 0
    var belowNaive = 0

    for partitionIndex in 0..<2 {
        let slices = try partition(partitionIndex)
        try check(slices.flatMap { $0.points }.sorted() == Array(0..<16), "Invalid partition")
        let normalRank = rank(slices.flatMap { slice in
            (0..<16).filter { xi in slice.basis.allSatisfy { parity($0 & xi) == 0 } }
        })

        for d in 2...3 {
            let addressMask = 15
            let domainSize = 1 << (4 + d)
            let dataSize = 1 << d
            var gpuInputs: [Int32] = []
            var cpuWalshRows: [Int32] = []
            var arbitraryCollisionCount = 0

            for trial in 0..<67 {
                let mode = (trial == 0 && d == 3 ? 0 : trial == 1 ? 1 : trial == 2 ? 3 : 2)
                var aggregate = [Int32](repeating: 0, count: domainSize)
                var aggregateTruth = [Int32](repeating: 0, count: domainSize)
                var multiplicities = [Int](repeating: 0, count: domainSize)
                var totalPre = 0
                var sumRank = 0
                var localSparsities: [Int] = []
                var liftedSparsities: [Int] = []
                var weights: [Int] = []

                for (index, slice) in slices.enumerated() {
                    let uSize = slice.points.count
                    let local = (0..<(uSize * dataSize)).map { point in
                        feature(mode: mode, piece: index,
                                u: point % uSize, y: point / uSize, dimensions: d, rng: &rng)
                    }
                    let localWalsh = directWalsh(local)
                    let localSupport = localWalsh.indices.filter { localWalsh[$0] != 0 }
                    let s = localSupport.count
                    var cut = [Int32](repeating: 0, count: domainSize)
                    for y in 0..<dataSize {
                        for (u, x) in slice.points.enumerated() {
                            let point = x + 16 * y
                            cut[point] = local[u + uSize * y]
                            aggregateTruth[point] = cut[point]
                        }
                    }
                    let directLift = directWalsh(cut)
                    let lifted = (0..<domainSize).map { freq -> Int32 in
                        let xi = freq & addressMask
                        let eta = freq >> 4
                        let intrinsic = slice.basis.enumerated().reduce(0) { value, entry in
                            value | (parity(entry.element & xi) << entry.offset)
                        }
                        return sign(parity(xi & slice.origin)) * localWalsh[intrinsic + eta * uSize]
                    }
                    try check(directLift == lifted, "Lift coefficient mismatch")
                    let liftSparsity = lifted.filter { $0 != 0 }.count
                    let weight = 1 << slice.codim
                    try check(liftSparsity == weight * s, "Lift sparsity mismatch")
                    weights.append(weight)
                    localSparsities.append(s)
                    liftedSparsities.append(liftSparsity)
                    totalPre += liftSparsity
                    sumRank += rank(localSupport)
                    for freq in lifted.indices {
                        aggregate[freq] += lifted[freq]
                        if lifted[freq] != 0 { multiplicities[freq] += 1 }
                    }
                    gpuInputs.append(contentsOf: cut)
                    cpuWalshRows.append(contentsOf: directLift)
                    exactChecks += domainSize
                }

                try check(weights.reduce(0) { $0 + 16 / $1 } == 16, "Kraft identity failed")
                let scales = zip(weights, localSparsities).map { $0 * $0 * $1 }
                let equality = Set(scales).count == 1
                let fixedScale: UInt64 = 1 << 20
                // An upper integer enclosure of sum sqrt(s_i); equality is handled algebraically.
                let sqrtSumUpper = localSparsities.reduce(UInt64(0)) { value, s in
                    value + upperSqrt(UInt64(s), fixedScale)
                }
                if !equality {
                    try check(sqrtSumUpper * sqrtSumUpper <= UInt64(totalPre) * fixedScale * fixedScale,
                              "Weighted Cauchy bound failed")
                }

                let directFull = directWalsh(aggregateTruth)
                try check(directFull == aggregate, "Spectral assembly mismatch")
                let support = directFull.indices.filter { directFull[$0] != 0 }
                try check(rank(support) <= normalRank + sumRank, "Global linear rank bound failed")
                let union = multiplicities.filter { $0 > 0 }.count
                let unique = multiplicities.filter { $0 == 1 }.count
                let M = multiplicities.max()!
                let canceled = union - support.count
                let repeats = totalPre - union
                try check(totalPre - repeats - canceled == support.count && support.count >= unique
                          && support.count >= 2 * union - totalPre && union * M >= totalPre,
                          "Collision counting inequality failed")
                let amplitude = Int(directFull.map(abs).max()!)
                try check(support.count * amplitude * amplitude >= domainSize * domainSize,
                          "Postcollision Parseval lower bound failed")
                if mode == 0 { try check(totalPre == support.count && M == 1, "Disjoint bands collided") }
                if canceled > 0 { comparedCollisionCases += 1 }
                if canceled > 0 && mode == 2 { arbitraryCollisionCount += 1 }
                if support.count < totalPre { belowNaive += 1 }
                gpuInputs.append(contentsOf: aggregateTruth)
                cpuWalshRows.append(contentsOf: directFull)
                count += 1

                if mode != 2 || trial == 3 {
                    exampleRows.append([
                        "partition": partitionIndex, "data_variables": d, "mode": mode,
                        "intrinsic_sparsities": localSparsities, "slice_codimensions": slices.map { $0.codim },
                        "lifted_sparsities": liftedSparsities, "precollision_count": totalPre,
                        "precollision_union": union, "postcollision_count": support.count,
                        "duplicate_occurrences": repeats, "complete_cancellations": canceled,
                        "unique_frequencies": unique, "maximum_frequency_multiplicity": M,
                        "global_spectral_linear_dimension": rank(support),
                        "sum_intrinsic_spectral_linear_dimensions": sumRank,
                        "normal_rank": normalRank, "cauchy_equality": equality,
                        "sum_sqrt_intrinsic_sparsities_upper": Double(sqrtSumUpper) / Double(fixedScale),
                    ])
                }
            }
            let gpuSpectrum = try metal.transform(gpuInputs, rowSize: domainSize)
            try check(cpuWalshRows == gpuSpectrum, "CPU/Metal spectra differ")
            totalGPUCoefficients += gpuSpectrum.count
            summaryRows.append([
                "partition": partitionIndex, "data_variables": d, "trials": 67,
                "address_slices": slices.map { $0.points.sorted() },
                "normal_rank": normalRank, "gpu_coefficients": gpuSpectrum.count,
                "arbitrary_truth_cases_with_complete_cancellation": arbitraryCollisionCount,
            ])
        }
    }
    let result: [String: Any] = [
        "schema": "a28-r041-slice-lift-metal-v1",
        "seed_hex": String(seed, radix: 16),
        "address_variables": 4,
        "trials": count,
        "trials_with_complete_cancellation": comparedCollisionCases,
        "trials_with_postcollision_count_below_lift_sum": belowNaive,
        "per_lift_identity_coefficient_checks": exactChecks,
        "cpu_gpu_coefficients_checked": totalGPUCoefficients,
        "cpu_gpu_mismatch_count": 0,
        "cpu_backend": "exact signed Int32 direct Walsh",
        "gpu_backend": "exact Int32 staged Metal FWHT",
        "metal_device": metal.device.name,
        "metal_command_seconds": metal.gpuSeconds,
        "metal_completed_command_buffers": metal.dispatches,
        "cauchy_check": "Kraft weights, exact equality condition and certified 20-bit sqrt upper bounds",
        "examples": exampleRows,
        "series": summaryRows,
        "all_checks_passed": true,
    ]
    let bytes = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try (bytes + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R41 checked \(count) cases, \(totalGPUCoefficients) CPU/Metal coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R41 verifier error: \(error)\n").utf8))
    exit(1)
}
