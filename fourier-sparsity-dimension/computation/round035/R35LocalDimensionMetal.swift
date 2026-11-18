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

func bitsOfHex(_ hex: String) throws -> [Int32] {
    try check(hex.hasPrefix("0x") && hex.count <= 66, "Invalid 256-bit truth table")
    let reversed = try hex.dropFirst(2).reversed().map { digit -> Int in
        guard let value = Int(String(digit), radix: 16) else {
            throw CheckFailure(reason: "Invalid hex digit")
        }
        return value
    }
    return (0..<256).map { point in
        let nibble = point / 4
        return nibble < reversed.count ? sign((reversed[nibble] >> (point % 4)) & 1) : 1
    }
}

func directWalsh(_ table: [Int32]) -> [Int32] {
    (0..<table.count).map { frequency in
        table.indices.reduce(into: Int32(0)) { sum, point in
            sum += table[point] * sign(parity(frequency & point))
        }
    }
}

func xorSpan(_ basis: [Int]) -> [Int] {
    var points = [0]
    for vector in basis { points.append(contentsOf: points.map { $0 ^ vector }) }
    return points.sorted()
}

func affineDimension(_ points: [Int]) -> Int {
    var basis: [Int: Int] = [:]
    for point in points {
        var vector = point ^ points[0]
        while vector != 0 {
            let leading = Int.bitWidth - vector.leadingZeroBitCount - 1
            if let existing = basis[leading] {
                vector ^= existing
            } else {
                basis[leading] = vector
                break
            }
        }
    }
    return basis.count
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

func readJSON(_ args: [String], _ index: Int) throws -> (Data, [String: Any]) {
    let bytes = try Data(contentsOf: URL(fileURLWithPath: args[index]))
    guard let object = try JSONSerialization.jsonObject(with: bytes) as? [String: Any] else {
        throw CheckFailure(reason: "Expected a JSON object")
    }
    return (bytes, object)
}

func subsetMask(_ points: [Int], _ condition: (Int) -> Bool) -> Int {
    points.enumerated().reduce(into: 0) { mask, entry in
        if condition(entry.element) { mask |= 1 << entry.offset }
    }
}

func spectrumOnBlock(_ base: [Int32], _ points: [Int]) -> [Int] {
    let kept = Set(points)
    return directWalsh(base.enumerated().map { kept.contains($0.offset) ? $0.element : 0 })
        .enumerated().compactMap { $0.element != 0 ? $0.offset : nil }
}

func exhaustiveSix(
    name: String,
    block: [Int],
    base: [Int32],
    expectedBasis: [Int],
    metal: MetalWalsh
) throws -> [String: Any] {
    try check(block.count == 16, "Exhaustive audit requires 16-point block")
    let start = ProcessInfo.processInfo.systemUptime
    let n = base.count
    let baseWalsh = directWalsh(base)
    let expected = xorSpan(expectedBasis)
    let localTruthMask = subsetMask(block) { base[$0] == -1 }
    let characterMasks = (0..<n).map { frequency in
        subsetMask(block) { parity(frequency & $0) == 1 } ^ localTruthMask
    }
    let oldPartial = characterMasks.map { Int32(block.count - 2 * $0.nonzeroBitCount) }
    let phaseCount = 1 << block.count
    let batchSize = 4096
    var compatible: [Int] = []
    var mismatches = 0

    for first in stride(from: 0, to: phaseCount, by: batchSize) {
        var deviceInput: [Int32] = []
        deviceInput.reserveCapacity(batchSize * n)
        for phase in first..<(first + batchSize) {
            var truth = base
            for (i, point) in block.enumerated() where phase & (1 << i) != 0 {
                truth[point] = -truth[point]
            }
            deviceInput.append(contentsOf: truth)
        }
        let deviceOutput = try metal.transform(deviceInput, rowSize: n)
        for phase in first..<(first + batchSize) {
            let rowStart = (phase - first) * n
            var bent = true
            for frequency in 0..<n {
                let newPartial = Int32(block.count - 2 * (phase ^ characterMasks[frequency]).nonzeroBitCount)
                let expectedCoefficient = baseWalsh[frequency] - oldPartial[frequency] + newPartial
                if expectedCoefficient != deviceOutput[rowStart + frequency] { mismatches += 1 }
                if abs(expectedCoefficient) != 8 { bent = false }
            }
            if bent { compatible.append(phase) }
        }
    }

    let full = phaseCount - 1
    let compatibleSet = Set(compatible)
    let complementary = compatible.filter { compatibleSet.contains($0 ^ full) }
    try check(mismatches == 0, "CPU/Metal spectra differ")
    try check(complementary == expected, "Complement-compatible masks differ from anticipated linear family")
    try check(compatible == expected, "Additional compatible single flips")
    return [
        "example": name,
        "physical_support": block,
        "affine_dimension_of_support": affineDimension(block),
        "partial_spectrum_support": spectrumOnBlock(base, block),
        "phase_count": phaseCount,
        "coefficients_checked": n * phaseCount,
        "coefficient_mismatches": mismatches,
        "compatible_masks": compatible,
        "complement_compatible_masks": complementary,
        "expected_linear_masks": expected,
        "expected_generators": expectedBasis,
        "exact_local_linear_dimension": expectedBasis.count,
        "elapsed_seconds": ProcessInfo.processInfo.systemUptime - start,
    ]
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 5, "usage: R35LocalDimensionMetal r32.json r33.json r34.json output.json")
    let (r32Bytes, r32) = try readJSON(args, 1)
    let (r33Bytes, r33) = try readJSON(args, 2)
    let (r34Bytes, r34) = try readJSON(args, 3)
    let supports = r32["points"] as! [Int]
    let solution = r33["solution"] as! [String: Any]
    let base8 = try bitsOfHex(solution["f_hex"] as! String)
    let flipped8 = try bitsOfHex(solution["f_plus_h_hex"] as! String)
    let difference = base8.indices.filter { base8[$0] != flipped8[$0] }
    let r34Rows = r34["dimensions"] as! [String: [String: Any]]
    let r34Physical = r34Rows["4"]!["S32"] as! [String: Any]
    let count4 = r34Physical["contained_flats_count"] as! Int
    try check(r34["signed_input_sha256"] as? String == hash(r32Bytes), "R34 checked a different physical input")
    try check(r34Rows["4"]!["directions"] as? Int == 11811, "Incomplete R34 dimension-four audit")
    try check(r34Rows["4"]!["affine_flats"] as? Int == 94488, "Incomplete R34 affine four-flat coverage")
    try check(count4 == 0 && difference == supports && difference.count == 32, "Eight-dimensional support check failed")

    let metal = try MetalWalsh()
    let cpu8 = directWalsh(base8) + directWalsh(flipped8)
    let gpu8 = try metal.transform(base8 + flipped8, rowSize: 256)
    try check(cpu8 == gpu8 && cpu8.allSatisfy { abs($0) == 16 }, "R33 pair no longer verified as bent")

    let domain = 0..<64
    let base6 = domain.map { point in sign(parity((point & 7) & (point >> 3))) }
    let a = domain.filter { $0 & 7 == 0 }
    let b = domain.filter { $0 & 35 == 1 } // x0=1, x1=0, y2=0
    try check(Set(a).isDisjoint(with: Set(b)), "Physical blocks intersect")
    let partialA = spectrumOnBlock(base6, a)
    let partialB = spectrumOnBlock(base6, b)
    try check(Set(partialA).isDisjoint(with: Set(partialB)), "Partial spectra intersect")
    let nonaffine = (a + b).sorted()
    try check(affineDimension(nonaffine) > 4, "The union of two minima should be nonaffine")
    let aSet = Set(a)
    let bSet = Set(b)
    let nonaffineAudit = try exhaustiveSix(
        name: "union of disjoint nonparallel minimal blocks",
        block: nonaffine, base: base6,
        expectedBasis: [subsetMask(nonaffine) { aSet.contains($0) }, subsetMask(nonaffine) { bSet.contains($0) }],
        metal: metal
    )
    let affine = domain.filter { $0 & 3 == 0 } // x0=x1=0
    try check(affineDimension(affine) == 4, "Expected affine dimension four")
    let affineAudit = try exhaustiveSix(
        name: "affine four-flat",
        block: affine, base: base6,
        expectedBasis: [
            subsetMask(affine) { $0 & 4 == 0 },
            subsetMask(affine) { $0 & 4 != 0 },
            subsetMask(affine) { $0 & 32 != 0 },
        ],
        metal: metal
    )

    let result: [String: Any] = [
        "schema": "a28-r035-local-dimension-metal-v1",
        "r32_input_sha256": hash(r32Bytes),
        "r33_input_sha256": hash(r33Bytes),
        "r34_flat_audit_input_sha256": hash(r34Bytes),
        "eight_variables": [
            "physical_support": supports,
            "reference_coefficients_checked": cpu8.count,
            "reference_pair_bent": true,
            "cpu_metal_mismatch_count": 0,
            "bent_difference_minimum_weight": 16,
            "minimum_weight_support": "affine 4-flat by the minimum-weight theorem for RM(4,8)",
            "contained_four_flats": count4,
            "exact_local_linear_dimension_containing_full_flip": 1,
            "derivation": "Two complementary nonzero masks in a local linear space would have weights summing to 32, forcing both to be 16.",
        ],
        "six_variables": [
            "base": "inner product of low and high 3-bit coordinates",
            "nonaffine_block_A": a,
            "nonaffine_block_B": b,
            "nonaffine_partial_spectrum_A": partialA,
            "nonaffine_partial_spectrum_B": partialB,
            "nonaffine_block": nonaffineAudit,
            "affine_block": affineAudit,
        ],
        "cpu_backend": "exact Int32 direct Walsh summation and independent popcount updates",
        "gpu_backend": "Metal Int32 staged FWHT",
        "device": metal.device.name,
        "metal4_supported": metal.device.supportsFamily(.metal4),
        "completed_command_buffers": metal.dispatches,
        "gpu_command_seconds": metal.gpuSeconds,
        "all_checks_passed": true,
    ]
    let output = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try (output + Data([10])).write(to: URL(fileURLWithPath: args[4]))
    print("R35-CPU-Metal: verified both complete six-variable 16-point phase sets and the R33 base pair")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R35 verifier error: \(error)\n").utf8))
    exit(1)
}
