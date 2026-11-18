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

struct MMSharpInput: Decodable {
    struct Extension: Decodable {
        let representative: String
        let family: [String]
        let exterior_weight: Int
        let e4_frequency_support_moves: Bool
        let total_union_weight: Int
    }

    let base_truth_table_hex: String
    let e4_physical_support_hex: String
    let e4_generators_hex: [String]
    let e4_family_hex: [String]
    let sharp_extensions: [Extension]
}

func parse64(_ text: String) throws -> UInt64 {
    guard text.count == 16, let value = UInt64(text, radix: 16) else {
        throw CheckFailure(reason: "Invalid hex truth table")
    }
    return value
}

func support(_ value: UInt64) -> [Int] {
    (0..<64).filter { value & (UInt64(1) << $0) != 0 }
}

func truthSigns(_ value: UInt64) -> [Int32] {
    (0..<64).map { sign(Int((value >> $0) & 1)) }
}

func partialSpectrum(_ truth: UInt64, _ block: UInt64) -> [Int] {
    let signTable = truthSigns(truth)
    let cut = signTable.enumerated().map { (block >> $0.offset) & 1 != 0 ? $0.element : 0 }
    return directWalsh(cut).enumerated().compactMap { $0.element != 0 ? $0.offset : nil }
}

func span(_ basis: [UInt64]) -> [UInt64] {
    var points: [UInt64] = [0]
    for vector in basis { points.append(contentsOf: points.map { $0 ^ vector }) }
    return points.sorted()
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 3, "usage: R36E4Metal enumeration.json output.json")
    let data = try Data(contentsOf: URL(fileURLWithPath: args[1]))
    let witness = try JSONDecoder().decode(MMSharpInput.self, from: data)
    let base = try parse64(witness.base_truth_table_hex)
    let block = try parse64(witness.e4_physical_support_hex)
    let e4 = try witness.e4_family_hex.map(parse64)
    let generators = try witness.e4_generators_hex.map(parse64)
    try check(e4 == span(generators) && Set(e4).count == 8, "Invalid E4 generators")
    try check(e4.reduce(0, |) == block && block.nonzeroBitCount == 16, "E4 support mismatch")
    try check(support(base) == (0..<64).filter {
        parity(($0 & 7) & ($0 >> 3)) == 1
    }, "Incorrect quadratic base")
    try check(witness.sharp_extensions.count == 90, "Wrong number of sharp examples")

    let canonicalPartialSpectrum = partialSpectrum(base, block)
    let characterMasks = (0..<64).map { frequency in
        (0..<64).reduce(into: UInt64(0)) { value, point in
            if parity(frequency & point) != 0 { value |= UInt64(1) << point }
        }
    }

    var phaseRows: [UInt64] = []
    var moving = 0
    var stationary = 0
    var distinctSupports = Set<UInt64>()
    var coveredCosets = Set<UInt64>()
    var examples: [String: Any] = [:]
    for row in witness.sharp_extensions {
        let representative = try parse64(row.representative)
        let family = try row.family.map(parse64)
        try check(family == span(generators + [representative]) && Set(family).count == 16,
                  "Invalid four-dimensional family")
        let coset = e4.map { $0 ^ representative }.sorted()
        try check(representative == coset[0], "Noncanonical coset")
        try check(coveredCosets.insert(representative).inserted, "Repeated extension modulo E4")
        let union = family.reduce(0, |)
        let outside = representative & ~block
        try check(outside.nonzeroBitCount == 8 && union.nonzeroBitCount == 24,
                  "Sharp extension does not use 3q points")
        try check(row.exterior_weight == 8 && row.total_union_weight == 24, "Recorded weights differ")
        distinctSupports.insert(union)
        let translated = partialSpectrum(base ^ representative, block)
        let moves = (translated != canonicalPartialSpectrum)
        let requiredIntersection = moves ? 0 : canonicalPartialSpectrum.count
        try check(Set(translated).intersection(canonicalPartialSpectrum).count == requiredIntersection,
                  "Partial spectra were unequal but not disjoint")
        try check(moves == row.e4_frequency_support_moves, "Partial spectrum motion mismatch")
        if moves { moving += 1 } else { stationary += 1 }
        let kind = moves ? "moving" : "stationary"
        if examples[kind] == nil {
            examples[kind] = [
                "representative": row.representative,
                "union_support": support(union),
                "outside_support": support(outside),
                "changed_e4_partial_support": translated,
            ]
        }
        phaseRows.append(contentsOf: family)
    }

    let metal = try MetalWalsh()
    let signs = phaseRows.flatMap { truthSigns(base ^ $0) }
    let device = try metal.transform(signs, rowSize: 64)
    var mismatches = 0
    var absViolations = 0
    for (row, mask) in phaseRows.enumerated() {
        let truth = base ^ mask
        for frequency in 0..<64 {
            let exact = Int32(64 - 2 * (truth ^ characterMasks[frequency]).nonzeroBitCount)
            if exact != device[row * 64 + frequency] { mismatches += 1 }
            if abs(exact) != 8 { absViolations += 1 }
        }
    }
    try check(mismatches == 0 && absViolations == 0, "CPU/GPU Walsh verification failed")
    let result: [String: Any] = [
        "schema": "a28-r036-e4-sharp-metal-v1",
        "input_sha256": hash(data),
        "device": metal.device.name,
        "metal4_supported": metal.device.supportsFamily(.metal4),
        "gpu_command_seconds": metal.gpuSeconds,
        "completed_command_buffers": metal.dispatches,
        "cpu_backend": "exact Int32 direct popcount Walsh",
        "gpu_backend": "Metal Int32 staged FWHT",
        "phase_family_count": moving + stationary,
        "stationary_family_count": stationary,
        "moving_family_count": moving,
        "phase_rows": phaseRows.count,
        "distinct_phase_masks": Set(phaseRows).count,
        "coefficients_checked": device.count,
        "gpu_cpu_mismatches": mismatches,
        "cpu_nonbent_coefficients": absViolations,
        "distinct_union_supports": distinctSupports.count,
        "minimum_union_size_of_these_extensions": 24,
        "e4_reference_partial_support": canonicalPartialSpectrum,
        "representative_examples": examples,
        "all_checks_passed": true,
    ]
    let output = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try (output + Data([10])).write(to: URL(fileURLWithPath: args[2]))
    print("R36 CPU+Metal verified \(phaseRows.count) phases in \(moving + stationary) sharp E4 extensions.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R36 verifier error: \(error)\n").utf8))
    exit(1)
}
