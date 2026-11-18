import CryptoKit
import Foundation
import Metal

// Each command encoder performs one FWHT stage in a separate output buffer.
let kernelSource = """
#include <metal_stdlib>
using namespace metal;

kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &domainSize [[buffer(2)]],
                       constant uint &stride [[buffer(3)]],
                       uint index [[thread_position_in_grid]]) {
    uint point = index % domainSize;
    int own = input[index];
    int other = input[index ^ stride];
    output[index] = (point & stride) != 0 ? other - own : own + other;
}
"""

struct CheckFailure: Error {
    let reason: String
}

func check(_ condition: Bool, _ reason: String) throws {
    guard condition else { throw CheckFailure(reason: reason) }
}

struct SignedMask: Decodable {
    let points: [Int]
    let solutions: [UInt32]
}

struct PairWitness: Decodable {
    struct Solution: Decodable {
        let f_hex: String
        let f_plus_h_hex: String
    }

    let solution: Solution
}

func hash(_ data: Data) -> String {
    SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func hashSpectrum(_ coefficients: [Int32]) -> String {
    let littleEndian = coefficients.map { $0.littleEndian }
    return littleEndian.withUnsafeBytes { hash(Data($0)) }
}

func parity(_ word: Int) -> Int {
    word.nonzeroBitCount & 1
}

func sign(_ bit: Int) -> Int32 {
    bit == 0 ? 1 : -1
}

func bitsOfHexWitness(_ hex: String) throws -> [Int] {
    try check(hex.hasPrefix("0x") && hex.count <= 66, "Expected at most 256 truth-table bits")
    let nibbles = try hex.dropFirst(2).reversed().map { digit -> Int in
        guard let value = Int(String(digit), radix: 16) else {
            throw CheckFailure(reason: "Invalid truth-table hex digit")
        }
        return value
    }
    return (0..<256).map { point in
        let nibble = point / 4
        return nibble < nibbles.count ? (nibbles[nibble] >> (point % 4)) & 1 : 0
    }
}

func cpuWalsh(_ points: [Int], _ signs: [Int32], _ domainSize: Int) -> [Int32] {
    // Direct summation is independent of the device butterfly algorithm.
    return (0..<domainSize).map { frequency in
        zip(points, signs).reduce(into: Int32(0)) { sum, entry in
            sum += entry.1 * sign(parity(frequency & entry.0))
        }
    }
}

func gpuWalsh(_ input: [Int32], _ domainSize: Int) throws -> ([Int32], [String: Any]) {
    guard let device = MTLCreateSystemDefaultDevice() else {
        throw CheckFailure(reason: "No Metal device")
    }
    let library = try device.makeLibrary(source: kernelSource, options: nil)
    guard let function = library.makeFunction(name: "fwht_stage"),
          let queue = device.makeCommandQueue(),
          let command = queue.makeCommandBuffer()
    else {
        throw CheckFailure(reason: "Unable to prepare the Metal command")
    }
    let pipeline = try device.makeComputePipelineState(function: function)

    let byteCount = input.count * MemoryLayout<Int32>.stride
    let inputBuffer = input.withUnsafeBytes {
        device.makeBuffer(bytes: $0.baseAddress!, length: byteCount, options: .storageModeShared)
    }
    guard var read = inputBuffer,
          var write = device.makeBuffer(length: byteCount, options: .storageModeShared)
    else {
        throw CheckFailure(reason: "Unable to allocate Metal buffers")
    }
    let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
    var stride = 1
    var sizeWord = UInt32(domainSize)
    var stages = 0
    while stride < domainSize {
        guard let encoder = command.makeComputeCommandEncoder() else {
            throw CheckFailure(reason: "Unable to make a Metal encoder")
        }
        var strideWord = UInt32(stride)
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(read, offset: 0, index: 0)
        encoder.setBuffer(write, offset: 0, index: 1)
        encoder.setBytes(&sizeWord, length: MemoryLayout<UInt32>.stride, index: 2)
        encoder.setBytes(&strideWord, length: MemoryLayout<UInt32>.stride, index: 3)
        encoder.dispatchThreads(
            MTLSize(width: input.count, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
        )
        encoder.endEncoding()
        swap(&read, &write)
        stride <<= 1
        stages += 1
    }
    command.commit()
    command.waitUntilCompleted()
    if let error = command.error { throw error }
    try check(command.status == .completed, "Metal command did not complete")
    let address = read.contents().bindMemory(to: Int32.self, capacity: input.count)
    let spectrum = Array(UnsafeBufferPointer(start: address, count: input.count))
    return (spectrum, [
        "backend": "Metal Int32 staged FWHT",
        "device": device.name,
        "metal4_supported": device.supportsFamily(.metal4),
        "stages": stages,
        "threads_per_stage": input.count,
        "command_buffer_completed": true,
        "gpu_command_seconds": command.gpuEndTime - command.gpuStartTime,
        "shared_buffer_bytes_total": 2 * byteCount,
    ])
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 4, "usage: R34S32AffineMetal signedmask.json mm_completion.json output.json")
    let signedData = try Data(contentsOf: URL(fileURLWithPath: args[1]))
    let witnessData = try Data(contentsOf: URL(fileURLWithPath: args[2]))
    let signed = try JSONDecoder().decode(SignedMask.self, from: signedData)
    let witness = try JSONDecoder().decode(PairWitness.self, from: witnessData)
    let domainSize = 128
    let points = signed.points
    try check(points.count == 32 && points == Array(Set(points)).sorted(), "Expected 32 sorted points")
    try check(points.allSatisfy { (0..<domainSize).contains($0) }, "Point outside F2^7")
    try check(signed.solutions.count == 1, "Expected the single archived normalized R32 mask")
    let referenceMask = signed.solutions[0]

    let f = try bitsOfHexWitness(witness.solution.f_hex)
    let flipped = try bitsOfHexWitness(witness.solution.f_plus_h_hex)
    try check((0..<256).filter { f[$0] != flipped[$0] } == points, "R33 difference support mismatch")
    try check(points.enumerated().allSatisfy {
        f[$0.element] == Int((referenceMask >> $0.offset) & 1)
    }, "R33 restriction disagrees with R32 reference signs")

    let beginning = ProcessInfo.processInfo.systemUptime
    let referenceSigns = points.indices.map { sign(Int((referenceMask >> $0) & 1)) }
    let referenceWalsh = cpuWalsh(points, referenceSigns, domainSize)
    let referenceSupport = referenceWalsh.indices.filter { referenceWalsh[$0] != 0 }
    try check(referenceSupport.count == 16 && referenceSupport.allSatisfy {
        abs(referenceWalsh[$0]) == 16
    }, "Reference spectrum must have exactly 16 nonzero values of absolute value 16")

    // Row 128*c+a is the signed mask multiplied by (-1)^(c + a dot x).
    let phaseCount = 2 * domainSize
    var cpu: [Int32] = []
    cpu.reserveCapacity(phaseCount * domainSize)
    var gpuInput = [Int32](repeating: 0, count: phaseCount * domainSize)
    var realizedMasks = Set<UInt32>()
    var translatedSupports = Set<[Int]>()
    var walshShiftMismatches = 0
    var translatedSupportMismatches = 0
    for phase in 0..<phaseCount {
        let slope = phase % domainSize
        let constant = phase / domainSize
        let rowStart = phase * domainSize
        let signs = points.enumerated().map { entry in
            referenceSigns[entry.offset] * sign(constant ^ parity(slope & entry.element))
        }
        var phaseMask: UInt32 = 0
        for (index, point) in points.enumerated() {
            gpuInput[rowStart + point] = signs[index]
            if signs[index] < 0 { phaseMask |= UInt32(1) << index }
        }
        realizedMasks.insert(phaseMask)

        let computed = cpuWalsh(points, signs, domainSize)
        cpu.append(contentsOf: computed)
        for frequency in 0..<domainSize {
            if computed[frequency] != sign(constant) * referenceWalsh[frequency ^ slope] {
                walshShiftMismatches += 1
            }
        }
        let computedSupport = computed.indices.filter { computed[$0] != 0 }
        translatedSupports.insert(computedSupport)
        if computedSupport != referenceSupport.map({ $0 ^ slope }).sorted() {
            translatedSupportMismatches += 1
        }
    }
    let (gpu, deviceRecord) = try gpuWalsh(gpuInput, domainSize)
    let gpuCpuMismatches = zip(gpu, cpu).reduce(into: 0) { count, pair in
        if pair.0 != pair.1 { count += 1 }
    }
    let stabilizer = (0..<domainSize).filter { shift in
        referenceSupport.map { $0 ^ shift }.sorted() == referenceSupport
    }
    let passed = realizedMasks.count == phaseCount && gpuCpuMismatches == 0
        && walshShiftMismatches == 0 && translatedSupportMismatches == 0
        && translatedSupports.count * stabilizer.count == domainSize
    let record: [String: Any] = [
        "schema": "a28-r034-s32-affine-metal-v1",
        "r032_signed_input_sha256": hash(signedData),
        "r033_witness_input_sha256": hash(witnessData),
        "r033_witness_sign_and_support_match": true,
        "normalization": "unnormalized real Walsh transform on F2^7, little-endian Int32",
        "domain_size": domainSize,
        "physical_support": points,
        "reference_negative_mask": referenceMask,
        "reference_walsh": referenceWalsh,
        "reference_nonzero_support": referenceSupport,
        "frequency_translation_stabilizer": stabilizer,
        "phase_count": phaseCount,
        "distinct_signed_masks": realizedMasks.count,
        "distinct_frequency_supports": translatedSupports.count,
        "coefficients_checked": cpu.count,
        "cpu_backend": "direct exact Int32 dot-product Walsh summation",
        "gpu": deviceRecord,
        "cpu_walsh_sha256": hashSpectrum(cpu),
        "gpu_walsh_sha256": hashSpectrum(gpu),
        "gpu_cpu_coefficient_mismatches": gpuCpuMismatches,
        "affine_coefficient_formula_mismatches": walshShiftMismatches,
        "frequency_support_translation_mismatches": translatedSupportMismatches,
        "elapsed_seconds_excluding_io_and_compile": ProcessInfo.processInfo.systemUptime - beginning,
        "all_checks_passed": passed,
    ]
    let output = try JSONSerialization.data(withJSONObject: record, options: [.prettyPrinted, .sortedKeys])
    try (output + Data([10])).write(to: URL(fileURLWithPath: args[3]))
    try check(passed, "Walsh cross-check failed")
    print("Verified \(phaseCount) affine phases and \(cpu.count) Walsh coefficients on CPU and Metal")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R34 verifier error: \(error)\n").utf8))
    exit(1)
}
