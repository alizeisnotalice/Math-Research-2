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

func pairTally(_ mask: Int) -> Int {
    ((mask & 255) ^ (mask >> 8)).nonzeroBitCount
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R37PhaseModesMetal output.json")
    let n = 16
    let phaseTables = 1 << n
    let characters = (0..<n).map { frequency in
        (0..<n).reduce(into: 0) { word, point in
            if parity(frequency & point) == 1 { word |= 1 << point }
        }
    }
    let metal = try MetalWalsh()
    var rows: [Int32] = []
    rows.reserveCapacity(n * phaseTables)
    for truth in 0..<phaseTables {
        for point in 0..<n {
            rows.append(sign((truth >> point) & 1))
        }
    }
    let gpu = try metal.transform(rows, rowSize: n)
    var rawHistogram: [String: Int] = [:]
    var q16Modes: [String: Int] = [:]
    var q8Allowed = 0
    var oneMinority = 0
    var twoMinorities = 0
    var nMod4Failures = 0
    var nMismatches = 0
    var signModeFailures = 0
    var distinctHighLabelSupports = Set<[Int]>()
    var h = SHA256()
    for truth in 0..<phaseTables {
        let w = characters.map { n - 2 * (truth ^ $0).nonzeroBitCount }
        for frequency in 0..<n {
            let value = w[frequency]
            if Int(gpu[truth * n + frequency]) != value { nMismatches += 1 }
            if (value - w[0]) % 4 != 0 { nMod4Failures += 1 }
            rawHistogram[String(value), default: 0] += 1
        }
        var littleEndian = w.map { Int16($0).littleEndian }
        littleEndian.withUnsafeMutableBytes { h.update(bufferPointer: UnsafeRawBufferPointer($0)) }

        // Distinguish two halves of the dual parameter cube using a=e3.
        let low = 1..<8
        let high = 8..<16
        let q16AtomBounds = low.allSatisfy { abs(w[$0]) <= 2 }
            && high.allSatisfy { abs(w[$0]) <= 4 }
        if q16AtomBounds {
            q16Modes[String(pairTally(truth)), default: 0] += 1
        }
        let q8AtomBoundsAndParity = low.allSatisfy { abs(w[$0]) == 2 }
            && high.allSatisfy { abs(w[$0]) <= 4 && w[$0] % 4 == 0 }
        if q8AtomBoundsAndParity { q8Allowed += 1 }

        let minorityWeight = min(truth.nonzeroBitCount, n - truth.nonzeroBitCount)
        if abs(w[0]) == 14 {
            oneMinority += 1
            if minorityWeight != 1 || !w.dropFirst().allSatisfy({ abs($0) == 2 }) {
                signModeFailures += 1
            }
        } else if abs(w[0]) == 12 {
            twoMinorities += 1
            if minorityWeight != 2 || !w.dropFirst().allSatisfy({ $0 == 0 || abs($0) == 4 }) {
                signModeFailures += 1
            }
            let minorityBits = (truth.nonzeroBitCount == 2 ? truth : truth ^ 65535)
            let points = (0..<16).filter { (minorityBits >> $0) & 1 != 0 }
            let labels = (1..<16).filter { w[$0] != 0 }
            if points.count != 2 || labels != (1..<16).filter({
                parity($0 & (points[0] ^ points[1])) == 0
            }) {
                signModeFailures += 1
            }
            distinctHighLabelSupports.insert(labels)
        }
    }

    try check(q16Modes == ["0": 2, "1": 32], "Unexpected q16 local mode split")
    try check(q8Allowed == 0, "Unexpected physically compatible q8 mode")
    try check(oneMinority == 32 && twoMinorities == 240 && distinctHighLabelSupports.count == 15,
              "Incorrect minority mode census")
    try check(nMismatches == 0 && nMod4Failures == 0 && signModeFailures == 0,
              "Coefficient or combinatorial checks failed")

    let result: [String: Any] = [
        "schema": "a28-r037-16-sign-walsh-mode-census-v1",
        "number_of_parameter_sign_tables": phaseTables,
        "coefficients_checked": n * phaseTables,
        "coefficient_cpu_mismatches": nMismatches,
        "all_table_coeffs_same_residue_mod4_violations": nMod4Failures,
        "minority_mode_support_violations": signModeFailures,
        "raw_walsh_histogram": rawHistogram,
        "raw_coefficients_sha256_cpu": h.finalize().map { String(format: "%02x", $0) }.joined(),
        "q16_atom_bound_modes_by_num_opposite_pairs": q16Modes,
        "q8_atom_bound_and_parity_modes": q8Allowed,
        "one_minority_sign_modes": oneMinority,
        "two_minority_sign_modes": twoMinorities,
        "two_minority_nonzero_supports": distinctHighLabelSupports.count,
        "device": metal.device.name,
        "gpu_backend": "Metal Int32 staged FWHT",
        "metal4_supported": metal.device.supportsFamily(.metal4),
        "gpu_completed_command_buffers": metal.dispatches,
        "gpu_command_seconds": metal.gpuSeconds,
        "cpu_backend": "exact direct character masks and popcounts",
        "all_checks_passed": true,
    ]
    let bytes = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try (bytes + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R37 verified all \(phaseTables) sign tables and \(n * phaseTables) Walsh coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R37 verifier error: \(error)\n").utf8))
    exit(1)
}
