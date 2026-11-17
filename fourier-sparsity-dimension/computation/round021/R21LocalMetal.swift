import Foundation
import Metal
import CryptoKit

// Independent finite audit for B={3,7,11,12,13,14}.  For m in [-4,4]^6,
// W_m(x)=sum_j m_j (-1)^(B_j dot x).  The CPU oracle below evaluates that
// defining dot product directly; Metal is only an independent staged FWHT.
let b: [Int] = [3, 7, 11, 12, 13, 14]
let n = 16
let coefficientCount = 531_441
let batchLimit = 4096

let kernelSource = #"""
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &stride [[buffer(2)]],
                       uint gid [[thread_position_in_grid]]) {
    uint local = gid & 15u;
    uint base = gid - local;
    uint mate = base + (local ^ stride);
    int left = input[gid];
    int right = input[mate];
    output[gid] = (local & stride) ? (right - left) : (left + right);
}
"""#

func sha256(_ path: String) -> String {
    let data = try! Data(contentsOf: URL(fileURLWithPath: path))
    return SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func parity(_ v: Int) -> Int { v.nonzeroBitCount & 1 }

func decode(_ index: Int) -> [Int] {
    var q = index
    var m = [Int](repeating: 0, count: 6)
    for j in 0..<6 { m[j] = q % 9 - 4; q /= 9 }
    return m
}

func directWalsh(_ m: [Int]) -> [Int32] {
    return (0..<n).map { x in
        var s = 0
        for j in 0..<6 { s += parity(x & b[j]) == 0 ? m[j] : -m[j] }
        return Int32(s)
    }
}

func metalFWHT(device: MTLDevice, pipeline: MTLComputePipelineState, rows: [[Int]]) throws -> [[Int32]] {
    let count = rows.count
    var initial = [Int32](repeating: 0, count: count * n)
    for r in 0..<count {
        for j in 0..<6 { initial[r * n + b[j]] = Int32(rows[r][j]) }
    }
    let length = initial.count * MemoryLayout<Int32>.stride
    let options: MTLResourceOptions = .storageModeShared
    let first = initial.withUnsafeBytes { device.makeBuffer(bytes: $0.baseAddress!, length: length, options: options)! }
    let second = device.makeBuffer(length: length, options: options)!
    var input = first, output = second
    let queue = device.makeCommandQueue()!
    let command = queue.makeCommandBuffer()!
    for h: UInt32 in [1, 2, 4, 8] {
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(input, offset: 0, index: 0)
        encoder.setBuffer(output, offset: 0, index: 1)
        var stride = h
        encoder.setBytes(&stride, length: MemoryLayout<UInt32>.stride, index: 2)
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: count * n, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding()
        swap(&input, &output)
    }
    command.commit()
    command.waitUntilCompleted()
    guard command.status == .completed else { throw NSError(domain: "metal", code: Int(command.status.rawValue)) }
    let result = input.contents().bindMemory(to: Int32.self, capacity: count * n)
    return (0..<count).map { Array(UnsafeBufferPointer(start: result + $0 * n, count: n)) }
}

let arguments = CommandLine.arguments
guard arguments.count == 2 else { fatalError("usage: R21CubeWalshMetal output.json") }
let outputPath = arguments[1]
guard let device = MTLCreateSystemDefaultDevice() else { fatalError("Metal unavailable") }
let library = try! device.makeLibrary(source: kernelSource, options: nil)
let pipeline = try! device.makeComputePipelineState(function: library.makeFunction(name: "fwht_stage")!)
let start = Date()
var signedCoefficientChecks = 0
var gpuCPUMismatches = 0
var invalidCPU = 0
var invalidGPU = 0
var solutions: [[String: Any]] = []
var cursor = 0
while cursor < coefficientCount {
    let end = min(cursor + batchLimit, coefficientCount)
    let rows = (cursor..<end).map(decode)
    let cpu = rows.map(directWalsh)
    let gpu = try! metalFWHT(device: device, pipeline: pipeline, rows: rows)
    for i in 0..<rows.count {
        for x in 0..<n {
            signedCoefficientChecks += 1
            if cpu[i][x] != gpu[i][x] { gpuCPUMismatches += 1 }
            if cpu[i][x] != 0 && abs(cpu[i][x]) != 16 { invalidCPU += 1 }
            if gpu[i][x] != 0 && abs(gpu[i][x]) != 16 { invalidGPU += 1 }
        }
        if cpu[i].allSatisfy({ $0 == 0 || abs($0) == 16 }) {
            solutions.append(["index_base9": cursor + i, "m": rows[i], "Wm": cpu[i]])
        }
    }
    cursor = end
}
let nonzeroSolutions = solutions.filter { ($0["m"] as! [Int]).allSatisfy { abs($0) == 4 } }
let allOddNegativeParity = nonzeroSolutions.allSatisfy {
    ($0["m"] as! [Int]).filter { $0 < 0 }.count % 2 == 1
}
let expected = solutions.count == 33 && solutions.filter { ($0["m"] as! [Int]).allSatisfy { $0 == 0 } }.count == 1 &&
    nonzeroSolutions.count == 32 && allOddNegativeParity
let result: [String: Any] = [
    "schema": "a28-r021-cube-walsh-metal-v1",
    "B": b, "coefficient_vectors_checked": coefficientCount,
    "signed_walsh_coefficients_checked": signedCoefficientChecks,
    "cpu_direct_oracle": "sum_j m_j*(-1)^(x dot B_j), explicit 16x6 dot products",
    "gpu_backend": "Metal int32 global staged ping-pong FWHT, four stages; CPU Int32 direct oracle",
    "batch_limit": batchLimit, "device": device.name, "registry_id": String(device.registryID),
    "gpu_cpu_mismatch_count": gpuCPUMismatches, "invalid_cpu_coefficients": invalidCPU,
    "invalid_gpu_coefficients": invalidGPU, "solution_count": solutions.count,
    "solutions": solutions, "nonzero_solutions_all_have_odd_negative_parity": allOddNegativeParity,
    "expected_33_zero_plus_32_all_abs4_odd_negative_parity": expected,
    "all_checks_passed": gpuCPUMismatches == 0 && invalidCPU == invalidGPU && expected,
    "elapsed_seconds": Date().timeIntervalSince(start)
]
try! JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys]).write(to: URL(fileURLWithPath: outputPath))
if !(result["all_checks_passed"] as! Bool) { exit(10) }
print("R21 passed: \(coefficientCount) vectors; \(signedCoefficientChecks) signed coefficients; \(solutions.count) solutions")
