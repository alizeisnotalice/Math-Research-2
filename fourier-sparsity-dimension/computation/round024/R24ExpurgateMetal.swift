import Foundation
import Metal
import CryptoKit

let n = 64
let affineDimension = 7

let kernelSource = """
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &stride [[buffer(2)]],
                       uint gid [[thread_position_in_grid]]) {
    uint local = gid & 63u;
    uint base = gid - local;
    uint mate = base + (local ^ stride);
    int left = input[gid];
    int right = input[mate];
    output[gid] = (local & stride) ? (right - left) : (left + right);
}
"""

func sha256(_ url: URL) -> String {
    let d = try! Data(contentsOf: url)
    return SHA256.hash(data: d).map { String(format: "%02x", $0) }.joined()
}
func parity(_ x: Int) -> Int { x.nonzeroBitCount & 1 }
func rankWords(_ originals: [UInt64]) -> Int {
    var pivot = [UInt64](repeating: 0, count: 64)
    for original in originals {
        var x = original
        while x != 0 {
            let p = 63 - x.leadingZeroBitCount
            if pivot[p] == 0 { pivot[p] = x; break }
            x ^= pivot[p]
        }
    }
    return pivot.filter { $0 != 0 }.count
}
func truthForCoordinate(_ coordinate: Int) -> UInt64 {
    var out: UInt64 = 0
    for x in 0..<n where ((x >> coordinate) & 1) == 1 { out |= UInt64(1) << UInt64(x) }
    return out
}
func directWalsh(_ phases: [UInt64]) -> [Int32] {
    var out = [Int32](repeating: 0, count: phases.count * n)
    for p in 0..<phases.count {
        let word = phases[p], base = p * n
        for xi in 0..<n {
            var total = 0
            for x in 0..<n {
                let bit = Int((word >> UInt64(x)) & 1) ^ parity(x & xi)
                total += bit == 0 ? 1 : -1
            }
            out[base + xi] = Int32(total)
        }
    }
    return out
}
func metalFWHT(device: MTLDevice, pipeline: MTLComputePipelineState, phases: [UInt64]) throws -> [Int32] {
    var input = [Int32](repeating: 0, count: phases.count * n)
    for p in 0..<phases.count { for x in 0..<n { input[p*n+x] = ((phases[p] >> UInt64(x)) & 1) == 0 ? 1 : -1 } }
    let length = input.count * MemoryLayout<Int32>.stride
    let options: MTLResourceOptions = .storageModeShared
    let first = input.withUnsafeBytes { device.makeBuffer(bytes: $0.baseAddress!, length: length, options: options)! }
    let second = device.makeBuffer(length: length, options: options)!
    var current = first, other = second
    let command = device.makeCommandQueue()!.makeCommandBuffer()!
    for h: UInt32 in [1, 2, 4, 8, 16, 32] {
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(current, offset: 0, index: 0); encoder.setBuffer(other, offset: 0, index: 1)
        var stride = h; encoder.setBytes(&stride, length: MemoryLayout<UInt32>.stride, index: 2)
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: phases.count*n, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding(); swap(&current, &other)
    }
    command.commit(); command.waitUntilCompleted()
    guard command.status == .completed else { throw NSError(domain: "metal", code: Int(command.status.rawValue)) }
    let ptr = current.contents().bindMemory(to: Int32.self, capacity: input.count)
    return Array(UnsafeBufferPointer(start: ptr, count: input.count))
}

let args = CommandLine.arguments
guard args.count == 3 else { fatalError("usage: R24ExpurgateMetal input.json output.json") }
let inputURL = URL(fileURLWithPath: args[1]), outputURL = URL(fileURLWithPath: args[2])
let root = try! JSONSerialization.jsonObject(with: Data(contentsOf: inputURL)) as! [String: Any]
let basis = (root["original_basis_hex"] as! [String]).map { UInt64($0.dropFirst(2), radix: 16)! }
let checks = root["check_rows"] as! [Int]
let suppliedRetained = root["retained_quotient_elements"] as! [Int]
precondition(basis.count == 15 && checks.count == 2)
let expectedAffine = [UInt64.max] + (0..<6).map(truthForCoordinate)
let basisRank = rankWords(basis), affineExact = Array(basis.prefix(7)) == expectedAffine
let checkRank = rankWords(checks.map(UInt64.init))
let retained = (0..<256).filter { q in checks.allSatisfy { parity(q & $0) == 0 } }
let inputCertificateValid = basisRank == 15 && affineExact && checkRank == 2 && retained == suppliedRetained && retained.count == 64

var phases = [UInt64](), descriptors = [[String: Int]]()
var minWeight = n + 1, minWitness: [String: Int]? = nil
for quotient in retained {
    var qword: UInt64 = 0
    for j in 0..<8 where ((quotient >> j) & 1) == 1 { qword ^= basis[7+j] }
    for affine in 0..<128 {
        var word = qword
        for j in 0..<7 where ((affine >> j) & 1) == 1 { word ^= basis[j] }
        phases.append(word); descriptors.append(["quotient": quotient, "affine": affine])
        if word != 0 && word.nonzeroBitCount < minWeight { minWeight = word.nonzeroBitCount; minWitness = descriptors.last }
    }
}
precondition(phases.count == 8192)
// V is a direction space around this fixed standard bent function; the zero
// direction itself is not a Boolean bent phase.  Reconstruct it from z·w.
var baseBent: UInt64 = 0
for point in 0..<n {
    let z = point & 7, w = point >> 3
    if parity(z & w) == 1 { baseBent |= UInt64(1) << UInt64(point) }
}
let bentPhases = phases.map { $0 ^ baseBent }
guard let device = MTLCreateSystemDefaultDevice() else { fatalError("Metal unavailable") }
let library = try! device.makeLibrary(source: kernelSource, options: nil)
let pipeline = try! device.makeComputePipelineState(function: library.makeFunction(name: "fwht_stage")!)
let started = Date()
let cpu = directWalsh(bentPhases)
let gpu = try! metalFWHT(device: device, pipeline: pipeline, phases: bentPhases)
var mismatch = 0, cpuBad = 0, gpuBad = 0
for i in 0..<cpu.count {
    if cpu[i] != gpu[i] { mismatch += 1 }
    if Int(abs(cpu[i])) != 8 { cpuBad += 1 }
    if Int(abs(gpu[i])) != 8 { gpuBad += 1 }
}
let output: [String: Any] = [
    "schema": "a28-r024-expurgate-metal-v1", "input_sha256": sha256(inputURL),
    "basis_rank": basisRank, "affine_basis_exact": affineExact, "check_rows": checks, "check_rank": checkRank,
    "retained_quotient_count": retained.count, "input_retained_quotients_match": retained == suppliedRetained,
    "retained_phase_count": phases.count, "base_bent_truth_hex": String(format: "0x%016llx", baseBent),
    "minimum_nonzero_direction_weight": minWeight,
    "minimum_weight_witness": minWitness!, "cpu_oracle": "direct 64x64 character sums",
    "gpu_backend": "Metal int32 global staged ping-pong FWHT, six stages", "device": device.name,
    "registry_id": String(device.registryID), "signed_walsh_coefficients_checked": cpu.count,
    "gpu_cpu_mismatch_count": mismatch, "cpu_bad_abs8": cpuBad, "gpu_bad_abs8": gpuBad,
    "all_checks_passed": inputCertificateValid && minWeight == 16 && mismatch == 0 && cpuBad == 0 && gpuBad == 0,
    "elapsed_seconds": Date().timeIntervalSince(started)
]
try! JSONSerialization.data(withJSONObject: output, options: [.prettyPrinted, .sortedKeys]).write(to: outputURL)
if !(output["all_checks_passed"] as! Bool) { exit(10) }
print("R24 passed: 8192 phases, 524288 coefficients, min direction weight 16")
