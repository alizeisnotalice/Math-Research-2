import Foundation
import Metal
import CryptoKit

// Independent audit.  B is rebuilt from its defining two 3-bit conditions;
// no support/mask data from the other implementation is consumed here.
let n = 64
let b = (0..<n).filter { (($0 & 7) == 7) != (($0 >> 3) == 7) }
let candidates = 1 << b.count
let batchLimit = 4096
precondition(b.count == 14 && candidates == 16_384)

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

func sha256(_ path: String) -> String {
    let data = try! Data(contentsOf: URL(fileURLWithPath: path))
    return SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func parity(_ v: Int) -> Int { v.nonzeroBitCount & 1 }
func sign(_ mask: Int, _ j: Int) -> Int { ((mask >> j) & 1) == 0 ? 1 : -1 }

// Defining formula, evaluated directly rather than by a CPU FWHT.
func directWalsh(_ mask: Int) -> [Int32] {
    (0..<n).map { x in
        Int32((0..<b.count).reduce(0) { sum, j in
            sum + (parity(x & b[j]) == 0 ? sign(mask, j) : -sign(mask, j))
        })
    }
}

func metalFWHT(device: MTLDevice, pipeline: MTLComputePipelineState, masks: [Int]) throws -> [[Int32]] {
    var initial = [Int32](repeating: 0, count: masks.count * n)
    for r in 0..<masks.count {
        for j in 0..<b.count { initial[r * n + b[j]] = Int32(sign(masks[r], j)) }
    }
    let length = initial.count * MemoryLayout<Int32>.stride
    let options: MTLResourceOptions = .storageModeShared
    let first = initial.withUnsafeBytes { device.makeBuffer(bytes: $0.baseAddress!, length: length, options: options)! }
    let second = device.makeBuffer(length: length, options: options)!
    var input = first, output = second
    let command = device.makeCommandQueue()!.makeCommandBuffer()!
    for h: UInt32 in [1, 2, 4, 8, 16, 32] {
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(input, offset: 0, index: 0)
        encoder.setBuffer(output, offset: 0, index: 1)
        var stride = h
        encoder.setBytes(&stride, length: MemoryLayout<UInt32>.stride, index: 2)
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: masks.count * n, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding()
        swap(&input, &output)
    }
    command.commit()
    command.waitUntilCompleted()
    guard command.status == .completed else { throw NSError(domain: "metal", code: Int(command.status.rawValue)) }
    let values = input.contents().bindMemory(to: Int32.self, capacity: masks.count * n)
    return (0..<masks.count).map { Array(UnsafeBufferPointer(start: values + $0 * n, count: n)) }
}

func affineRestrictionMasks() -> Set<Int> {
    var out = Set<Int>()
    for a in 0..<n {
        for constant in 0...1 {
            var mask = 0
            for j in 0..<b.count {
                if (constant ^ parity(a & b[j])) == 1 { mask |= 1 << j }
            }
            out.insert(mask)
        }
    }
    return out
}

let arguments = CommandLine.arguments
guard arguments.count == 2 else { fatalError("usage: R22Type2CoreMetal output.json") }
let outputPath = arguments[1]
guard let device = MTLCreateSystemDefaultDevice() else { fatalError("Metal unavailable") }
let library = try! device.makeLibrary(source: kernelSource, options: nil)
let pipeline = try! device.makeComputePipelineState(function: library.makeFunction(name: "fwht_stage")!)
let start = Date()
var coefficientChecks = 0
var mismatches = 0
var invalidCPU = 0
var invalidGPU = 0
var legal: [Int] = []
var startMask = 0
while startMask < candidates {
    let end = min(startMask + batchLimit, candidates)
    let masks = Array(startMask..<end)
    let cpu = masks.map(directWalsh)
    let gpu = try! metalFWHT(device: device, pipeline: pipeline, masks: masks)
    for r in 0..<masks.count {
        var isLegal = true
        for x in 0..<n {
            coefficientChecks += 1
            if cpu[r][x] != gpu[r][x] { mismatches += 1 }
            if cpu[r][x] != 0 && abs(cpu[r][x]) != 8 { invalidCPU += 1; isLegal = false }
            if gpu[r][x] != 0 && abs(gpu[r][x]) != 8 { invalidGPU += 1 }
        }
        if isLegal { legal.append(masks[r]) }
    }
    startMask = end
}
let affine = affineRestrictionMasks()
let legalSet = Set(legal)
// The classification is an affine coset of the restricted affine-sign code,
// not necessarily that code through zero.  Anchor it at the first legal mask.
let cosetBase = legal.first ?? -1
let affineCoset = Set(affine.map { $0 ^ cosetBase })
let expected = legal.count == 128 && affine.count == 128 && legalSet == affineCoset
let result: [String: Any] = [
    "schema": "a28-r022-type2-core-metal-v1", "B": b,
    "sign_convention": "bit j equals 1 iff coefficient at B[j] is -1",
    "candidate_sign_masks_checked": candidates, "signed_walsh_coefficients_checked": coefficientChecks,
    "cpu_direct_oracle": "sum over 14 support points sign_j*(-1)^(x dot B_j), explicit dot products",
    "gpu_backend": "Metal int32 global staged ping-pong FWHT, six stages; CPU direct oracle",
    "batch_limit": batchLimit, "device": device.name, "registry_id": String(device.registryID),
    "gpu_cpu_mismatch_count": mismatches, "invalid_cpu_coefficients": invalidCPU,
    "invalid_gpu_coefficients": invalidGPU, "legal_sign_masks_14bit": legal,
    "affine_restriction_masks_14bit": affine.sorted(), "affine_coset_base_14bit": cosetBase,
    "affine_restriction_coset_masks_14bit": affineCoset.sorted(), "legal_count": legal.count,
    "affine_restriction_count": affine.count, "legal_equals_affine_restriction_coset": legalSet == affineCoset,
    "all_checks_passed": mismatches == 0 && invalidCPU == invalidGPU && expected,
    "elapsed_seconds": Date().timeIntervalSince(start)
]
try! JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys]).write(to: URL(fileURLWithPath: outputPath))
if !(result["all_checks_passed"] as! Bool) { exit(10) }
print("R22 passed: \(candidates) masks; \(coefficientChecks) signed coefficients; \(legal.count) legal masks")
