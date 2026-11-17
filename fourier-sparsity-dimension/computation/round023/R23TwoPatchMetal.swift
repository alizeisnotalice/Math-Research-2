import Foundation
import Metal
import CryptoKit

// Independent finite verification of the two selected partial-Walsh patches.
// Coordinates are x|y|z|w, with x,y each mu bits and z,w each r bits.

let kernelSource = """
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &n [[buffer(2)]],
                       constant uint &stride [[buffer(3)]],
                       uint gid [[thread_position_in_grid]]) {
    uint local = gid % n;
    uint base = gid - local;
    uint mate = base + (local ^ stride);
    int left = input[gid];
    int right = input[mate];
    output[gid] = (local & stride) ? (right - left) : (left + right);
}
"""

func parity(_ x: Int) -> Int { x.nonzeroBitCount & 1 }
func sha256(_ path: String) -> String {
    let d = try! Data(contentsOf: URL(fileURLWithPath: path))
    return SHA256.hash(data: d).map { String(format: "%02x", $0) }.joined()
}

func rank(_ vectors: [Int]) -> Int {
    var pivots: [Int: Int] = [:]
    for original in vectors {
        var x = original
        while x != 0 {
            let p = x.bitWidth - 1 - x.leadingZeroBitCount
            if let v = pivots[p] { x ^= v } else { pivots[p] = x; break }
        }
    }
    return pivots.count
}

func rankWords(_ vectors: [[UInt64]]) -> Int {
    var pivots: [Int: [UInt64]] = [:]
    for original in vectors {
        var v = original
        var found: Int? = nil
        for w in stride(from: v.count - 1, through: 0, by: -1) where v[w] != 0 {
            found = w * 64 + (63 - v[w].leadingZeroBitCount); break
        }
        while let p = found {
            if let q = pivots[p] {
                for i in 0..<v.count { v[i] ^= q[i] }
                found = nil
                for w in stride(from: v.count - 1, through: 0, by: -1) where v[w] != 0 {
                    found = w * 64 + (63 - v[w].leadingZeroBitCount); break
                }
            } else { pivots[p] = v; break }
        }
    }
    return pivots.count
}

// Definition-level O(N^2) character sums; deliberately not CPU FWHT.
func directWalsh(_ values: [Int32], phases: Int, n: Int) -> [Int32] {
    var out = [Int32](repeating: 0, count: phases * n)
    for p in 0..<phases {
        let offset = p * n
        for xi in 0..<n {
            var total = 0
            for x in 0..<n {
                total += parity(x & xi) == 0 ? Int(values[offset + x]) : -Int(values[offset + x])
            }
            out[offset + xi] = Int32(total)
        }
    }
    return out
}

func metalFWHT(device: MTLDevice, pipeline: MTLComputePipelineState, values: [Int32], phases: Int, n: Int) throws -> [Int32] {
    let length = values.count * MemoryLayout<Int32>.stride
    let options: MTLResourceOptions = .storageModeShared
    let first = values.withUnsafeBytes { device.makeBuffer(bytes: $0.baseAddress!, length: length, options: options)! }
    let second = device.makeBuffer(length: length, options: options)!
    var input = first, output = second
    let command = device.makeCommandQueue()!.makeCommandBuffer()!
    var h = 1
    while h < n {
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(input, offset: 0, index: 0)
        encoder.setBuffer(output, offset: 0, index: 1)
        var n32 = UInt32(n), h32 = UInt32(h)
        encoder.setBytes(&n32, length: MemoryLayout<UInt32>.stride, index: 2)
        encoder.setBytes(&h32, length: MemoryLayout<UInt32>.stride, index: 3)
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: phases * n, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding()
        swap(&input, &output); h <<= 1
    }
    command.commit(); command.waitUntilCompleted()
    guard command.status == .completed else { throw NSError(domain: "metal", code: Int(command.status.rawValue)) }
    let ptr = input.contents().bindMemory(to: Int32.self, capacity: values.count)
    return Array(UnsafeBufferPointer(start: ptr, count: values.count))
}

func makePatches(mu: Int, r: Int) -> [[Int]] {
    let k = mu + r, n = 1 << (2 * k), mask = (1 << mu) - 1
    var patches = [[Int](), [Int]()] 
    for xi in 0..<n {
        let x = xi & mask, y = (xi >> mu) & mask, w = xi >> (2 * mu + r)
        if w == 0 && x == 0 && (y & 1) == 1 { patches[0].append(xi) }
        if w == 0 && y == 0 && (x & 1) == 1 { patches[1].append(xi) }
    }
    precondition(Set(patches[0]).isDisjoint(with: Set(patches[1])))
    precondition(patches.allSatisfy { $0.count == (1 << k) / 2 })
    return patches
}

func normalSpace(_ support: [Int], n: Int) -> [Int] {
    let base = support.min()!
    let directions = support.map { $0 ^ base }
    return (0..<n).filter { z in directions.allSatisfy { parity(z & $0) == 0 } }
}

func audit(mu: Int, r: Int, device: MTLDevice, pipeline: MTLComputePipelineState) throws -> [String: Any] {
    let k = mu + r, n = 2 * k, N = 1 << n, q = 1 << k, mask = (1 << mu) - 1, rmask = (1 << r) - 1
    var f = [Int](repeating: 0, count: N), a = [Int](repeating: 0, count: N)
    for point in 0..<N {
        let x = point & mask, y = (point >> mu) & mask, z = (point >> (2 * mu)) & rmask, w = point >> (2 * mu + r)
        f[point] = (parity(x & y) ^ parity(z & w) ^ ((x == 0 && z == 0) ? 1 : 0))
        a[point] = (z == 0 && ((x == 0) != (y == 0))) ? 1 : 0
    }
    let patches = makePatches(mu: mu, r: r)
    var normals = [[Int]]()
    var directionRanks = [Int]()
    for patch in patches {
        let directions = patch.map { $0 ^ patch.min()! }
        let directionRank = rank(directions)
        precondition(directionRank == k - 1)
        directionRanks.append(directionRank)
        let normal = normalSpace(patch, n: N)
        precondition(rank(normal) == k + 1)
        normals.append(normal)
    }
    let combinedNormals = normals.flatMap { $0 }
    let projectionRank = rank(combinedNormals)
    precondition(projectionRank == k + mu)
    let parameterCount = 1 << (n + 2)
    let kernel = (0..<parameterCount).filter { parameter in combinedNormals.allSatisfy { parity(parameter & $0) == 0 } }
    precondition(kernel.count == (1 << (r + 2)))
    var restrictionRows = [[UInt64]]()
    for parameter in kernel {
        let lin = parameter & (N - 1), c = (parameter >> n) & 1, t = (parameter >> (n + 1)) & 1
        var row = [UInt64](repeating: 0, count: (N + 63) / 64)
        for x in 0..<N where a[x] == 1 {
            if (parity(lin & x) ^ c ^ (t * a[x])) == 1 { row[x >> 6] |= UInt64(1) << UInt64(x & 63) }
        }
        restrictionRows.append(row)
    }
    let restrictionRank = rankWords(restrictionRows)
    precondition(restrictionRank == 1)

    let phaseCount = N * 4
    var full = [Int32](repeating: 0, count: phaseCount * N)
    var partial = [Int32](repeating: 0, count: phaseCount * N)
    var linears = [Int](repeating: 0, count: phaseCount)
    var phase = 0
    for lin in 0..<N { for c in 0...1 { for t in 0...1 {
        linears[phase] = lin
        let base = phase * N
        for x in 0..<N {
            let bit = f[x] ^ parity(lin & x) ^ c ^ (t * a[x])
            full[base + x] = bit == 0 ? 1 : -1
            partial[base + x] = a[x] == 0 ? 0 : full[base + x]
        }
        phase += 1
    } } }
    let cpuFull = directWalsh(full, phases: phaseCount, n: N)
    let cpuPartial = directWalsh(partial, phases: phaseCount, n: N)
    let gpuFull = try metalFWHT(device: device, pipeline: pipeline, values: full, phases: phaseCount, n: N)
    let gpuPartial = try metalFWHT(device: device, pipeline: pipeline, values: partial, phases: phaseCount, n: N)
    var mismatches = 0, fullBadCPU = 0, fullBadGPU = 0, patchBadCPU = 0, patchBadGPU = 0
    for i in 0..<cpuFull.count {
        if cpuFull[i] != gpuFull[i] { mismatches += 1 }
        if cpuPartial[i] != gpuPartial[i] { mismatches += 1 }
        if Int(abs(cpuFull[i])) != q { fullBadCPU += 1 }
        if Int(abs(gpuFull[i])) != q { fullBadGPU += 1 }
    }
    for p in 0..<phaseCount {
        for patch in patches { for xi in patch {
            let index = p * N + (xi ^ linears[p])
            if Int(abs(cpuPartial[index])) != q { patchBadCPU += 1 }
            if Int(abs(gpuPartial[index])) != q { patchBadGPU += 1 }
        } }
    }
    return [
        "mu": mu, "free_dimension": r, "k": k, "N": N, "q": q, "phase_count": phaseCount,
        "support_A_size": a.filter { $0 == 1 }.count, "base_patches": patches,
        "patch_sizes": patches.map { $0.count }, "patch_direction_ranks": directionRanks,
        "joint_patch_codimensions": [k + 1, k + 1], "parameter_projection_rank": projectionRank,
        "parameter_kernel_size": kernel.count, "K_dimension": r + 2, "K_restriction_rank": restrictionRank,
        "full_cpu_bad_absq": fullBadCPU, "full_gpu_bad_absq": fullBadGPU,
        "partial_patch_cpu_bad_absq": patchBadCPU, "partial_patch_gpu_bad_absq": patchBadGPU,
        "gpu_cpu_mismatch_count": mismatches,
        "all_checks_passed": mismatches == 0 && fullBadCPU == 0 && fullBadGPU == 0 && patchBadCPU == 0 && patchBadGPU == 0
    ]
}

let args = CommandLine.arguments
guard args.count == 2 else { fatalError("usage: R23TwoPatchMetal output.json") }
guard let device = MTLCreateSystemDefaultDevice() else { fatalError("Metal unavailable") }
let library = try! device.makeLibrary(source: kernelSource, options: nil)
let pipeline = try! device.makeComputePipelineState(function: library.makeFunction(name: "fwht_stage")!)
let started = Date()
let rows = try! [(3, 0), (3, 1), (4, 0)].map { try audit(mu: $0.0, r: $0.1, device: device, pipeline: pipeline) }
let result: [String: Any] = [
    "schema": "a28-r023-two-patch-metal-v1", "device": device.name, "registry_id": String(device.registryID),
    "cpu_oracle": "explicit definition-level O(N^2) character sums for full and A-weighted partial spectra",
    "gpu_backend": "Metal int32 global staged ping-pong FWHT, elementwise compared to CPU oracle",
    "rows": rows, "all_checks_passed": rows.allSatisfy { $0["all_checks_passed"] as! Bool },
    "elapsed_seconds": Date().timeIntervalSince(started)
]
try! JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys]).write(to: URL(fileURLWithPath: args[1]))
if !(result["all_checks_passed"] as! Bool) { exit(10) }
print("R23 passed all three parameter sets")
