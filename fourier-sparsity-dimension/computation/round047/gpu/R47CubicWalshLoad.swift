import Foundation
import Metal

let metalSource = """
#include <metal_stdlib>
using namespace metal;

uint mix(uint x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    return x ^ (x >> 16);
}

int eval(device const uint *pointMask, uint id, uint c1, uint c2, uint c3, uint c4) {
    uint n = popcount(pointMask[0] & id)
           + popcount(pointMask[1] & c1)
           + popcount(pointMask[2] & c2)
           + popcount(pointMask[3] & c3)
           + popcount(pointMask[4] & c4);
    return (n & 1u) ? -1 : 1;
}

struct Params { uint first; uint seed; uint count; };

kernel void cubic_scan(device const uint *monomialMasks [[buffer(0)]],
                       device uint *candidateAbsCounts [[buffer(1)]],
                       device uint *sparsities [[buffer(2)]],
                       device int *samples [[buffer(3)]],
                       constant Params &param [[buffer(4)]],
                       uint3 g [[threadgroup_position_in_grid]],
                       uint lane [[thread_index_in_threadgroup]]) {
    threadgroup int walsh[512];
    threadgroup atomic_uint histogram[257];
    uint id = param.first + g.x;
    uint c1 = mix(id ^ param.seed ^ 0x9e3779b9u);
    uint c2 = mix(id ^ param.seed ^ 0x3c6ef372u);
    uint c3 = mix(id ^ param.seed ^ 0xdaa66d2bu);
    uint c4 = mix(id ^ param.seed ^ 0x78dde6e4u);
    walsh[lane] = eval(monomialMasks + lane * 5u, id, c1, c2, c3, c4);
    walsh[lane+256u] = eval(monomialMasks + (lane+256u) * 5u, id, c1, c2, c3, c4);
    atomic_store_explicit(histogram+lane, 0u, memory_order_relaxed);
    if (lane == 0) atomic_store_explicit(histogram+256u, 0u, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint h=1; h<512; h*=2) {
        uint left = ((lane & ~(h-1u)) << 1) | (lane & (h-1u));
        uint right = left+h;
        int a = walsh[left], b = walsh[right];
        walsh[left] = a+b;
        walsh[right] = a-b;
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    int x = walsh[lane], y = walsh[lane+256u];
    atomic_fetch_add_explicit(histogram + uint(abs(x))/2u, 1u, memory_order_relaxed);
    atomic_fetch_add_explicit(histogram + uint(abs(y))/2u, 1u, memory_order_relaxed);
    if (g.x == 0 || g.x + 1u == param.count) {
        uint sample = g.x == 0 ? 0u : 512u;
        samples[sample+lane] = x;
        samples[sample+lane+256u] = y;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint row = g.x * 257u;
    candidateAbsCounts[row+lane] = atomic_load_explicit(histogram+lane, memory_order_relaxed);
    if (lane == 0) {
        candidateAbsCounts[row+256u] = atomic_load_explicit(histogram+256u, memory_order_relaxed);
        sparsities[g.x] = 512u - atomic_load_explicit(histogram, memory_order_relaxed);
    }
}

kernel void reduce_abs(device const uint *matrix [[buffer(0)]],
                       device uint *binTotals [[buffer(1)]],
                       constant uint &count [[buffer(2)]],
                       uint3 bin [[threadgroup_position_in_grid]],
                       uint lane [[thread_index_in_threadgroup]]) {
    threadgroup uint tmp[256];
    uint total = 0;
    for (uint i=lane; i<count; i+=256) total += matrix[257u*i+bin.x];
    tmp[lane] = total;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint h=128; h>0; h/=2) {
        if (lane < h) tmp[lane] += tmp[lane+h];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (lane == 0) binTotals[bin.x] = tmp[0];
}
"""

struct ScanError: Error { let reason: String }
func require(_ condition: Bool, _ message: String) throws {
    guard condition else { throw ScanError(reason: message) }
}

struct Params {
    var first: UInt32
    var seed: UInt32
    var count: UInt32
}

func mix(_ value: UInt32) -> UInt32 {
    var x = value
    x ^= x >> 16
    x = x &* 0x7feb352d
    x ^= x >> 15
    x = x &* 0x846ca68b
    return x ^ (x >> 16)
}

func monomialMasks() -> [UInt32] {
    let masks3 = (0..<512).filter { $0.nonzeroBitCount == 3 }
    let lower = (0..<512).filter { $0.nonzeroBitCount < 3 }
    let all = masks3 + lower
    precondition(all.count == 130)
    var masks = [UInt32](repeating: 0, count: 512*5)
    for point in 0..<512 {
        for (bit, monomial) in all.enumerated() where point & monomial == monomial {
            masks[5*point+bit/32] |= (1 << (bit & 31))
        }
    }
    return masks
}

func reference(_ id: UInt32, _ seed: UInt32, _ masks: [UInt32]) -> [Int32] {
    let words = [id, mix(id ^ seed ^ 0x9e3779b9), mix(id ^ seed ^ 0x3c6ef372),
                 mix(id ^ seed ^ 0xdaa66d2b), mix(id ^ seed ^ 0x78dde6e4)]
    let truth = (0..<512).map { point -> Int32 in
        let parity = (0..<5).reduce(0) { $0 + (words[$1] & masks[point*5+$1]).nonzeroBitCount }
        return parity & 1 == 0 ? 1 : -1
    }
    // Direct Walsh sum is independent of the GPU butterfly algorithm.
    return (0..<512).map { alpha in
        truth.indices.reduce(into: Int32(0)) { sum, x in
            sum += (alpha & x).nonzeroBitCount & 1 == 0 ? truth[x] : -truth[x]
        }
    }
}

final class Slot {
    let matrix: MTLBuffer
    let sparsity: MTLBuffer
    let samples: MTLBuffer
    let sums: MTLBuffer
    var command: MTLCommandBuffer?
    var first: UInt32 = 0
    var number: Int = 0

    init(_ device: MTLDevice, _ batch: Int) {
        matrix = device.makeBuffer(length: batch * 257 * 4, options: .storageModePrivate)!
        sparsity = device.makeBuffer(length: batch * 4, options: .storageModeShared)!
        samples = device.makeBuffer(length: 1024 * 4, options: .storageModeShared)!
        sums = device.makeBuffer(length: 257 * 4, options: .storageModeShared)!
    }
}

func run() throws {
    let args = CommandLine.arguments
    try require(args.count == 7,
        "usage: R47CubicWalshLoad output.json first_id count seed batch_size sample_stride")
    guard let first = UInt32(args[2]), let count = Int(args[3]), count > 0,
          let seed = UInt32(args[4], radix: 16), let batch = Int(args[5]),
          batch >= 256 && batch <= 32768, let strideSample = Int(args[6]), strideSample >= 1 else {
        throw ScanError(reason: "Invalid sweep parameters")
    }
    try require(UInt64(first) + UInt64(count) < UInt64(UInt32.max), "ANF counter range overflows")
    guard let device = MTLCreateSystemDefaultDevice() else { throw ScanError(reason: "No Metal device") }
    let library = try device.makeLibrary(source: metalSource, options: nil)
    let scan = try device.makeComputePipelineState(function: library.makeFunction(name: "cubic_scan")!)
    let reduce = try device.makeComputePipelineState(function: library.makeFunction(name: "reduce_abs")!)
    let queues = (0..<2).map { _ in device.makeCommandQueue()! }
    let ring = (0..<4).map { _ in Slot(device, batch) }
    let masks = monomialMasks()
    let maskBuffer = masks.withUnsafeBytes {
        device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
    }
    var hist = [UInt64](repeating: 0, count: 257)
    var sparsityHist = [UInt64](repeating: 0, count: 513)
    var kernelTimes: [(Double,Double)] = []
    var sampled: [UInt32] = []
    var sampledCoef = 0
    var samplesSkipped = 0
    var reaped = 0
    let uptimeStart = ProcessInfo.processInfo.systemUptime
    let epochOffset = Date().timeIntervalSince1970 - uptimeStart
    print("Started candidate scan on \(device.name), pid=\(ProcessInfo.processInfo.processIdentifier)")

    func collect(_ slot: Slot) throws {
        guard let command = slot.command else { return }
        command.waitUntilCompleted()
        if let error = command.error { throw error }
        try require(command.status == .completed, "GPU command did not complete")
        kernelTimes.append((command.gpuStartTime, command.gpuEndTime))
        let absolute = slot.sums.contents().bindMemory(to: UInt32.self, capacity: 257)
        let sparsities = slot.sparsity.contents().bindMemory(to: UInt32.self, capacity: slot.number)
        let samples = slot.samples.contents().bindMemory(to: Int32.self, capacity: 1024)
        var countAbs: UInt64 = 0
        var energy: UInt64 = 0
        for i in hist.indices {
            let count = UInt64(absolute[i])
            hist[i] += count
            countAbs += count
            energy += count * UInt64(4*i*i)
        }
        try require(countAbs == UInt64(slot.number)*512
                    && energy == UInt64(slot.number)*512*512, "Histogram count or Parseval fails")
        for i in 0..<slot.number {
            let value = Int(sparsities[i])
            try require((1...512).contains(value), "Invalid sparsity")
            sparsityHist[value] += 1
        }
        if reaped % strideSample == 0 {
            for sampleNumber in 0..<2 {
                let offset = sampleNumber == 0 ? 0 : slot.number-1
                if offset == 0 && sampleNumber == 1 { continue }
                let id = slot.first + UInt32(offset)
                let cpu = reference(id, seed, masks)
                try require(cpu.indices.allSatisfy { cpu[$0] == samples[512*sampleNumber+$0] },
                            "Exact sampled CPU/GPU coefficient mismatch")
                try require(cpu.filter { $0 != 0 }.count == Int(sparsities[offset]),
                            "Sampled sparsity mismatch")
                sampled.append(id)
                sampledCoef += 512
            }
        } else { samplesSkipped += 1 }
        reaped += 1
        slot.command = nil
    }

    var next = 0
    var submitted = 0
    while next < count {
        let slot = ring[submitted % ring.count]
        try collect(slot)
        let items = min(batch, count-next)
        let command = queues[submitted % queues.count].makeCommandBuffer()!
        let compute = command.makeComputeCommandEncoder()!
        compute.setComputePipelineState(scan)
        compute.setBuffer(maskBuffer, offset: 0, index: 0)
        compute.setBuffer(slot.matrix, offset: 0, index: 1)
        compute.setBuffer(slot.sparsity, offset: 0, index: 2)
        compute.setBuffer(slot.samples, offset: 0, index: 3)
        var params = Params(first: first + UInt32(next), seed: seed, count: UInt32(items))
        compute.setBytes(&params, length: MemoryLayout<Params>.stride, index: 4)
        let width = MTLSize(width: 256, height: 1, depth: 1)
        compute.dispatchThreadgroups(MTLSize(width: items, height: 1, depth: 1), threadsPerThreadgroup: width)
        compute.endEncoding()
        let reduction = command.makeComputeCommandEncoder()!
        reduction.setComputePipelineState(reduce)
        reduction.setBuffer(slot.matrix, offset: 0, index: 0)
        reduction.setBuffer(slot.sums, offset: 0, index: 1)
        var n = UInt32(items)
        reduction.setBytes(&n, length: 4, index: 2)
        reduction.dispatchThreadgroups(MTLSize(width: 257, height: 1, depth: 1), threadsPerThreadgroup: width)
        reduction.endEncoding()
        slot.first = params.first
        slot.number = items
        slot.command = command
        command.commit()
        next += items
        submitted += 1
        if submitted % 64 == 0 { print("Submitted \(next)/\(count) cubic ANFs") }
    }
    for slot in ring { try collect(slot) }
    let duration = ProcessInfo.processInfo.systemUptime - uptimeStart
    try require(reaped == submitted && sparsityHist.reduce(0,+) == UInt64(count)
                && hist.reduce(0,+) == UInt64(count)*512, "Lost candidates")
    let data: [String: Any] = [
        "schema": "a28-r047-cubic-anf-gpu-sweep-v1",
        "device": device.name, "pid": ProcessInfo.processInfo.processIdentifier,
        "first_candidate_id": first, "candidate_count": count,
        "seed_hex": String(seed,radix:16), "batch_size": batch, "sample_stride": strideSample,
        "distinctness": "first 32 cubic ANF coefficients equal the unique nonzero candidate counter",
        "gpu_dispatch_wall_seconds": duration,
        "gpu_command_intervals_utc_unix": kernelTimes.map { [$0.0 + epochOffset,$0.1 + epochOffset] },
        "gpu_command_buffer_seconds_total": kernelTimes.reduce(0) { $0 + ($1.1-$1.0) },
        "gpu_command_buffers": submitted, "metal_command_queues": queues.count,
        "inflight_command_buffers": ring.count,
        "candidate_histogram_buffer_bytes_per_slot": batch*257*4,
        "walsh_coefficients_evaluated": UInt64(count)*512,
        "absolute_walsh_histogram_even_bins_0_to_512": hist,
        "sparsity_histogram": sparsityHist,
        "sampled_ids": sampled, "sampled_exact_cpu_gpu_coefficients": sampledCoef,
        "batches_without_cpu_sample_by_design": samplesSkipped,
        "parseval_and_sample_checks_passed": true,
    ]
    let bytes = try JSONSerialization.data(withJSONObject: data, options: [.prettyPrinted,.sortedKeys])
    try (bytes + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("Completed \(count) cubic ANFs, \(count*512) coefficients, \(sampledCoef) independently checked.")
}

do { try run() }
catch {
    FileHandle.standardError.write(Data(("GPU cubic scan error: \(error)\n").utf8))
    exit(1)
}
