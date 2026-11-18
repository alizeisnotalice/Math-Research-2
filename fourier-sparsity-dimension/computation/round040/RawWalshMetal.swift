import Foundation
import Metal

let shader = """
#include <metal_stdlib>
using namespace metal;

kernel void census(
    device const ushort *seeds [[buffer(0)]],
    device const ushort *characters [[buffer(1)]],
    device uint *output [[buffer(2)]],
    uint2 index [[thread_position_in_grid]]) {
    const ushort phase = ushort(index.x) ^ seeds[index.y];
    const uint dc = popcount(uint(phase));
    uint maximum = 0;
    uint active = 0;
    uint abs4 = (dc == 6u || dc == 10u) ? 1u : 0u;
    uint energy = (16u - 2u * dc) * (16u - 2u * dc);
    for (uint a = 1; a < 16; ++a) {
        const int walsh = 16 - 2 * int(popcount(uint(phase ^ characters[a])));
        const uint absolute = uint(abs(walsh));
        maximum = max(maximum, absolute);
        active += (absolute != 0);
        abs4 += (absolute == 4);
        energy += absolute * absolute;
    }
    output[index.y * 65536u + index.x] =
        dc | (maximum << 5) | (active << 10) | (abs4 << 14) | ((energy == 256u) << 19);
}
"""

func truth(_ predicate: (Int) -> Bool) -> UInt16 {
    var out: UInt16 = 0
    for x in 0..<16 where predicate(x) {
        out |= 1 << x
    }
    return out
}

let characters = (0..<16).map { a in truth { x in (a & x).nonzeroBitCount == 1 || (a & x).nonzeroBitCount == 3 } }
let quadratics = (0..<4).flatMap { i in ((i + 1)..<4).map { j in characters[1 << i] & characters[1 << j] } }
let candidates = (0..<64).map { choice in
    quadratics.enumerated().reduce(UInt16(0)) { word, term in
        (choice >> term.offset) & 1 == 1 ? word ^ term.element : word
    }
}
let seeds = [UInt16(0)] + candidates.filter { phase in
    characters.allSatisfy { c in
        let count = (phase ^ c).nonzeroBitCount
        return count == 6 || count == 10
    }
}
precondition(seeds.count == 29)
let device = MTLCreateSystemDefaultDevice()!
let options = MTLCompileOptions()
let library = try device.makeLibrary(source: shader, options: options)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "census")!)
let queue = device.makeCommandQueue()!
let input = device.makeBuffer(bytes: seeds, length: seeds.count * 2, options: .storageModeShared)!
let chars = device.makeBuffer(bytes: characters, length: characters.count * 2, options: .storageModeShared)!
let outputCount = seeds.count * 65536
let output = device.makeBuffer(length: outputCount * 4, options: .storageModeShared)!
let command = queue.makeCommandBuffer()!
let encoder = command.makeComputeCommandEncoder()!
encoder.setComputePipelineState(pipeline)
encoder.setBuffer(input, offset: 0, index: 0)
encoder.setBuffer(chars, offset: 0, index: 1)
encoder.setBuffer(output, offset: 0, index: 2)
encoder.dispatchThreads(
    MTLSize(width: 65536, height: seeds.count, depth: 1),
    threadsPerThreadgroup: MTLSize(width: min(256, pipeline.maxTotalThreadsPerThreadgroup), height: 1, depth: 1)
)
encoder.endEncoding()
command.commit()
command.waitUntilCompleted()
guard command.status == .completed else {
    fatalError("Metal kernel failed: \(String(describing: command.error))")
}

func directSummary(_ phase: UInt16) -> UInt32 {
    let dc = phase.nonzeroBitCount
    var maximum = 0
    var active = 0
    var absoluteFour = (dc == 6 || dc == 10) ? 1 : 0
    var energy = (16 - 2 * dc) * (16 - 2 * dc)
    for a in 1..<16 {
        let absolute = abs(16 - 2 * (phase ^ characters[a]).nonzeroBitCount)
        maximum = max(maximum, absolute)
        active += absolute == 0 ? 0 : 1
        absoluteFour += absolute == 4 ? 1 : 0
        energy += absolute * absolute
    }
    return UInt32(dc | (maximum << 5) | (active << 10) | (absoluteFour << 14) | ((energy == 256 ? 1 : 0) << 19))
}

let results = output.contents().bindMemory(to: UInt32.self, capacity: outputCount)
let cpu = (0..<65536).map { directSummary(UInt16($0)) }
var bentCounts: [Int] = []
var restrictedCounts: [String: Int] = [:]
var mismatches = 0
for i in 0..<seeds.count {
    var bentCount = 0
    for x in 0..<65536 {
        let summary = results[i * 65536 + x]
        if summary != cpu[x ^ Int(seeds[i])] { mismatches += 1 }
        if (summary >> 14) & 31 == 16 { bentCount += 1 }
        let offzeroMax = (summary >> 5) & 31
        if i == 0 && offzeroMax <= 4 {
            let weight = summary & 31
            let smallerWeight = min(weight, 16 - weight)
            restrictedCounts[String(smallerWeight), default: 0] += 1
        }
    }
    bentCounts.append(bentCount)
}
precondition(mismatches == 0)
precondition(bentCounts.allSatisfy { $0 == 896 })
precondition(restrictedCounts == ["0": 2, "1": 32, "2": 240, "6": 896])
let record: [String: Any] = [
    "device": device.name,
    "completed_Metal_commands": 1,
    "gpu_seconds": command.gpuEndTime - command.gpuStartTime,
    "summary_rows": outputCount,
    "walsh_coefficient_count": outputCount * 16,
    "phase_base_seeds": seeds.map { Int($0) },
    "bent_counts": bentCounts,
    "cpu_gpu_mismatches": mismatches,
    "parameter_profiles_with_max_nonzero_at_most_4_by_min_weight": restrictedCounts,
]
let json = try JSONSerialization.data(withJSONObject: record, options: [.prettyPrinted, .sortedKeys])
print(String(decoding: json, as: UTF8.self))
