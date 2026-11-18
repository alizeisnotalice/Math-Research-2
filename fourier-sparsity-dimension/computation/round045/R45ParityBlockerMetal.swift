import CryptoKit
import Dispatch
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
    var commandSpans: [(Double, Double)] = []

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
        commandSpans.append((command.gpuStartTime, command.gpuEndTime))
        gpuSeconds += command.gpuEndTime - command.gpuStartTime
        let address = read.contents().bindMemory(to: Int32.self, capacity: input.count)
        return Array(UnsafeBufferPointer(start: address, count: input.count))
    }
}

struct RNG {
    var state: UInt64
    mutating func next() -> UInt64 {
        state &+= 0x9e3779b97f4a7c15
        var output = state
        output = (output ^ (output >> 30)) &* 0xbf58476d1ce4e5b9
        output = (output ^ (output >> 27)) &* 0x94d049bb133111eb
        return output ^ (output >> 31)
    }
}

func cpuFWHT(_ input: [Int32]) -> [Int32] {
    var result = input
    var width = 1
    while width < result.count {
        for start in stride(from: 0, to: result.count, by: 2 * width) {
            for offset in start..<(start + width) {
                let left = result[offset]
                let right = result[offset + width]
                result[offset] = left + right
                result[offset + width] = left - right
            }
        }
        width <<= 1
    }
    return result
}

final class CpuTask: @unchecked Sendable {
    private let lock = NSLock()
    private var output: [Int32] = []
    private var interval: (Double, Double) = (0, 0)

    func execute(_ inputs: [[Int32]]) {
        let start = ProcessInfo.processInfo.systemUptime
        let transformed = inputs.flatMap(cpuFWHT)
        let end = ProcessInfo.processInfo.systemUptime
        lock.lock()
        output = transformed
        interval = (start, end)
        lock.unlock()
    }

    func result() -> ([Int32], (Double, Double)) {
        lock.lock()
        defer { lock.unlock() }
        return (output, interval)
    }
}

struct PathInfo {
    let nodes: [Int]
    let choices: [Int]
}

struct ShiftAudit {
    let delta: Int
    let chosen: Int
    let count: Int
    let largest: Int
    let feasible: Int
}

struct Fixture {
    let addressBits: Int
    let dataBits: Int
    let q: [Int]
    let leaves: [[Int32]]
    let description: String
    let components: Int
}

func paths(_ internalNodes: Int) -> [PathInfo] {
    var result = [PathInfo(nodes: [], choices: [])]
    for v in 1..<internalNodes {
        let parent = (v - 1) / 2
        let choice = (v & 1) == 0 ? 1 : 0
        result.append(PathInfo(nodes: result[parent].nodes + [parent],
                               choices: result[parent].choices + [choice]))
    }
    return result
}

func frequencySpan(_ masks: [Int], _ choices: [Int]) -> [(Int, Int)] {
    var span = [(0,0)]
    for (mask, bit) in zip(masks, choices) {
        span.append(contentsOf: span.map { ($0.0 ^ mask, $0.1 ^ bit) })
    }
    return span
}

func siblingBuckets(_ q: [Int], _ path: PathInfo) -> [[Int]] {
    func subtreeMasks(_ node: Int) -> [Int] {
        if node >= q.count { return [] }
        return [q[node]] + subtreeMasks(2 * node + 1) + subtreeMasks(2 * node + 2)
    }
    return zip(path.nodes, path.choices).map { parent, choice in
        subtreeMasks(2 * parent + 1 + (choice ^ 1))
    }
}

func minimumBlocker(
    _ q: [Int], _ v: Int, _ path: PathInfo, _ addressCount: Int
) throws -> ShiftAudit {
    let buckets = siblingBuckets(q, path)
    let ancestorLabels = path.nodes.map { q[$0] }
    var minimum = Int.max
    var maximum = 0
    var bestDelta = 0
    var chosen = 0
    var count = 0
    for delta in 1..<addressCount where parity(delta & q[v]) == 1 {
        if ancestorLabels.contains(where: { parity(delta & $0) != 0 }) { continue }
        count += 1
        var places = 0
        for (i, masks) in buckets.enumerated() where masks.contains(where: { parity(delta & $0) != 0 }) {
            places |= 1 << i
        }
        let size = places.nonzeroBitCount
        if size < minimum {
            minimum = size
            chosen = places
            bestDelta = delta
        }
        maximum = max(maximum, size)
    }
    try check(count == addressCount >> (path.nodes.count + 1), "Wrong number of feasible dual shifts")
    return ShiftAudit(delta: bestDelta, chosen: chosen, count: minimum, largest: maximum, feasible: count)
}

func rank(_ values: Set<Int>) -> Int {
    var pivots: [Int: Int] = [:]
    for x in values {
        var y = x
        while y != 0 {
            let bit = Int.bitWidth - y.leadingZeroBitCount - 1
            if let previous = pivots[bit] { y ^= previous }
            else {
                pivots[bit] = y
                break
            }
        }
    }
    return pivots.count
}

func randomLabels(_ bits: Int, _ internalNodes: Int, _ rng: inout RNG) -> [Int] {
    let path = paths(internalNodes)
    var q: [Int] = []
    for v in 0..<internalNodes {
        let ancestors = path[v].nodes.map { q[$0] }
        let span = Set(frequencySpan(ancestors, path[v].choices).map { $0.0 })
        let available = (1..<(1 << bits)).filter { !span.contains($0) }
        q.append(available[Int(rng.next() % UInt64(available.count))])
    }
    return q
}

func noncharacterLeaves(_ leaves: Int, _ dataBits: Int, _ rng: inout RNG) -> [[Int32]] {
    (0..<leaves).map { _ in
        while true {
            let table = (0..<(1 << dataBits)).map { _ in sign(Int(rng.next() & 1)) }
            if directWalsh(table).filter({ $0 != 0 }).count > 1 { return table }
        }
    }
}

func independentComponents(_ t: Int) -> Fixture {
    var q = [Int](repeating: 0, count: 8 * t - 1)
    for index in 0..<(t - 1) { q[index] = 1 << index }
    for i in 0..<t {
        let root = t - 1 + i
        let left = 2 * root + 1
        let right = left + 1
        let privateRoot = 1 << (t - 1 + 3 * i)
        let b = privateRoot << 1
        let c = b << 1
        q[root] = privateRoot
        q[left] = b
        q[right] = c
        q[2 * left + 1] = c
        q[2 * left + 2] = b ^ c
        q[2 * right + 1] = b
        q[2 * right + 2] = b ^ c
    }
    let leaves = (0..<(8*t)).map { index -> [Int32] in
        let component = index / 8
        let local = index % 8
        let negative = [3,5,6,7].contains(local)
        let table = (0..<4).map { data in sign((data ^ (component & 3)) == 3 ? 1 : 0) }
        return negative ? table.map { -$0 } : table
    }
    return Fixture(addressBits: 4 * t - 1, dataBits: 2, q: q, leaves: leaves,
                   description: "independent selector components with reversed and dependent labels",
                   components: t)
}

struct Counters {
    var feasibleShifts = 0
    var localBounds = 0
    var confinementChecks = 0
    var packingChecks = 0
    var packingLabelsAtCapacity = 0
    var skewedNodes = 0
    var zeroBlockerActive = 0
    var zeroBlockerAll = 0
    var highestMin = 0
    var globallyDependent = 0
    var orderingCycles = 0
    var bandsLosingCoefficients = 0
}

func audit(_ f: Fixture, _ totals: inout Counters) throws -> ([Int32], [Int32], [String: Any]) {
    let path = paths(f.q.count)
    let addressCount = 1 << f.addressBits
    let dataCount = 1 << f.dataBits
    let n = addressCount * dataCount
    let scale: Int32 = 64
    var allNodes: [[String: Any]] = []
    var allShifts: [ShiftAudit] = []
    var bands: [[(Int, Int)]] = []
    var edges = Set<Int>()
    try check(f.leaves.allSatisfy { directWalsh($0).filter { $0 != 0 }.count > 1 },
              "Component fixture includes a character leaf")
    for (v, info) in path.enumerated() {
        let ancestorLabels = info.nodes.map { f.q[$0] }
        let localSpan = frequencySpan(ancestorLabels, info.choices)
        try check(Set(localSpan.map { $0.0 }).count == 1 << info.nodes.count,
                  "Ancestral queries are dependent")
        try check(!localSpan.contains { $0.0 == f.q[v] }, "Path query is dependent")
        bands.append(localSpan.map { ($0.0 ^ f.q[v], $0.1) })
        allShifts.append(try minimumBlocker(f.q, v, info, addressCount))
        for label in ancestorLabels { edges.insert(label * addressCount + f.q[v]) }
    }
    let cycle = edges.contains { edges.contains(($0 % addressCount) * addressCount + $0 / addressCount) }
    let distinct = Set(f.q)
    let rankOfLabels = rank(distinct)
    if cycle { totals.orderingCycles += 1 }
    if rankOfLabels < distinct.count { totals.globallyDependent += 1 }
    let selector = (0..<addressCount).map { address -> Int in
        var v = 0
        while v < f.q.count { v = 2 * v + 1 + parity(f.q[v] & address) }
        return v - f.q.count
    }
    let truth = (0..<n).map { point in f.leaves[selector[point % addressCount]][point / addressCount] }
    var means = [Int32](repeating: 0, count: dataCount)
    var signals = f.q.map { _ in [Int32](repeating: 0, count: dataCount) }
    for w in 0..<dataCount {
        var rows = [Int32](repeating: 0, count: f.q.count) + f.leaves.map { scale * $0[w] }
        for v in f.q.indices.reversed() {
            let l = rows[2 * v + 1], r = rows[2 * v + 2]
            try check((l+r) % 2 == 0, "Local scale insufficient")
            rows[v] = (l + r) / 2
            signals[v][w] = (l - r) / 2
        }
        means[w] = rows[0]
        let direct = truth[(w * addressCount)..<((w+1) * addressCount)].reduce(Int32(0), +)
        try check(direct * scale == means[w] * Int32(addressCount), "Root subtree mean mismatch")
    }
    let localSpectra = signals.map(directWalsh)
    let meanSpectrum = directWalsh(means)
    var full = [Int32](repeating: 0, count: n)
    for w in 0..<dataCount {
        let numerator = meanSpectrum[w] * Int32(addressCount)
        try check(numerator % scale == 0, "Mean Walsh not integral")
        full[w * addressCount] = numerator / scale
    }
    for v in f.q.indices {
        let coefficientScale = Int32(addressCount >> path[v].nodes.count)
        for (alpha, twist) in bands[v] {
            for (beta, value) in localSpectra[v].enumerated() {
                let local = coefficientScale * sign(twist) * value
                try check(local % scale == 0, "Local band Walsh not integral")
                full[alpha + addressCount * beta] += local / scale
            }
        }
    }
    let s = full.filter { $0 != 0 }.count
    let room = 1 << (path.last!.nodes.count)
    var capacities = [Int](repeating: 0, count: addressCount)
    var weightedNumerator = 0
    var localZeroActive = 0
    var localZeroAll = 0
    var localSkewed = 0
    var locallyVanished = 0
    for v in f.q.indices {
        let record = allShifts[v]
        let p = path[v]
        let depth = p.nodes.count
        let ell = record.count
        let local = localSpectra[v]
        let sDelta = local.filter { $0 != 0 }.count
        let localSize = sDelta << depth
        let outerScale = Int32(addressCount >> depth)
        var sBand = 0
        let selected = p.nodes.indices.filter { (record.chosen >> $0) & 1 == 1 }
        let gate = frequencySpan(selected.map { f.q[p.nodes[$0]] }, selected.map { p.choices[$0] })
        for (alpha, twist) in bands[v] {
            capacities[alpha] += room >> ell
            for (beta, localValue) in local.enumerated() {
                let index = alpha + addressCount * beta
                if full[index] != 0 { sBand += 1 }
                if localValue != 0 && full[index] == 0 { locallyVanished += 1 }
                let afterGate = gate.reduce(into: Int32(0)) { sum, shift in
                    sum += sign(shift.1) * full[(alpha ^ shift.0) + addressCount * beta]
                }
                let expectedNumerator = outerScale * sign(twist) * localValue * Int32(1 << ell)
                try check(afterGate * scale == expectedNumerator,
                          "Blocker multiplication fails to exactly recover the local spectrum")
                totals.confinementChecks += 1
            }
        }
        try check((sBand << ell) >= localSize, "True band sparsity below local bound")
        weightedNumerator += (sDelta << depth) * ((room * room) >> (2 * ell))
        totals.localBounds += 1
        totals.feasibleShifts += record.feasible
        if record.largest > ell { localSkewed += 1 }
        if ell == 0 {
            localZeroAll += 1
            if sDelta > 0 { localZeroActive += 1 }
        }
        totals.highestMin = max(totals.highestMin, ell)
        if f.components != 0 {
            allNodes.append([
                "node": v, "query": f.q[v], "depth": depth, "lca_blockers_minimum": ell,
                "lca_blockers_maximum_over_shifts": record.largest, "minimizing_shift": record.delta,
                "minimizing_blocker_ancestors": selected.map { p.nodes[$0] },
                "band_addresses": bands[v].map { $0.0 }.sorted(),
                "delta_sparsity": sDelta, "true_band_sparsity": sBand,
            ])
        }
    }
    try check(capacities.allSatisfy { $0 <= room }, "Weighted band capacities exceed one")
    try check(weightedNumerator <= s * room * room, "Weighted global sparsity bound fails")
    totals.packingChecks += addressCount
    totals.packingLabelsAtCapacity += capacities.filter { $0 == room }.count
    totals.skewedNodes += localSkewed
    totals.zeroBlockerAll += localZeroAll
    totals.zeroBlockerActive += localZeroActive
    totals.bandsLosingCoefficients += locallyVanished
    let info: [String: Any] = [
        "description": f.description, "components": f.components,
        "address_bits": f.addressBits, "data_bits": f.dataBits,
        "query_rank": rankOfLabels, "number_of_distinct_query_labels": distinct.count,
        "ordering_cycle_between_two_labels": cycle,
        "zero_blocker_all_nodes": localZeroAll, "zero_blocker_active_nodes": localZeroActive,
        "nodes_with_shift_dependent_blocker_count": localSkewed,
        "global_sparsity": s, "weighted_bound_numerator": weightedNumerator,
        "weighted_bound_denominator": room * room,
        "nodes": allNodes,
    ]
    return (truth, full, info)
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R45ParityBlockerMetal output.json")
    let seed: UInt64 = 0xa280450ddbb12345
    var rng = RNG(state: seed)
    let metal = try MetalWalsh()
    var cases: [Fixture] = []
    for root in 1..<8 {
        for left in 1..<8 where left != root {
            for right in 1..<8 where right != root {
                cases.append(Fixture(addressBits: 3, dataBits: 2, q: [root, left, right],
                                     leaves: noncharacterLeaves(4, 2, &rng),
                                     description: "exhaustive depth-two query labels in F2^3", components: 0))
            }
        }
    }
    for _ in 0..<64 {
        let q = randomLabels(5, 7, &rng)
        let leaves = noncharacterLeaves(8, 3, &rng)
        cases.append(Fixture(addressBits: 5, dataBits: 3, q: q, leaves: leaves,
                             description: "seeded arbitrary path-independent depth-three labels", components: 0))
    }
    cases += [1, 2, 4].map(independentComponents)
    var totals = Counters()
    var details: [[String: Any]] = []
    var batches: [Int: [([Int32], [Int32])]] = [:]
    for f in cases {
        let (truth, spectrum, info) = try audit(f, &totals)
        if f.components != 0 {
            try check(info["ordering_cycle_between_two_labels"] as? Bool == true,
                      "Fixture accidentally globally orderable")
            try check((info["query_rank"] as! Int) < (info["number_of_distinct_query_labels"] as! Int),
                      "Component labels are distinct but linearly independent")
            try check((info["zero_blocker_active_nodes"] as! Int) == f.components,
                      "Expected one active zero-blocker root per component")
            details.append(info)
        }
        batches[truth.count, default: []].append((truth, spectrum))
    }

    var coefficients = 0
    var cpuGpuOverlap: Double = 0
    for domain in batches.keys.sorted() {
        let batch = batches[domain]!
        let cpuWorker = CpuTask()
        let group = DispatchGroup()
        let input = batch.map { $0.0 }
        DispatchQueue.global(qos: .userInitiated).async(group: group) { cpuWorker.execute(input) }
        let gpu = try metal.transform(input.flatMap { $0 }, rowSize: domain)
        group.wait()
        let (cpu, span) = cpuWorker.result()
        let gpuSpan = metal.commandSpans.last!
        try check(cpu == gpu && cpu == batch.flatMap { $0.1 }, "CPU/Metal/Haar coefficients differ")
        cpuGpuOverlap += max(0, min(span.1, gpuSpan.1) - max(span.0, gpuSpan.0))
        coefficients += cpu.count
    }
    let output: [String: Any] = [
        "schema": "a28-r045-parity-query-lca-blocker-metal-v1",
        "seed_hex": String(seed, radix: 16),
        "number_of_trees": cases.count,
        "exhaustive_three_address_variable_label_configurations": 252,
        "random_five_address_variable_cases": 64,
        "all_feasible_dual_shift_evaluations": totals.feasibleShifts,
        "nodes_where_blocker_count_depends_on_shift": totals.skewedNodes,
        "maximum_minimal_blocker_count": totals.highestMin,
        "local_band_lower_bound_checks": totals.localBounds,
        "exact_blocker_multiplication_coefficient_checks": totals.confinementChecks,
        "address_frequency_capacity_checks": totals.packingChecks,
        "address_frequencies_capacity_equals_one": totals.packingLabelsAtCapacity,
        "trees_with_dependence_of_distinct_query_labels": totals.globallyDependent,
        "trees_with_two_label_order_cycles": totals.orderingCycles,
        "zero_blocker_nodes_all": totals.zeroBlockerAll,
        "zero_blocker_nodes_active": totals.zeroBlockerActive,
        "local_nonzero_coefficient_globally_vanishing_occurrences": totals.bandsLosingCoefficients,
        "full_coefficients_cpu_metal_and_haar": coefficients,
        "full_coefficient_mismatches": 0,
        "component_examples": details,
        "device": metal.device.name,
        "command_buffers": metal.dispatches,
        "gpu_command_seconds": metal.gpuSeconds,
        "measured_cpu_gpu_overlap_seconds": cpuGpuOverlap,
        "cpu_backend": "Int32 direct local Walsh and full FWHT",
        "gpu_backend": "Int32 Metal FWHT",
        "all_checks_passed": true,
    ]
    let bytes = try JSONSerialization.data(withJSONObject: output, options: [.prettyPrinted, .sortedKeys])
    try (bytes + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R45 verified \(cases.count) parity trees and \(coefficients) full integer coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R45 verifier error: \(error)\n").utf8))
    exit(1)
}
