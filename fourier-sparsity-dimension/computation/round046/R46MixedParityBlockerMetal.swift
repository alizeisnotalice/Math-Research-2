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

struct Label {
    let q: Int
    let b: Int
    let constant: Int

    var isPrivate: Bool { q != 0 }

    func bit(_ x: Int, _ w: Int) -> Int {
        parity(q & x) ^ parity(b & w) ^ constant
    }

    func jointMask(_ xBits: Int) -> Int {
        q | (b << xBits)
    }
}

struct Fixture {
    let kind: String
    let xBits: Int
    let wBits: Int
    let labels: [Label]
    let leaves: [[Int32]]
}

func blockerBuckets(_ labels: [Label], _ path: PathInfo) -> [[Int]] {
    func privateQueries(_ v: Int) -> [Int] {
        if v >= labels.count { return [] }
        let current = labels[v].isPrivate ? [labels[v].q] : []
        return current + privateQueries(2 * v + 1) + privateQueries(2 * v + 2)
    }
    return zip(path.nodes, path.choices).map { v, decision in
        privateQueries(2 * v + 1 + (decision ^ 1))
    }
}

func shifts(_ labels: [Label], _ v: Int, _ path: PathInfo, _ xSize: Int) throws -> ShiftAudit {
    let privateAncestors = path.nodes.filter { labels[$0].isPrivate }.map { labels[$0].q }
    let buckets = blockerBuckets(labels, path)
    var minimum = Int.max
    var maximum = 0
    var chosen = 0
    var optimizing = 0
    var tested = 0
    for delta in 1..<xSize where parity(delta & labels[v].q) == 1 {
        if privateAncestors.contains(where: { parity(delta & $0) != 0 }) { continue }
        tested += 1
        var mask = 0
        for (i, others) in buckets.enumerated() where others.contains(where: { parity(delta & $0) != 0 }) {
            mask |= 1 << i
        }
        let count = mask.nonzeroBitCount
        maximum = max(maximum, count)
        if count < minimum {
            minimum = count
            chosen = mask
            optimizing = delta
        }
    }
    try check(tested == xSize >> (privateAncestors.count + 1), "Private path shift count mismatch")
    return ShiftAudit(delta: optimizing, chosen: chosen, count: minimum, largest: maximum, feasible: tested)
}

func leaves(_ count: Int, _ size: Int, _ rng: inout RNG) -> [[Int32]] {
    (0..<count).map { _ in
        while true {
            let table = (0..<size).map { _ in sign(Int(rng.next() & 1)) }
            if directWalsh(table).filter({ $0 != 0 }).count > 1 { return table }
        }
    }
}

func fixedFixtures(_ wBits: Int) -> [Fixture] {
    let size = 1 << wBits
    let a = (0..<size).map { sign(($0 & 3) == 3 ? 1 : 0) }
    let b = (0..<size).map { sign(($0 & 3) == 2 ? 1 : 0) }
    let negA = a.map { -$0 }
    let negB = b.map { -$0 }
    func fixed(_ q: Int, _ b: Int, _ c: Int = 0) -> Label { Label(q: q, b: b, constant: c) }
    let mixed = Fixture(
        kind: "shared and offset-private LCA blockers together", xBits: 2, wBits: wBits,
        labels: [fixed(0,1), fixed(1,2), fixed(2,2), fixed(2,3),
                 fixed(3,2), fixed(1,1), fixed(3,3,1)],
        leaves: [a,negA,b,negB,a,negA,b,negB]
    )
    let unreachable = Fixture(
        kind: "repeated shared parity produces a dead private path", xBits: 2, wBits: wBits,
        labels: [fixed(0,1),fixed(0,1),fixed(1,2),fixed(2,1),fixed(2,2),fixed(2,3),fixed(3,0)],
        leaves: [a,negA,b,negB,a,negA,b,b]
    )
    let dependency = Fixture(
        kind: "cross-branch reversal and dependent distinct private masks", xBits: 3, wBits: wBits,
        labels: [fixed(4,1),fixed(1,2),fixed(2,3),fixed(2,0),
                 fixed(3,3),fixed(1,1),fixed(3,2,1)],
        leaves: [a,negA,b,negB,b,negB,a,negA]
    )
    let sharedOnly = Fixture(
        kind: "all queries shared", xBits: 1, wBits: wBits,
        labels: [fixed(0,1),fixed(0,3),fixed(0,1)],
        leaves: [a,negA,b,negB]
    )
    return [mixed, unreachable, dependency, sharedOnly]
}

func randomFixture(_ xBits: Int, _ wBits: Int, _ depth: Int, _ rng: inout RNG) -> Fixture {
    let nodeCount = (1 << depth) - 1
    let path = paths(nodeCount)
    var all: [Label] = []
    for v in 0..<nodeCount {
        let old = path[v].nodes.map { all[$0].q }.filter { $0 != 0 }
        let mask = Int(rng.next() & UInt64((1 << wBits) - 1))
        let kind = rng.next() & 3
        if kind < 2 && old.count < xBits {
            let span = Set(frequencySpan(old, Array(repeating: 0, count: old.count)).map { $0.0 })
            let choices = (1..<(1 << xBits)).filter { !span.contains($0) }
            let form = choices[Int(rng.next() % UInt64(choices.count))]
            all.append(Label(q: form, b: mask, constant: Int(rng.next() & 1)))
        } else {
            let selector = kind == 3 ? 1 : mask
            all.append(Label(q: 0, b: selector, constant: 0))
        }
    }
    return Fixture(kind: "seeded path-independent repeated-shared parity tree",
                   xBits: xBits, wBits: wBits, labels: all,
                   leaves: leaves(nodeCount + 1, 1 << wBits, &rng))
}

func gcd(_ x: Int, _ y: Int) -> Int {
    var a = x
    var b = y
    while b != 0 { (a,b) = (b,a % b) }
    return a
}

struct Stats {
    var markers = 0
    var activePrivateSharedOffsetBlockers = 0
    var coefficientIsolation = 0
    var shiftChecks = 0
    var inactiveUnreachable = 0
    var markedBlockerVarying = 0
    var exactGateCapacityChecks = 0
    var densityRankChecks = 0
    var markedPairwiseDualChecks = 0
    var inactiveButFormalPrivate = 0
    var largestL = 0
    var strictRankImprovement = 0
}

func audit(_ f: Fixture, _ stats: inout Stats) throws -> ([Int32], [Int32], [String: Any]) {
    let xSize = 1 << f.xBits
    let wSize = 1 << f.wBits
    let labels = f.labels
    let path = paths(labels.count)
    let scale: Int32 = 32
    let maxDepth = path.last!.nodes.count
    let room = 1 << maxDepth
    let probabilityDenominator = (1...wSize).reduce(1) { a, b in a / gcd(a,b) * b }
    let privateAncestors = path.map { entry in entry.nodes.indices.filter { labels[entry.nodes[$0]].isPrivate } }
    var volumes = [Int](repeating: 0, count: labels.count)
    var means = [Int32](repeating: 0, count: wSize)
    var signal = labels.map { _ in [Int32](repeating: 0, count: wSize) }
    var truth = [Int32](repeating: 0, count: xSize*wSize)
    try check(f.leaves.allSatisfy { directWalsh($0).filter { $0 != 0 }.count > 1 },
              "Test leaves must be noncharacters")

    for w in 0..<wSize {
        var hits = [Int](repeating: 0, count: labels.count)
        for x in 0..<xSize {
            var v = 0
            while v < labels.count {
                hits[v] += 1
                v = 2*v + 1 + labels[v].bit(x,w)
            }
            truth[x + xSize*w] = f.leaves[v - labels.count][w]
        }
        var values = [Int32](repeating: 0, count: labels.count) + f.leaves.map { scale*$0[w] }
        for v in labels.indices.reversed() {
            let indices = path[v].nodes.indices.filter { !labels[path[v].nodes[$0]].isPrivate }
            let reachable = indices.allSatisfy { i in
                labels[path[v].nodes[i]].bit(0,w) == path[v].choices[i]
            }
            let k = privateAncestors[v].count
            try check(hits[v] == (reachable ? xSize >> k : 0), "Direct traversal and E_v disagree")
            if reachable { volumes[v] += 1 }
            let l = 2*v + 1
            let r = l + 1
            if labels[v].isPrivate {
                try check((values[l] + values[r]) % 2 == 0, "Haar scale insufficient")
                values[v] = (values[l] + values[r]) / 2
                if reachable { signal[v][w] = (values[l] - values[r]) / 2 }
            } else {
                values[v] = values[l + labels[v].bit(0,w)]
            }
        }
        means[w] = values[0]
        let direct = truth[(w*xSize)..<((w+1)*xSize)].reduce(Int32(0), +)
        try check(direct * scale == means[w] * Int32(xSize), "Root mean disagrees with truth table")
    }

    let root = directWalsh(means)
    let rootSparsity = root.filter { $0 != 0 }.count
    let small = signal.map(directWalsh)
    let smallS = small.map { $0.filter { $0 != 0 }.count }
    var formula = [Int32](repeating: 0, count: truth.count)
    var spectra: [[(Int,Int)]] = Array(repeating: [], count: labels.count)
    var best: [Int:ShiftAudit] = [:]
    for (beta,value) in root.enumerated() {
        let numerator = value * Int32(xSize)
        try check(numerator % scale == 0, "Root normalization failed")
        formula[beta * xSize] = numerator / scale
    }
    for v in labels.indices where labels[v].isPrivate {
        let activeParentIndices = privateAncestors[v]
        let pp = activeParentIndices.map { labels[path[v].nodes[$0]].q }
        let projection = frequencySpan(pp, Array(repeating: 0, count: pp.count))
        try check(Set(projection.map { $0.0 }).count == 1 << pp.count, "Dependent private ancestor")
        try check(!projection.contains { $0.0 == labels[v].q }, "Own private query depends on ancestors")
        let nodeShifts = try shifts(labels, v, path[v], xSize)
        stats.shiftChecks += nodeShifts.feasible
        if volumes[v] == 0 { stats.inactiveUnreachable += 1 }
        if smallS[v] == 0 { continue }
        best[v] = nodeShifts
        let rows = activeParentIndices.map { labels[path[v].nodes[$0]].jointMask(f.xBits) }
        let phases = activeParentIndices.map { index in
            path[v].choices[index] ^ labels[path[v].nodes[index]].constant
        }
        let own = labels[v]
        spectra[v] = frequencySpan(rows, phases).map { ($0.0 ^ own.jointMask(f.xBits), $0.1 ^ own.constant) }
        let liftScale = Int32(xSize >> pp.count)
        for (base,twist) in spectra[v] {
            for (beta,value) in small[v].enumerated() {
                let numerator = sign(twist) * liftScale * value
                try check(numerator % scale == 0, "Local parity offset normalization failed")
                formula[base ^ (beta << f.xBits)] += numerator / scale
            }
        }
    }

    var theoreticalCapacity = [Int](repeating: 0, count: xSize)
    var actualCapacity = theoreticalCapacity
    var globalBound = 0
    var globalRankBound = 0
    var rankMax = 0
    var L = 0
    var details: [[String: Any]] = []
    var subtreeTerm = [Int](repeating: 0, count: 2 * labels.count + 1)
    for v in best.keys.sorted() {
        let checked = best[v]!
        let ell = checked.count
        let k = privateAncestors[v].count
        let selected = path[v].nodes.indices.filter { (checked.chosen >> $0) & 1 == 1 }
        let blockingNodes = selected.map { path[v].nodes[$0] }
        let privateShifted = blockingNodes.filter { labels[$0].isPrivate && labels[$0].b != 0 }
        let shared = blockingNodes.filter { !labels[$0].isPrivate }
        if !privateShifted.isEmpty && !shared.isEmpty { stats.activePrivateSharedOffsetBlockers += 1 }
        if checked.largest > ell { stats.markedBlockerVarying += 1 }
        L = max(L, ell)
        stats.markers += 1
        let volume = volumes[v]
        try check(volume > 0 && smallS[v] * volume >= wSize, "Marked uncertainty bound failed")
        let countGate = (0..<truth.count).filter { point in
            let x = point % xSize
            let w = point / xSize
            return selected.allSatisfy { i in
                labels[path[v].nodes[i]].bit(x,w) == path[v].choices[i]
            }
        }.count
        try check((countGate << ell) >= truth.count, "Actual gate measure below dyadic charge")
        let gate = frequencySpan(
            blockingNodes.map { labels[$0].jointMask(f.xBits) },
            selected.map { i in path[v].choices[i] ^ labels[path[v].nodes[i]].constant }
        )
        try check(countGate > 0, "Active gate is infeasible")
        let quotient = truth.count / countGate
        let gateRank = quotient.trailingZeroBitCount
        try check(gateRank <= ell && countGate * (1 << gateRank) == truth.count,
                  "Gate cardinality does not give a dyadic row rank")
        var gateCoefficients: [Int:Int] = [:]
        for (frequency, twist) in gate {
            gateCoefficients[frequency, default: 0] += twist == 0 ? 1 : -1
        }
        try check(gateCoefficients.values.filter { $0 != 0 }.count == 1 << gateRank,
                  "Gate rank does not match exact gate spectrum")
        if gateRank < ell { stats.strictRankImprovement += 1 }
        rankMax = max(rankMax, gateRank)
        let liftScale = Int32(xSize >> k)
        var bandSize = 0
        for (base, twist) in spectra[v] {
            let alpha = base & (xSize - 1)
            theoreticalCapacity[alpha] += room >> ell
            actualCapacity[alpha] += countGate
            for (beta, value) in small[v].enumerated() {
                let index = base ^ (beta << f.xBits)
                if formula[index] != 0 { bandSize += 1 }
                let convolution = gate.reduce(into: Int32(0)) { $0 += formula[index ^ $1.0] * sign($1.1) }
                let expected = value * liftScale * sign(twist) * Int32(1 << ell)
                try check(convolution * scale == expected, "LCA gate convolution fails to isolate")
                stats.coefficientIsolation += 1
            }
        }
        try check((bandSize << ell) >= (smallS[v] << k), "True band too small")
        let markedTerm = wSize * (probabilityDenominator / volume)
        subtreeTerm[v] = markedTerm
        let rankCharge = (room * room) >> (2 * gateRank)
        let blockerCharge = (room * room) >> (2 * ell)
        try check(blockerCharge > 0 && rankCharge >= blockerCharge, "Scaled dyadic charge must be positive")
        globalBound += blockerCharge * (1 << k) * markedTerm
        globalRankBound += rankCharge * (1 << k) * markedTerm
        try check((bandSize << gateRank) >= (smallS[v] << k), "True band below refined gate-rank bound")
        if f.kind != "seeded path-independent repeated-shared parity tree" {
            details.append([
                "node": v, "q": labels[v].q, "data_offset_mask": labels[v].b,
                "private_ancestors": k, "shared_gate_cardinality": volume, "shared_data_domain": wSize,
                "delta": checked.delta, "minimal_blockers": ell, "blocking_nodes": blockingNodes,
                "shared_blocking_nodes": shared, "offset_private_blocking_nodes": privateShifted,
                "true_band_sparsity": bandSize, "marked_signal_sparsity": smallS[v],
                "actual_blocker_gate_cardinality": countGate, "blocker_gate_rank": gateRank,
            ])
        }
    }
    try check(theoreticalCapacity.allSatisfy { $0 <= room }, "Dyadic frequency charge exceeds one")
    try check(actualCapacity.allSatisfy { $0 <= truth.count }, "Actual gate measures exceed one")
    stats.exactGateCapacityChecks += xSize
    let s = formula.filter { $0 != 0 }.count
    let denominator = room * room * probabilityDenominator
    try check(best.isEmpty || globalBound > 0, "Nonempty markers gave zero weighted bound")
    try check(globalBound <= globalRankBound && globalRankBound <= (s - rootSparsity) * denominator,
              "Weighted sums exceed allowed non-root spectrum")

    var formal = [Int](repeating: 0, count: 2 * labels.count + 1)
    var splits = formal
    for v in labels.indices.reversed() {
        let left = 2*v + 1
        let right = left + 1
        formal[v] = formal[left] + formal[right]
        splits[v] = splits[left] + splits[right]
        if labels[v].isPrivate {
            if formal[v] > 0 || subtreeTerm[v] > 0 {
                if subtreeTerm[v] == 0 { stats.inactiveButFormalPrivate += 1 }
                formal[v] += 1
            }
            subtreeTerm[v] += 2 * (subtreeTerm[left] + subtreeTerm[right])
        } else {
            if formal[left] > 0 && formal[right] > 0 { splits[v] += 1 }
            subtreeTerm[v] = subtreeTerm[left] + subtreeTerm[right]
        }
        if formal[v] != 0 {
            let total = formal[v] + splits[v] + 1
            try check(4 * subtreeTerm[v] * volumes[v] >= total * total * wSize * probabilityDenominator,
                      "Density-weighted formal subtree rank bound fails")
            stats.densityRankChecks += 1
        }
    }
    let threshold = formal[0] + splits[0] + 1
    try check((threshold * threshold) <= s << (2 * L + 2), "Formal global rank+sparsity inequality fails")
    try check((threshold * threshold) <= s << (2 * rankMax + 2), "Refined formal gate-rank inequality fails")
    stats.largestL = max(stats.largestL, L)
    let info: [String:Any] = [
        "kind": f.kind, "x_bits": f.xBits, "w_bits": f.wBits,
        "sparsity": s, "root_mean_sparsity": rootSparsity,
        "weighted_lower_numerator": globalBound, "weighted_lower_denominator": denominator,
        "weighted_gate_rank_lower_numerator": globalRankBound, "max_gate_rank_among_marked": rankMax,
        "max_blocker_among_marked": L,
        "n_form": formal[0], "g": splits[0],
        "marked_node_details": details,
        "never_reached_private_nodes": labels.indices.filter { labels[$0].isPrivate && volumes[$0] == 0 },
    ]
    return (truth, formula, info)
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R46MixedParityBlockerMetal output.json")
    let seed: UInt64 = 0xa280460edabc2345
    var rng = RNG(state: seed)
    let metal = try MetalWalsh()
    var cases: [Fixture] = []
    for w in 2...4 {
        cases += fixedFixtures(w)
        for kind in 0..<64 {
            cases.append(randomFixture(3 + kind % 3, w, (kind & 1) == 0 ? 3 : 4, &rng))
        }
    }
    var results: [[String: Any]] = []
    var stats = Stats()
    var batches: [Int: [([Int32], [Int32])]] = [:]
    for f in cases {
        let (truth, prediction, info) = try audit(f, &stats)
        if f.kind != "seeded path-independent repeated-shared parity tree" { results.append(info) }
        if f.kind == "shared and offset-private LCA blockers together" {
            let details = info["marked_node_details"] as! [[String: Any]]
            let target = details.first { $0["node"] as? Int == 3 }
            let sharedNodes = target?["shared_blocking_nodes"] as? [Int]
            let privateNodes = target?["offset_private_blocking_nodes"] as? [Int]
            try check(sharedNodes == [0] && privateNodes == [1],
                      "Expected the same node to need both types of blocker")
        }
        if f.kind == "repeated shared parity produces a dead private path" {
            try check(info["never_reached_private_nodes"] as? [Int] == [4], "Expected a dead path")
        }
        batches[truth.count, default: []].append((truth, prediction))
    }
    var coefficientCount = 0
    var overlap: Double = 0
    for n in batches.keys.sorted() {
        let batch = batches[n]!
        let worker = CpuTask()
        let group = DispatchGroup()
        let input = batch.map { $0.0 }
        DispatchQueue.global(qos: .userInitiated).async(group: group) { worker.execute(input) }
        let gpu = try metal.transform(input.flatMap { $0 }, rowSize: n)
        group.wait()
        let (cpu, time) = worker.result()
        let gpuTime = metal.commandSpans.last!
        try check(cpu == gpu && cpu == batch.flatMap { $0.1 }, "CPU/Metal/Haar coefficient mismatch")
        coefficientCount += cpu.count
        overlap += max(0, min(time.1, gpuTime.1) - max(time.0, gpuTime.0))
    }
    let record: [String: Any] = [
        "schema": "a28-r046-mixed-parity-lca-blocker-metal-v1",
        "seed_hex": String(seed, radix: 16),
        "total_trees": cases.count,
        "random_trees": 192,
        "structured_trees": 12,
        "marked_nodes": stats.markers,
        "marked_nodes_with_shared_and_shifted_private_blockers": stats.activePrivateSharedOffsetBlockers,
        "unreachable_private_nodes": stats.inactiveUnreachable,
        "inactive_nodes_still_in_formal_private_tree": stats.inactiveButFormalPrivate,
        "max_l_among_marked": stats.largestL,
        "marked_nodes_with_gate_rank_strictly_smaller_than_blocker_count": stats.strictRankImprovement,
        "marked_nodes_with_shift_dependent_count": stats.markedBlockerVarying,
        "feasible_dual_shift_count": stats.shiftChecks,
        "gate_isolation_coefficient_checks": stats.coefficientIsolation,
        "both_exact_gate_capacity_checks": stats.exactGateCapacityChecks,
        "density_weighted_subtree_rank_checks": stats.densityRankChecks,
        "cpu_gpu_haar_coefficients": coefficientCount,
        "coefficient_mismatches": 0,
        "device": metal.device.name,
        "gpu_command_buffers": metal.dispatches,
        "gpu_command_seconds": metal.gpuSeconds,
        "cpu_gpu_measured_overlap_seconds": overlap,
        "all_checks_passed": true,
        "structured_results": results,
    ]
    let output = try JSONSerialization.data(withJSONObject: record, options: [.prettyPrinted, .sortedKeys])
    try (output + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R46 verified \(cases.count) mixed parity trees and \(coefficientCount) exact coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R46 verifier error: \(error)\n").utf8))
    exit(1)
}
