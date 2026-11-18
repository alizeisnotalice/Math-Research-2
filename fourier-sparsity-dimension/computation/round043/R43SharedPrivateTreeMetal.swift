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

enum Query {
    case privateBit(Int)
    case shared(parityMask: Int, quadraticA: Int, quadraticB: Int, invert: Int)

    func branch(address: Int, shared: Int) -> Int {
        switch self {
        case .privateBit(let bit):
            return (address >> bit) & 1
        case .shared(let mask, let a, let b, let invert):
            return parity(shared & mask) ^ (((shared >> a) & (shared >> b)) & 1) ^ invert
        }
    }

    var bit: Int? {
        if case .privateBit(let bit) = self { return bit }
        return nil
    }
}

struct Fixture {
    let kind: String
    let queries: [Query]
    let leaves: [[Int32]]
}

struct Reachability {
    let privateAncestors: Int
    let privateChoices: Int
    let shared: [(node: Int, choice: Int)]
}

func pathGeometry(_ queries: [Query]) -> [Reachability] {
    var paths = [Reachability(privateAncestors: 0, privateChoices: 0, shared: [])]
    for child in 1..<queries.count {
        let parent = (child - 1) / 2
        let right = ((child & 1) == 0 ? 1 : 0)
        let path = paths[parent]
        if let bit = queries[parent].bit {
            paths.append(Reachability(privateAncestors: path.privateAncestors | (1 << bit),
                                      privateChoices: path.privateChoices | (right << bit), shared: path.shared))
        } else {
            paths.append(Reachability(privateAncestors: path.privateAncestors, privateChoices: path.privateChoices,
                                      shared: path.shared + [(parent, right)]))
        }
    }
    return paths
}

func fixtures(_ trial: Int, _ dataBits: Int, _ rng: inout RNG) -> Fixture {
    let dataCount = 1 << dataBits
    // Shift the quadratic beyond the available data bits to obtain a pure shared parity.
    let pureParity = Query.shared(parityMask: 3, quadraticA: dataBits, quadraticB: dataBits, invert: 0)
    let nonlinear = Query.shared(parityMask: 1, quadraticA: 0, quadraticB: 1, invert: 0)
    let a = (0..<dataCount).map { sign(($0 & 3) == 3 ? 1 : 0) }
    let b = (0..<dataCount).map { sign(($0 & 3) == 2 ? 1 : 0) }
    let c = (0..<dataCount).map { sign(($0 & 3) == 1 ? 1 : 0) }
    let d = (0..<dataCount).map { sign(($0 & 3) == 0 ? 1 : 0) }

    if trial == 0 {
        // The left child repeats the root parity, hence its right branch is impossible.
        return Fixture(kind: "repeated parity plus nonlinear shared query",
                       queries: [pureParity, pureParity, nonlinear, .privateBit(0), .privateBit(1),
                                 .privateBit(2), .privateBit(3)],
                       leaves: [a, b, b, a, c, d, d, c])
    }
    if trial == 1 {
        // Root mean difference vanishes, but the private descendants keep its bit relevant.
        return Fixture(kind: "zero private root difference with active descendants",
                       queries: [.privateBit(0), pureParity, pureParity, .privateBit(1),
                                 .privateBit(2), .privateBit(3), .privateBit(4)],
                       leaves: [a, b, a, b, b, a, b, a])
    }
    if trial == 2 {
        let leaf = (0..<dataCount).map { sign($0 & 1) }
        return Fixture(kind: "exact n-bound witness",
                       queries: (0..<7).map(Query.privateBit),
                       leaves: (0..<8).map { $0 & 1 == 0 ? leaf : leaf.map { -$0 } })
    }
    if trial == 3 {
        return Fixture(kind: "no private queries",
                       queries: [pureParity, pureParity, nonlinear, nonlinear, pureParity, nonlinear, pureParity],
                       leaves: [d, a, a, c, b, d, b, c])
    }
    if trial == 4 {
        let plus = [Int32](repeating: 1, count: dataCount)
        return Fixture(kind: "shared root multiplexes two independent private queries",
                       queries: [pureParity, .privateBit(0), .privateBit(1)],
                       leaves: [plus, plus.map { -$0 }, plus, plus.map { -$0 }])
    }

    let depth = (trial & 1 == 0 ? 4 : 3)
    let leafCount = 1 << depth
    var privateCount = 0
    var queries: [Query] = []
    for i in 0..<(leafCount - 1) {
        let decision = rng.next()
        if decision & 3 < 2 && privateCount < 8 {
            queries.append(.privateBit(privateCount))
            privateCount += 1
        } else if decision & 3 == 2 {
            queries.append(pureParity)
        } else {
            let mask = Int(rng.next() & UInt64(dataCount - 1))
            let a = i % dataBits
            let b = (a + 1) % dataBits
            queries.append(.shared(parityMask: mask, quadraticA: a, quadraticB: b,
                                   invert: Int(rng.next() & 1)))
        }
    }
    let leaves = (0..<leafCount).map { _ in
        (0..<dataCount).map { _ in sign(Int(rng.next() & 1)) }
    }
    return Fixture(kind: "seeded random mix", queries: queries, leaves: leaves)
}

struct Audit {
    let truth: [Int32]
    let spectrum: [Int32]
    let s: Int
    let rootS: Int
    let localS: [Int]
    let c: Int
    let effective: Int
    let effectiveSharedSplits: Int
    let privateCount: Int
    let unreachablePrivate: [Int]
    let inactiveSignalEffectivePrivate: [Int]
    let unreachablePrivateSharedPairs: Int
    let activeSubtreeChecks: Int
}

func audit(_ fixture: Fixture, dataBits: Int) throws -> Audit {
    let dataCount = 1 << dataBits
    let q = fixture.queries
    let paths = pathGeometry(q)
    let privateCount = q.compactMap { $0.bit }.count
    let addressCount = 1 << privateCount
    let n = dataCount * addressCount
    let addressMask = addressCount - 1
    let scale: Int32 = 16
    var truth = [Int32](repeating: 0, count: n)
    var means = [Int32](repeating: 0, count: dataCount)
    var signals = q.map { _ in [Int32](repeating: 0, count: dataCount) }
    var sizes = [Int](repeating: 0, count: q.count)
    var unreachable = 0

    for w in 0..<dataCount {
        var numbers = [Int](repeating: 0, count: q.count)
        for x in 0..<addressCount {
            var v = 0
            while v < q.count {
                numbers[v] += 1
                let decision = q[v].branch(address: x, shared: w)
                v = 2 * v + 1 + decision
            }
            truth[x | (w << privateCount)] = fixture.leaves[v - q.count][w]
        }
        var values = [Int32](repeating: 0, count: q.count) + fixture.leaves.map { scale * $0[w] }
        for v in q.indices.reversed() {
            let left = 2 * v + 1
            let right = left + 1
            let path = paths[v]
            let reachable = path.shared.allSatisfy { q[$0.node].branch(address: 0, shared: w) == $0.choice }
            if reachable { sizes[v] += 1 }
            let expected = reachable ? addressCount >> path.privateAncestors.nonzeroBitCount : 0
            try check(numbers[v] == expected, "Collapsed shared reachability count mismatch")
            if q[v].bit != nil {
                try check((values[left] + values[right]) % 2 == 0, "Working scale insufficient")
                values[v] = (values[left] + values[right]) / 2
                if reachable { signals[v][w] = (values[left] - values[right]) / 2 }
                else { unreachable += 1 }
            } else {
                let branch = q[v].branch(address: 0, shared: w)
                values[v] = values[left + branch]
            }
        }
        means[w] = values[0]
        let directMean = truth[(w * addressCount)..<((w + 1) * addressCount)].reduce(Int32(0), +)
        try check(directMean * scale == means[w] * Int32(addressCount), "Collapsed private mean is wrong")
    }

    let root = directWalsh(means)
    let rootS = root.filter { $0 != 0 }.count
    var predicted = [Int32](repeating: 0, count: n)
    for eta in 0..<dataCount {
        let unscaled = root[eta] * Int32(addressCount)
        try check(unscaled % scale == 0, "Mean lifting nonintegral")
        predicted[eta << privateCount] = unscaled / scale
    }
    var localS = [Int](repeating: 0, count: q.count)
    var closure = 0
    for v in q.indices {
        guard let bit = q[v].bit else { continue }
        let path = paths[v]
        let spectral = directWalsh(signals[v])
        localS[v] = spectral.filter { $0 != 0 }.count
        if localS[v] == 0 { continue }
        closure |= path.privateAncestors | (1 << bit)
        var mask = path.privateAncestors
        let addressFactor = Int32(addressCount >> path.privateAncestors.nonzeroBitCount)
        while true {
            let xFreq = (1 << bit) | mask
            for eta in 0..<dataCount {
                let position = xFreq | (eta << privateCount)
                try check(predicted[position] == 0, "Private Haar components overlap")
                let numerator = spectral[eta] * addressFactor * sign(parity(mask & path.privateChoices))
                try check(numerator % scale == 0, "Private lifting nonintegral")
                predicted[position] = numerator / scale
            }
            if mask == 0 { break }
            mask = (mask - 1) & path.privateAncestors
        }
    }

    let spectrum = cpuFWHT(truth)
    try check(predicted == spectrum, "Haar formula differs from full CPU spectrum")
    let support = spectrum.indices.filter { spectrum[$0] != 0 }
    let effective = support.reduce(0) { $0 | ($1 & addressMask) }
    try check(effective == closure, "Fourier effective private variables and closure differ")
    for b in 0..<privateCount {
        let sensitive = truth.indices.contains { truth[$0] != truth[$0 ^ (1 << b)] }
        try check(sensitive == (effective & (1 << b) != 0), "Truth-table relevance mismatch")
    }

    var c = [Int](repeating: 0, count: 2 * q.count + 1)
    var nEff = c
    var gEff = c
    var subtreeChecks = 0
    for v in q.indices.reversed() {
        let left = 2 * v + 1
        let right = left + 1
        let children = nEff[left] + nEff[right]
        gEff[v] = gEff[left] + gEff[right]
        if q[v].bit != nil {
            c[v] = localS[v] + 2 * (c[left] + c[right])
            nEff[v] = children + (children > 0 || localS[v] > 0 ? 1 : 0)
        } else {
            c[v] = c[left] + c[right]
            nEff[v] = children
            if nEff[left] > 0 && nEff[right] > 0 { gEff[v] += 1 }
        }
        if nEff[v] > 0 {
            let plusOne = nEff[v] + gEff[v] + 1
            try check(4 * c[v] * sizes[v] >= plusOne * plusOne * dataCount,
                      "Gated subtree Cauchy-uncertainty bound fails")
            subtreeChecks += 1
        }
    }
    let unionSparsity = support.count
    try check(unionSparsity == rootS + c[0] && nEff[0] == effective.nonzeroBitCount,
              "Global sparsity or dimension formula failed")
    let rankBudget = nEff[0] + gEff[0] + 1
    try check(rankBudget * rankBudget <= 4 * unionSparsity, "Effective private/shared bound fails")
    let neverReached = q.indices.filter { q[$0].bit != nil && sizes[$0] == 0 }
    let noSignalButEffective = q.indices.filter { node in
        guard let bit = q[node].bit else { return false }
        return localS[node] == 0 && effective & (1 << bit) != 0
    }
    return Audit(truth: truth, spectrum: spectrum, s: unionSparsity, rootS: rootS, localS: localS, c: c[0],
                 effective: nEff[0], effectiveSharedSplits: gEff[0], privateCount: privateCount,
                 unreachablePrivate: neverReached,
                 inactiveSignalEffectivePrivate: noSignalButEffective,
                 unreachablePrivateSharedPairs: unreachable, activeSubtreeChecks: subtreeChecks)
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R43SharedPrivateTreeMetal output.json")
    let metal = try MetalWalsh()
    let seed: UInt64 = 0xa280430d6ab12345
    var rng = RNG(state: seed)
    var totalCases = 0
    var zeroSignalEffective = 0
    var impossiblePrivateNodes = 0
    var impossiblePrivateSharedPairs = 0
    var activeSubtreeChecks = 0
    var equalityCases = 0
    var sharedSplitTotal = 0
    var positiveSharedSplitCases = 0
    var byDomain: [Int: [(truth: [Int32], spectrum: [Int32])]] = [:]
    var series: [[String: Any]] = []
    var examples: [[String: Any]] = []
    for dataBits in 2...5 {
        var localCoefficients = 0
        for trial in 0..<32 {
            let fixture = fixtures(trial, dataBits, &rng)
            let result = try audit(fixture, dataBits: dataBits)
            if trial == 0 {
                try check(result.unreachablePrivate == [4], "Repeated-parity dead path not found")
            }
            if trial == 1 {
                try check(result.inactiveSignalEffectivePrivate == [0] && result.effective == 5,
                          "Root should remain an effective variable when its difference vanishes")
            }
            if trial == 2 {
                try check(result.effective == 7 && result.s == 16, "Sharp fixture differs")
            }
            if trial == 4 {
                try check(result.effective == 2 && result.effectiveSharedSplits == 1 && result.s == 4,
                          "Mixed sharp fixture differs")
            }
            if trial < 5 {
                examples.append([
                    "kind": fixture.kind, "data_bits": dataBits, "allocated_private_queries": result.privateCount,
                    "effective_private_queries": result.effective, "sparsity": result.s,
                    "shared_nodes_with_two_effective_private_branches": result.effectiveSharedSplits,
                    "root_mean_sparsity": result.rootS, "private_haar_cost": result.c,
                    "signal_sparsity_by_node": result.localS, "never_reached_private_nodes": result.unreachablePrivate,
                    "also_effective_private_nodes_with_zero_own_signal": result.inactiveSignalEffectivePrivate,
                ])
            }
            byDomain[result.truth.count, default: []].append((result.truth, result.spectrum))
            totalCases += 1
            localCoefficients += result.truth.count
            zeroSignalEffective += result.inactiveSignalEffectivePrivate.count
            impossiblePrivateNodes += result.unreachablePrivate.count
            impossiblePrivateSharedPairs += result.unreachablePrivateSharedPairs
            activeSubtreeChecks += result.activeSubtreeChecks
            sharedSplitTotal += result.effectiveSharedSplits
            if result.effectiveSharedSplits > 0 { positiveSharedSplitCases += 1 }
            let rankBudget = result.effective + result.effectiveSharedSplits + 1
            if rankBudget * rankBudget == 4 * result.s { equalityCases += 1 }
        }
        series.append(["data_bits": dataBits, "cases": 32, "coefficients": localCoefficients])
    }

    var coefficients = 0
    for size in byDomain.keys.sorted() {
        let batch = byDomain[size]!
        let gpu = try metal.transform(batch.flatMap { $0.truth }, rowSize: size)
        try check(gpu == batch.flatMap { $0.spectrum }, "Metal differs from independently checked CPU spectrum")
        coefficients += gpu.count
    }
    let record: [String: Any] = [
        "schema": "a28-r043-shared-private-haar-tree-metal-v1",
        "seed_hex": String(seed, radix: 16),
        "shared_data_range": [2, 5],
        "private_coordinate_limit": 8,
        "cases": totalCases,
        "equality_cases_for_effective_private_shared_bound": equalityCases,
        "effective_shared_split_total": sharedSplitTotal,
        "cases_with_effective_shared_split": positiveSharedSplitCases,
        "total_never_reached_private_node_instances": impossiblePrivateNodes,
        "number_of_private_node_shared_input_pairs_unreachable": impossiblePrivateSharedPairs,
        "effective_private_node_instances_with_zero_own_signal": zeroSignalEffective,
        "active_subtree_gated_lower_bound_checks": activeSubtreeChecks,
        "coefficients_cpu_gpu_and_haar": coefficients,
        "cpu_gpu_and_haar_disagreements": 0,
        "cpu_backend": "exact Int32 FWHT with direct integer intrinsic Walsh coefficients",
        "gpu_backend": "exact Int32 staged Metal FWHT",
        "device": metal.device.name,
        "metal4_supported": metal.device.supportsFamily(.metal4),
        "command_buffers_completed": metal.dispatches,
        "gpu_command_seconds": metal.gpuSeconds,
        "series": series,
        "examples": examples,
        "all_checks_passed": true,
    ]
    let data = try JSONSerialization.data(withJSONObject: record, options: [.prettyPrinted, .sortedKeys])
    try (data + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R43 verified \(totalCases) mixed trees and \(coefficients) exact coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R43 verifier error: \(error)\n").utf8))
    exit(1)
}
