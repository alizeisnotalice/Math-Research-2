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

indirect enum Tree {
    case leaf
    case query(Int, Tree, Tree)
}

func orderedTrees(minimum: Int, labels: Int) -> [Tree] {
    var result: [Tree] = [.leaf]
    for label in minimum..<labels {
        let smaller = orderedTrees(minimum: label + 1, labels: labels)
        for left in smaller {
            for right in smaller { result.append(.query(label, left, right)) }
        }
    }
    return result
}

struct Node {
    let label: Int?
    let pathMask: Int
    let pathChoices: Int
    let left: Int
    let right: Int
    let leafIndex: Int?
}

struct FlatTree {
    let nodes: [Node]
    let selector: [Int]
    let leafCount: Int
}

func flattened(_ tree: Tree, _ labels: Int) -> FlatTree {
    var nodes: [Node] = []
    var leafCount = 0
    func visit(_ subtree: Tree, _ mask: Int, _ choices: Int) -> Int {
        let index = nodes.count
        nodes.append(Node(label: nil, pathMask: mask, pathChoices: choices, left: 0, right: 0, leafIndex: nil))
        switch subtree {
        case .leaf:
            nodes[index] = Node(label: nil, pathMask: mask, pathChoices: choices,
                                left: 0, right: 0, leafIndex: leafCount)
            leafCount += 1
        case .query(let j, let left, let right):
            let bit = 1 << j
            let l = visit(left, mask | bit, choices)
            let r = visit(right, mask | bit, choices | bit)
            nodes[index] = Node(label: j, pathMask: mask, pathChoices: choices,
                                left: l, right: r, leafIndex: nil)
        }
        return index
    }
    _ = visit(tree, 0, 0)
    let selector = (0..<(1 << labels)).map { address in
        var v = 0
        while let j = nodes[v].label {
            v = (address & (1 << j)) == 0 ? nodes[v].left : nodes[v].right
        }
        return nodes[v].leafIndex!
    }
    return FlatTree(nodes: nodes, selector: selector, leafCount: leafCount)
}

func leafTables(
    _ flat: FlatTree, _ dataSize: Int, _ labels: Int, _ mode: Int, _ rng: inout RNG
) -> [[Int32]] {
    let common = (0..<dataSize).map { sign(($0 & 3) == 3 ? 1 : 0) }
    if mode == 0 {
        return Array(repeating: common, count: flat.leafCount)
    }
    if mode == 3 {
        return (0..<flat.leafCount).map { _ in
            while true {
                let values = (0..<dataSize).map { _ in sign(Int(rng.next() & 1)) }
                if directWalsh(values).filter({ $0 != 0 }).count > 1 { return values }
            }
        }
    }
    let nodesByLabel = (0..<labels).map { j in flat.nodes.indices.filter { flat.nodes[$0].label == j } }
    let target = nodesByLabel.indices.max { nodesByLabel[$0].count < nodesByLabel[$1].count }!
    let targetBits = 1 << target
    let picked = nodesByLabel[target].enumerated().compactMap { i, v in (mode == 1 || i % 2 == 0) ? v : nil }
    var leaves = Array(repeating: common, count: flat.leafCount)
    for node in flat.nodes {
        guard let leaf = node.leafIndex, node.pathMask & targetBits != 0 else { continue }
        let owner = picked.contains { i in
            node.pathMask & flat.nodes[i].pathMask == flat.nodes[i].pathMask
                && node.pathChoices & flat.nodes[i].pathMask == flat.nodes[i].pathChoices
        }
        if owner && node.pathChoices & targetBits != 0 { leaves[leaf] = common.map { -$0 } }
    }
    return leaves
}

struct BandAudit {
    let label: Int
    let physicalNodes: Int
    let activeNodes: Int
    let rawOccurrenceCount: Int
    let preUnion: Int
    let afterCollision: Int
    let liftedMax: Int
    let isolatorSizes: [Int]
}

struct TreeAudit {
    let truth: [Int32]
    let formula: [Int32]
    let bands: [BandAudit]
    let shapeLeaves: Int
}

func smallestIsolator(
    _ v: Int, _ active: [Int], _ nodes: [Node]
) throws -> Int {
    let own = nodes[v]
    let conflicts: [Int] = active.filter { $0 != v }.map { other in
        let rival = nodes[other]
        return own.pathMask & rival.pathMask & (own.pathChoices ^ rival.pathChoices)
    }
    try check(conflicts.allSatisfy { $0 != 0 }, "Repeated label nodes have overlapping query paths")
    var best = own.pathMask
    var subset = own.pathMask
    while true {
        if subset.nonzeroBitCount < best.nonzeroBitCount
            && conflicts.allSatisfy({ subset & $0 != 0 }) { best = subset }
        if subset == 0 { break }
        subset = (subset - 1) & own.pathMask
    }
    return best
}

func audit(
    _ tree: FlatTree, _ leaves: [[Int32]], _ labels: Int, _ dataSize: Int
) throws -> TreeAudit {
    let n = dataSize << labels
    let addressCount = 1 << labels
    let nodes = tree.nodes
    let truth = (0..<n).map { i in leaves[tree.selector[i % addressCount]][i >> labels] }
    try check(leaves.allSatisfy { directWalsh($0).filter { $0 != 0 }.count > 1 },
              "Expected each leaf be noncharacter")
    let scale: Int32 = 16
    var signals = nodes.map { _ in [Int32](repeating: 0, count: dataSize) }
    var root = [Int32](repeating: 0, count: dataSize)
    for y in 0..<dataSize {
        var means = [Int32](repeating: 0, count: nodes.count)
        for (i, node) in nodes.enumerated().reversed() {
            if let leaf = node.leafIndex { means[i] = scale * leaves[leaf][y] }
            else {
                try check((means[node.left] + means[node.right]) % 2 == 0, "Subtree scale insufficient")
                means[i] = (means[node.left] + means[node.right]) / 2
                signals[i][y] = (means[node.left] - means[node.right]) / 2
            }
        }
        root[y] = means[0]
        let start = y << labels
        let direct = truth[start..<(start + addressCount)].reduce(Int32(0), +)
        try check(root[y] * Int32(addressCount) == direct * scale, "Root mean not uniform")
    }
    let transforms = signals.map(directWalsh)
    let sizes = transforms.map { $0.filter { $0 != 0 }.count }
    var predicted = [Int32](repeating: 0, count: n)
    var occurrences = [Int](repeating: 0, count: n)
    for (eta, value) in directWalsh(root).enumerated() {
        let numerator = value * Int32(addressCount)
        try check(numerator % scale == 0, "Root lift not integral")
        predicted[eta << labels] = numerator / scale
    }
    for (i,node) in nodes.enumerated() {
        guard let j = node.label, sizes[i] != 0 else { continue }
        try check(node.pathMask < 1 << j, "Tree not ordered on a path")
        let depth = node.pathMask.nonzeroBitCount
        let weight = Int32(addressCount >> depth)
        var subset = node.pathMask
        while true {
            let frequency = subset | (1 << j)
            let phase = sign(parity(subset & node.pathChoices))
            for (eta, value) in transforms[i].enumerated() where value != 0 {
                let index = frequency | (eta << labels)
                let numerator = value * phase * weight
                try check(numerator % scale == 0, "Signal lift not integral")
                predicted[index] += numerator / scale
                occurrences[index] += 1
            }
            if subset == 0 { break }
            subset = (subset - 1) & node.pathMask
        }
    }

    var bands: [BandAudit] = []
    for j in 0..<labels {
        let physical = nodes.indices.filter { nodes[$0].label == j }
        if physical.isEmpty { continue }
        let active = physical.filter { sizes[$0] != 0 }
        let positions = (0..<dataSize).flatMap { eta in
            ((1 << j)..<(2 << j)).map { (eta << labels) | $0 }
        }
        let s = positions.filter { predicted[$0] != 0 }.count
        let before = positions.reduce(0) { value, idx in value + occurrences[idx] }
        let beforeUnion = positions.filter { occurrences[$0] > 0 }.count
        try check(active.isEmpty == (s == 0), "All signals vanished in a nonempty band")
        var isolatorSizes: [Int] = []
        var maxLift = 0
        for v in active {
            let node = nodes[v]
            let liftSize = (1 << node.pathMask.nonzeroBitCount) * sizes[v]
            let choose = try smallestIsolator(v, active, nodes)
            let distinct = choose.nonzeroBitCount
            try check(distinct <= active.count - 1, "Minimal isolator exceeds band multiplicity")
            try check((s << distinct) >= liftSize, "Refined sparsity bound fails")
            try check((s << (active.count - 1)) >= liftSize, "Band multiplicity bound fails")
            for x in 0..<(1 << j) where x & choose == node.pathChoices & choose {
                for w in 0..<dataSize {
                    let selected = active.reduce(into: Int32(0)) { value, other in
                        if x & nodes[other].pathMask == nodes[other].pathChoices { value += signals[other][w] }
                    }
                    let isolated = x & node.pathMask == node.pathChoices ? signals[v][w] : 0
                    try check(selected == isolated, "Isolation does not recover one complete signal")
                }
            }
            isolatorSizes.append(distinct)
            maxLift = max(maxLift, liftSize)
        }
        bands.append(BandAudit(label: j, physicalNodes: physical.count, activeNodes: active.count,
                               rawOccurrenceCount: before, preUnion: beforeUnion, afterCollision: s,
                               liftedMax: maxLift, isolatorSizes: isolatorSizes))
    }
    return TreeAudit(truth: truth, formula: predicted, bands: bands, shapeLeaves: tree.leafCount)
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

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R44OrderedCollisionMetal output.json")
    let labels = 4
    let shapes = orderedTrees(minimum: 0, labels: labels).map { flattened($0, labels) }
    try check(shapes.count == 1806, "Wrong number of ordered shapes")
    let seed: UInt64 = 0xa280440dcaf12345
    var rng = RNG(state: seed)
    let metal = try MetalWalsh()
    var coefficients = 0
    var cases = 0
    var multiplicityChecks = 0
    var bandsWithRepeatedActiveNodes = 0
    var bandsWithCancellations = 0
    var simultaneousInactiveAndRepeated = 0
    var boundsTightWithRepeatedNodes = 0
    var inactiveNodes = 0
    var shapeExamples: [String: [String: Any]] = [:]
    var series: [[String: Any]] = []
    let modes = 4
    for dataBits in 2...3 {
        let dataSize = 1 << dataBits
        var cpuInputs: [[Int32]] = []
        var formula: [Int32] = []
        for (shapeId, tree) in shapes.enumerated() {
            for mode in 0..<modes {
                let leaves = leafTables(tree, dataSize, labels, mode, &rng)
                let checked = try audit(tree, leaves, labels, dataSize)
                if mode == 0 {
                    try check(checked.bands.allSatisfy { $0.activeNodes == 0 }, "Constant label mode did not deactivate queries")
                }
                for band in checked.bands {
                    multiplicityChecks += band.activeNodes
                    inactiveNodes += band.physicalNodes - band.activeNodes
                    if band.activeNodes > 1 { bandsWithRepeatedActiveNodes += 1 }
                    if band.preUnion > band.afterCollision { bandsWithCancellations += 1 }
                    if band.activeNodes > 1 && band.physicalNodes > band.activeNodes { simultaneousInactiveAndRepeated += 1 }
                    let tight = band.activeNodes > 1 && band.liftedMax == band.afterCollision << (band.activeNodes - 1)
                    if tight { boundsTightWithRepeatedNodes += 1 }
                    let key: String?
                    if tight && band.activeNodes == 3 && band.preUnion > band.afterCollision { key = "three-active-tight" }
                    else if band.activeNodes > 1 && band.physicalNodes > band.activeNodes { key = "repeated-with-inactive" }
                    else { key = nil }
                    if let key, shapeExamples[key] == nil {
                        shapeExamples[key] = [
                            "data_bits": dataBits, "shape_index": shapeId, "leaf_mode": mode,
                            "maximum_frequency_label": band.label, "physical_query_sites": band.physicalNodes,
                            "active_query_sites": band.activeNodes, "number_of_leaves": checked.shapeLeaves,
                            "precollision_occurrences": band.rawOccurrenceCount, "precollision_union": band.preUnion,
                            "postcollision_band_sparsity": band.afterCollision, "largest_local_lift": band.liftedMax,
                            "minimum_isolator_sizes": band.isolatorSizes,
                            "nodes": tree.nodes.map { node in
                                ["label": node.label ?? -1, "path_mask": node.pathMask, "path_values": node.pathChoices,
                                 "left": node.left, "right": node.right, "leaf_index": node.leafIndex ?? -1]
                            },
                            "leaf_truth_tables": leaves,
                        ]
                    }
                }
                cases += 1
                cpuInputs.append(checked.truth)
                formula.append(contentsOf: checked.formula)
            }
        }
        let cpuWorker = CpuTask()
        let group = DispatchGroup()
        let batch = cpuInputs
        DispatchQueue.global(qos: .userInitiated).async(group: group) { cpuWorker.execute(batch) }
        let gpu = try metal.transform(batch.flatMap { $0 }, rowSize: dataSize << labels)
        group.wait()
        let (cpu, span) = cpuWorker.result()
        let gpuSpan = metal.commandSpans.last!
        try check(cpu == gpu && cpu == formula, "Parallel CPU / Metal / local Haar formula disagree")
        coefficients += cpu.count
        series.append([
            "data_bits": dataBits, "shapes": shapes.count, "modes": modes, "cases": cpuInputs.count,
            "coefficients": cpu.count,
            "cpu_fwht_wall_seconds": span.1 - span.0,
            "gpu_fwht_command_seconds": gpuSpan.1 - gpuSpan.0,
            "measured_cpu_gpu_overlap_seconds": max(0, min(span.1, gpuSpan.1) - max(span.0, gpuSpan.0)),
        ])
    }
    try check(shapeExamples["three-active-tight"] != nil, "No exact exponent-two collision fixture")
    try check(shapeExamples["repeated-with-inactive"] != nil, "No simultaneous repeated and inactive nodes")
    let result: [String: Any] = [
        "schema": "a28-r044-ordered-tree-collision-bands-metal-v1",
        "ordered_shapes_on_four_query_labels": shapes.count,
        "cases": cases,
        "data_variable_counts": [2, 3],
        "leaf_functions": "only noncharacters",
        "seed_hex": String(seed, radix: 16),
        "batches": series,
        "coefficient_checks_cpu_metal_formula": coefficients,
        "active_site_multiplicity_bound_checks": multiplicityChecks,
        "bands_with_repeated_active_label": bandsWithRepeatedActiveNodes,
        "bands_with_complete_coefficient_cancellation": bandsWithCancellations,
        "bands_with_two_or_more_active_and_some_inactive_sites": simultaneousInactiveAndRepeated,
        "bounds_tight_bands_with_two_or_more_active_sites": boundsTightWithRepeatedNodes,
        "individual_inactive_node_occurrences": inactiveNodes,
        "examples": shapeExamples,
        "command_buffers_completed": metal.dispatches,
        "gpu_device": metal.device.name,
        "gpu_backend": "Metal integer FWHT",
        "cpu_backend": "Int32 CPU FWHT with small integer direct Walsh and branch analysis",
        "cpu_gpu_formula_mismatches": 0,
        "all_checks_passed": true,
    ]
    let bytes = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try (bytes + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R44 verified \(cases) ordered tree cases and \(coefficients) integer coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R44 verifier error: \(error)\n").utf8))
    exit(1)
}
