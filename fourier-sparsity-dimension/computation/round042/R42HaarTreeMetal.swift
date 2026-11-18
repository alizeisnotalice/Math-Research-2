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

struct Node {
    let depth: Int
    let ancestors: Int
    let choices: Int
    let firstLeaf: Int
    let leafCount: Int
}

func nodeGeometry() -> [Node] {
    var nodes = [Node(depth: 0, ancestors: 0, choices: 0, firstLeaf: 0, leafCount: 8)]
    for index in 1..<7 {
        let parentIndex = (index - 1) / 2
        let parent = nodes[parentIndex]
        let isRight = (index & 1) == 0
        let halfSize = parent.leafCount / 2
        nodes.append(Node(
            depth: parent.depth + 1,
            ancestors: parent.ancestors | (1 << parentIndex),
            choices: parent.choices | (isRight ? 1 << parentIndex : 0),
            firstLeaf: parent.firstLeaf + (isRight ? halfSize : 0),
            leafCount: halfSize
        ))
    }
    return nodes
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

func leafTables(trial: Int, dataCount: Int, rng: inout RNG) -> [[Int32]] {
    // Translated quadratics have four Walsh coefficients, independently of extra data bits.
    let a = (0..<dataCount).map { sign(($0 & 3) == 3 ? 1 : 0) }
    let b = (0..<dataCount).map { sign(($0 & 3) == 2 ? 1 : 0) }
    let c = (0..<dataCount).map { sign(($0 & 3) == 1 ? 1 : 0) }
    let d = (0..<dataCount).map { sign(($0 & 3) == 0 ? 1 : 0) }
    switch trial {
    case 0:
        return [a, b, b, a, c, d, d, c]
    case 1:
        return [a, b, a, c, b, a, c, a]
    case 2:
        return [a, a.map { -$0 }, b, b.map { -$0 }, c, c.map { -$0 }, d, d.map { -$0 }]
    default:
        return (0..<8).map { _ in
            while true {
                let truth = (0..<dataCount).map { _ in sign(Int(rng.next() & 1)) }
                if directWalsh(truth).filter({ $0 != 0 }).count > 1 { return truth }
            }
        }
    }
}

func leafIndex(for address: Int) -> Int {
    var node = 0
    for _ in 0..<3 {
        node = 2 * node + 1 + ((address >> node) & 1)
    }
    return node - 7
}

func liftSpectrum(
    _ leaves: [[Int32]], _ nodes: [Node], _ dataCount: Int
) throws -> ([Int32], Int, [Int], [Int], [Int]) {
    let domainSize = dataCount << 7
    // Store eight times each subtree mean and eight times each Haar difference.
    let rootMean = directWalsh((0..<dataCount).map { y in
        leaves.reduce(into: Int32(0)) { $0 += $1[y] }
    })
    let rootSparsity = rootMean.filter { $0 != 0 }.count
    var result = [Int32](repeating: 0, count: domainSize)
    for eta in 0..<dataCount { result[eta << 7] = 16 * rootMean[eta] }
    var deltaCounts: [Int] = []
    var zeroNonconstant: [Int] = []
    var zeroWithNonconstantCommonMean: [Int] = []

    for (index, node) in nodes.enumerated() {
        let middle = node.firstLeaf + node.leafCount / 2
        let end = node.firstLeaf + node.leafCount
        let left = node.firstLeaf..<middle
        let right = middle..<end
        let leftSum = (0..<dataCount).map { y in
            left.reduce(into: Int32(0)) { value, leaf in value += leaves[leaf][y] }
        }
        let rightSum = (0..<dataCount).map { y in
            right.reduce(into: Int32(0)) { value, leaf in value += leaves[leaf][y] }
        }
        let scale = Int32(1 << node.depth)
        let delta = zip(leftSum, rightSum).map { scale * ($0 - $1) }
        let intrinsic = directWalsh(delta)
        let count = intrinsic.filter { $0 != 0 }.count
        deltaCounts.append(count)
        if count == 0 && Set(leaves[node.firstLeaf..<end]).count > 1 {
            zeroNonconstant.append(index)
            if Set(leftSum).count > 1 { zeroWithNonconstantCommonMean.append(index) }
        }

        let dataFactor = Int32(1 << (4 - node.depth))
        var ancestorFrequency = node.ancestors
        while true {
            let freq = ancestorFrequency | (1 << index)
            let phaseSign = sign(parity(ancestorFrequency & node.choices))
            for eta in 0..<dataCount {
                let position = freq | (eta << 7)
                try check(result[position] == 0, "Two Haar spectra overlapped")
                result[position] = dataFactor * phaseSign * intrinsic[eta]
            }
            if ancestorFrequency == 0 { break }
            ancestorFrequency = (ancestorFrequency - 1) & node.ancestors
        }
    }
    return (result, rootSparsity, deltaCounts, zeroNonconstant, zeroWithNonconstantCommonMean)
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 2, "usage: R42HaarTreeMetal output.json")
    let metal = try MetalWalsh()
    let nodes = nodeGeometry()
    let seed: UInt64 = 0xA280420BCDE91234
    var rng = RNG(state: seed)
    var cases = 0
    var results: [[String: Any]] = []
    var examples: [[String: Any]] = []
    var countCoefficients = 0
    var noncharacterLeafChecks = 0
    var zeroNonconstantCount = 0
    var zeroNonconstantCommonMeanCount = 0
    let addresses = (0..<128).map(leafIndex(for:))

    for dataBits in 2...5 {
        let dataCount = 1 << dataBits
        let domainSize = 128 * dataCount
        var gpuInputs: [Int32] = []
        var cpuSpectra: [Int32] = []
        for trial in 0..<40 {
            let leaves = leafTables(trial: trial, dataCount: dataCount, rng: &rng)
            try check(leaves.allSatisfy { directWalsh($0).filter { $0 != 0 }.count > 1 },
                      "Unexpected character leaf")
            noncharacterLeafChecks += leaves.count
            let truth = (0..<domainSize).map { point in leaves[addresses[point & 127]][point >> 7] }
            let cpu = cpuFWHT(truth)
            let (predicted, rootS, deltas, special, variedMean) = try liftSpectrum(leaves, nodes, dataCount)
            try check(cpu == predicted, "Haar coefficient formula disagrees with full CPU transform")
            let count = cpu.filter { $0 != 0 }.count
            let expected = nodes.enumerated().reduce(rootS) { value, entry in
                value + ((1 << entry.element.depth) * deltas[entry.offset])
            }
            try check(count == expected, "Haar sparsity formula disagrees")
            if trial == 0 { try check(variedMean.contains(1) && variedMean.contains(2), "Depth-one fixture failed") }
            if trial == 1 { try check(variedMean.contains(0), "Nonconstant common mean at root not exercised") }
            if trial == 2 { try check(rootS == 0 && special == [0, 1, 2], "Zero-root fixture failed") }
            zeroNonconstantCount += special.count
            zeroNonconstantCommonMeanCount += variedMean.count
            if trial < 4 {
                examples.append([
                    "data_bits": dataBits, "trial": trial,
                    "total_sparsity": count, "root_mean_sparsity": rootS, "delta_sparsities_by_node": deltas,
                    "zero_delta_query_nonconstant_nodes": special,
                    "also_common_mean_nonconstant_nodes": variedMean,
                ])
            }
            cases += 1
            gpuInputs.append(contentsOf: truth)
            cpuSpectra.append(contentsOf: cpu)
        }
        let gpu = try metal.transform(gpuInputs, rowSize: domainSize)
        try check(gpu == cpuSpectra, "Complete CPU and Metal spectra disagree")
        countCoefficients += gpu.count
        results.append(["data_bits": dataBits, "cases": 40, "domain_size": domainSize,
                        "coefficients_checked": gpu.count])
    }
    let record: [String: Any] = [
        "schema": "a28-r042-haar-tree-metal-v1",
        "number_of_leaves": 8,
        "independent_query_bits": 7,
        "seed_hex": String(seed, radix: 16),
        "random_cases_per_data_dimension": 37,
        "deterministic_fixtures_per_data_dimension": 3,
        "cases": cases,
        "coefficients_cpu_gpu_and_haar": countCoefficients,
        "cpu_gpu_mismatches": 0,
        "haar_coefficient_mismatches": 0,
        "noncharacter_leaf_checks": noncharacterLeafChecks,
        "node_instances_with_zero_delta_and_query_nonconstant_subtree": zeroNonconstantCount,
        "of_these_with_nonconstant_common_data_mean": zeroNonconstantCommonMeanCount,
        "normalization": "local means and deltas scaled by 8, Walsh coefficients unnormalized",
        "cpu_backend": "exact Int32 CPU FWHT and integer direct intrinsic Walsh",
        "gpu_backend": "exact Int32 staged Metal FWHT",
        "device": metal.device.name,
        "metal4_supported": metal.device.supportsFamily(.metal4),
        "gpu_command_seconds": metal.gpuSeconds,
        "completed_command_buffers": metal.dispatches,
        "dimension_series": results,
        "example_rows": examples,
        "all_checks_passed": true,
    ]
    let output = try JSONSerialization.data(withJSONObject: record, options: [.prettyPrinted, .sortedKeys])
    try (output + Data([10])).write(to: URL(fileURLWithPath: args[1]))
    print("R42 verified \(cases) Haar trees and \(countCoefficients) coefficients.")
}

do {
    try run()
} catch {
    FileHandle.standardError.write(Data(("R42 verifier error: \(error)\n").utf8))
    exit(1)
}
