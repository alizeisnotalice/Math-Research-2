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


let rows = try String(contentsOfFile: CommandLine.arguments[1], encoding: .utf8).split(separator: "\n")
let gpu = try MetalWalsh()
var count = 0
for line in rows {
    let col = line.split(separator: " ")
    if Int(col[3]) == 0 { continue }
    let w = UInt64(col[5], radix: 16)!
    var signs = [Int32]()
    for t in 0..<64 {
        for x in 0..<64 {
            let a = x & 1, b = (x >> 1) & 1, c = (x >> 2) & 1
            let h = (x >> 3) & 1, u = (x >> 4) & 1, v = (x >> 5) & 1
            let p = parity(t & (1 | (a << 1) | (b << 2)))
            let z = parity(t & ((a << 3) | (b << 4)))
            let g = (1 ^ v) & (h == 0 ? p : z)
            let extra = t >= 32 ? Int((w >> x) & 1) : 0
            signs.append(sign((a & b) ^ (c & h) ^ (u & v) ^ g ^ extra))
        }
    }
    let values = try gpu.transform(signs,rowSize:64)
    try check(values.allSatisfy{abs($0)==8},"Walsh magnitude")
    count += values.count
}
print("device=\(gpu.device.name) coefficients=\(count) commandBuffers=\(gpu.dispatches) gpuSeconds=\(gpu.gpuSeconds)")
