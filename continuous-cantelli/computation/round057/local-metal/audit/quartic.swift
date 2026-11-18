import Foundation
import Metal
import Darwin

struct Point: Codable {
    let r, m, edge, b, eps, h, k, A, B: Float
    let tuned: UInt32
}

struct Row: Codable {
    let point: Point
    let exact, leading, gap, offset: Float
    let index, salt: UInt32
}

func gaussLegendre(_ count: Int) -> [SIMD2<Float>] {
    var nodes = [SIMD2<Float>](repeating: SIMD2.zero, count: count)
    for i in 0..<(count/2) {
        var z = cos(Double.pi * (Double(i) + 0.75) / (Double(count) + 0.5))
        var derivative = 0.0
        for _ in 0..<12 {
            var previous = 0.0
            var value = 1.0
            for j in 1...count {
                let secondPrevious = previous
                previous = value
                value = ((2 * Double(j) - 1) * z * previous
                         - (Double(j) - 1) * secondPrevious) / Double(j)
            }
            derivative = Double(count) * (z * value - previous) / (z * z - 1)
            let adjustment = value / derivative
            z -= adjustment
            if abs(adjustment) < 2e-16 { break }
        }
        let weight = Float(1 / ((1-z*z) * derivative*derivative))
        nodes[i] = SIMD2(Float((1-z)/2), weight)
        nodes[count-1-i] = SIMD2(Float((1+z)/2), weight)
    }
    return nodes
}

let sampleCount = 8192
let indices = Array(UInt32(4_500_000)..<UInt32(4_500_000 + sampleCount)) + [984450415]
let order = 96
precondition(MemoryLayout<Row>.stride == 64)
guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
    fatalError("Missing Metal device")
}
let source = try String(contentsOfFile: "quartic.metal", encoding: .utf8)
let options = MTLCompileOptions()
options.mathMode = .safe
let library = try device.makeLibrary(source: source, options: options)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "quartic")!)
let nodes = gaussLegendre(order)
let nodeBuffer = nodes.withUnsafeBytes {
    device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
}
let indexBuffer = indices.withUnsafeBytes {
    device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
}
let output = device.makeBuffer(length: indices.count * MemoryLayout<Row>.stride, options: .storageModeShared)!
let command = queue.makeCommandBuffer()!
let encoder = command.makeComputeCommandEncoder()!
var params = SIMD4(UInt32(indices.count), UInt32(order), 0, 0)
encoder.setComputePipelineState(pipeline)
encoder.setBuffer(nodeBuffer, offset: 0, index: 0)
encoder.setBuffer(output, offset: 0, index: 1)
encoder.setBytes(&params, length: MemoryLayout<SIMD4<UInt32>>.stride, index: 2)
encoder.setBuffer(indexBuffer, offset: 0, index: 3)
encoder.dispatchThreads(MTLSize(width: indices.count, height: 1, depth: 1),
                        threadsPerThreadgroup: MTLSize(width: 256, height: 1, depth: 1))
encoder.endEncoding()
command.commit()
command.waitUntilCompleted()
guard command.status == .completed else {
    fatalError(String(describing: command.error))
}
let rows = output.contents().bindMemory(to: Row.self, capacity: indices.count)
let result = Array(UnsafeBufferPointer(start: rows, count: indices.count))
try JSONEncoder().encode(result).write(to: URL(fileURLWithPath: "audit-input.json"), options: .atomic)
print("PID=\(getpid()) rows=\(indices.count) command_seconds=\(command.gpuEndTime-command.gpuStartTime)")
