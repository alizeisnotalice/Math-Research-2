import Foundation
import Metal
import Darwin

struct Point: Codable {
    let r, m, edge, b, eps, h, k, A, B: Float
    let tuned: UInt32
}
struct Profile: Codable {
    let origin: String
    let index: UInt32
    let point: Point
}
struct Result: Codable {
    let exact, leading, A, B: Float
}
struct Report: Codable {
    let order: Int
    let commandSeconds: Double
    let values: [Result]
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

let profiles = try JSONDecoder().decode([Profile].self, from: Data(contentsOf: URL(fileURLWithPath: "points.json")))
let points = profiles.map(\.point)
precondition(!points.isEmpty && MemoryLayout<Point>.stride == 40 && MemoryLayout<Result>.stride == 16)
guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
    fatalError("Metal unavailable")
}
let base = try String(contentsOfFile: "../quartic.metal", encoding: .utf8)
let new = try String(contentsOfFile: "kernel.metal", encoding: .utf8)
let shared = base.components(separatedBy: "kernel void quartic")[0]
let options = MTLCompileOptions()
options.mathMode = .safe
let library = try device.makeLibrary(source: shared + "\n" + new, options: options)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "matched")!)
let input = points.withUnsafeBytes {
    device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
}
let output = device.makeBuffer(length: points.count*MemoryLayout<Result>.stride, options: .storageModeShared)!
var reports: [Report] = []
for order in [96,128,160] {
    let nodes = gaussLegendre(order)
    let buffer = nodes.withUnsafeBytes {
        device.makeBuffer(bytes: $0.baseAddress!, length: $0.count, options: .storageModeShared)!
    }
    let command = queue.makeCommandBuffer()!
    let encoder = command.makeComputeCommandEncoder()!
    var nodeCount = UInt32(order)
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(buffer, offset: 0, index: 0)
    encoder.setBuffer(output, offset: 0, index: 1)
    encoder.setBytes(&nodeCount, length: MemoryLayout<UInt32>.stride, index: 2)
    encoder.setBuffer(input, offset: 0, index: 3)
    encoder.dispatchThreads(MTLSize(width: points.count, height: 1, depth: 1),
                            threadsPerThreadgroup: MTLSize(width: 256, height: 1, depth: 1))
    encoder.endEncoding()
    command.commit()
    command.waitUntilCompleted()
    guard command.status == .completed else { fatalError(String(describing: command.error)) }
    let data = output.contents().bindMemory(to: Result.self, capacity: points.count)
    reports.append(Report(order: order, commandSeconds: command.gpuEndTime-command.gpuStartTime,
                          values: Array(UnsafeBufferPointer(start: data, count: points.count))))
}
let json = JSONEncoder()
json.outputFormatting = [.sortedKeys, .prettyPrinted]
try json.encode(reports).write(to: URL(fileURLWithPath: "matched-gpu.json"), options: .atomic)
print("PID=\(getpid()) cases=\(points.count) orders=96,128,160")
