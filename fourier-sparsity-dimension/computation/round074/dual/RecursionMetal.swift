import Foundation
import Metal

struct Config { var first: UInt32; var length: UInt32 }
struct AuditFailure: Error { let details: String }

func check(_ condition: Bool, _ message: String) throws {
    guard condition else { throw AuditFailure(details: message) }
}

func hashCode(_ input: UInt64) -> UInt64 {
    var x = input
    x = (x ^ (x >> 30)) &* 0xbf58476d1ce4e5b9
    x = (x ^ (x >> 27)) &* 0x94d049bb133111eb
    return x ^ (x >> 31)
}

final class Slot {
    let memory: MTLBuffer
    var command: MTLCommandBuffer?
    var first = 0
    var size = 0
    init(device: MTLDevice, limit: Int, width: Int) {
        memory = device.makeBuffer(length: limit * width * 4, options: .storageModeShared)!
    }
}

func main() throws {
    let args = CommandLine.arguments
    try check(args.count == 5, "usage: recursion_metal result.json first count batch")
    guard let first = Int(args[2]), let count = Int(args[3]), let batch = Int(args[4]),
          first >= 0, count > 0, batch > 0, batch <= 65536 else {
        throw AuditFailure(details: "Invalid range")
    }
    let names = ["RecursionCommon.h", "RecursionKernel.metal"]
    let source = try names.map { try String(contentsOfFile: $0, encoding: .utf8) }.joined(separator: "\n")
    guard let device = MTLCreateSystemDefaultDevice() else { throw AuditFailure(details: "Metal unavailable") }
    let library = try device.makeLibrary(source: source, options: nil)
    let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "audit")!)
    let width = 40
    let buffers = (0..<4).map { _ in Slot(device: device, limit: batch, width: width) }
    let queues = (0..<2).map { _ in device.makeCommandQueue()! }
    let start = ProcessInfo.processInfo.systemUptime
    let epoch = Date().timeIntervalSince1970 - start
    var totals = Array(repeating: UInt64(0), count: width - 1)
    var timings = [[Double]]()
    var deviceSeconds: Double = 0
    var checksum: UInt64 = 0
    var controls: [String: [UInt64]] = [:]

    func collect(_ slot: Slot) throws {
        guard let command = slot.command else { return }
        command.waitUntilCompleted()
        if let err = command.error { throw err }
        try check(command.status == .completed, "Metal command did not complete")
        let start = command.gpuStartTime
        let end = command.gpuEndTime
        deviceSeconds += end - start
        timings.append([start + epoch, end + epoch])
        let raw = slot.memory.contents().bindMemory(to: UInt32.self, capacity: slot.size * width)
        for local in 0..<slot.size {
            let base = local * width
            let id = slot.first + local
            try check(raw[base+width-1] == 1, "Metal audit failed for id \(id)")
            var sig = hashCode(UInt64(id))
            for j in 0..<(width-1) {
                let item = UInt64(raw[base+j])
                totals[j] += item
                sig = hashCode(sig ^ item)
            }
            checksum ^= sig
            if id < 6 { controls[String(id)] = (0..<(width-1)).map { UInt64(raw[base+$0]) } }
        }
        slot.command = nil
    }

    let taskStart = ProcessInfo.processInfo.systemUptime
    var batchIndex = 0
    for position in stride(from: first, to: first+count, by: batch) {
        let slot = buffers[batchIndex % buffers.count]
        try collect(slot)
        let queue = queues[batchIndex % queues.count]
        let command = queue.makeCommandBuffer()!
        let encoder = command.makeComputeCommandEncoder()!
        let length = min(batch, first+count-position)
        var config = Config(first: UInt32(position), length: UInt32(length))
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(slot.memory, offset: 0, index: 0)
        encoder.setBytes(&config, length: MemoryLayout<Config>.size, index: 1)
        encoder.dispatchThreads(MTLSize(width: length, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: 64, height: 1, depth: 1))
        encoder.endEncoding()
        slot.command = command
        slot.first = position
        slot.size = length
        command.commit()
        batchIndex += 1
    }
    for slot in buffers { try collect(slot) }
    let result: [String: Any] = [
        "first_id": first, "count": count, "device": device.name,
        "pid": ProcessInfo.processInfo.processIdentifier,
        "metal_wall_seconds": ProcessInfo.processInfo.systemUptime-taskStart,
        "gpu_command_seconds": deviceSeconds, "gpu_command_buffers": timings.count,
        "gpu_command_intervals_epoch": timings, "positional_checksum64": String(format:"%016llx", checksum),
        "sums": totals, "controls": controls, "all_checks_passed": true
    ]
    let output = try JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys])
    try output.write(to: URL(fileURLWithPath: args[1]))
}

try main()
