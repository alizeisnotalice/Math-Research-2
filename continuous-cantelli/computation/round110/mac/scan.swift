import Foundation
import Metal
import Darwin

struct Settings { var seed, panels, mode, count: UInt32 }
struct Candidate { var roots, diagnostics, margins: SIMD4<Float> }
struct Hit: Codable {
    var roots, diagnostics, margins: [Float]
    var cycle, mode, index: Int
}

func four(_ a: SIMD4<Float>) -> [Float] { [a.x, a.y, a.z, a.w] }

func run() throws {
    let args = CommandLine.arguments
    guard args.count == 6, let panels = UInt32(args[1]), let seconds = Double(args[2]),
          let batch = Int(args[3]), let seed = UInt32(args[4]), batch > 0, panels > 0,
          let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue()
    else { throw NSError(domain: "panels gpu_seconds batch seed output", code: 1) }
    let source = try String(contentsOfFile: "multiwell.metal", encoding: .utf8)
    let settings = MTLCompileOptions()
    settings.mathMode = .safe
    let library = try device.makeLibrary(source: source, options: settings)
    let pipe = try device.makeComputePipelineState(function: library.makeFunction(name: "multiwell")!)
    let buffer = device.makeBuffer(length: batch * MemoryLayout<Candidate>.stride, options: .storageModeShared)!
    let threads = min(pipe.maxTotalThreadsPerThreadgroup, 4 * pipe.threadExecutionWidth)
    let utc = ISO8601DateFormatter()
    utc.formatOptions = [.withInternetDateTime, .withFractionalSeconds]
    var hits = [String: Hit]()
    var counts = [String:Int]()
    var total = 0, valid = 0, cycle = 0
    var gpu = 0.0
    let start = ProcessInfo.processInfo.systemUptime
    let startUTC = utc.string(from: Date())
    print("start \(startUTC) pid=\(getpid()) device=\(device.name) panels=\(panels) batch=\(batch) seed=\(seed)")
    fflush(stdout)
    repeat {
        let mode = cycle % 3
        var cfg = Settings(seed: seed &+ UInt32(cycle) &* 0x9e3779b9, panels: panels, mode: UInt32(mode), count: UInt32(batch))
        let command = queue.makeCommandBuffer()!, encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipe)
        encoder.setBytes(&cfg, length: MemoryLayout<Settings>.stride, index: 0)
        encoder.setBuffer(buffer, offset: 0, index: 1)
        encoder.dispatchThreads(MTLSize(width: batch, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: threads, height: 1, depth: 1))
        encoder.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed, command.gpuEndTime > command.gpuStartTime
        else { throw command.error ?? NSError(domain: "Metal data or timing unavailable", code: 2) }
        gpu += command.gpuEndTime - command.gpuStartTime
        total += batch
        let output = buffer.contents().bindMemory(to: Candidate.self, capacity: batch)
        for i in 0..<batch {
            let a = output[i]
            if a.diagnostics.x < 0 { continue }
            guard four(a.margins).allSatisfy(\.isFinite) else { counts["nonfinite", default:0] += 1; continue }
            valid += 1
            if a.margins.y < -1e-6 { counts["jensen_negative", default:0] += 1 }
            if a.diagnostics.z > 0 { counts["rate_below", default:0] += 1 }
            else { counts["rate_above", default:0] += 1 }
            let hit = Hit(roots: four(a.roots), diagnostics: four(a.diagnostics), margins: four(a.margins),
                          cycle: cycle, mode: mode, index: i)
            if a.margins.x < (hits["min_scaled_grid_gap"]?.margins[0] ?? .infinity) { hits["min_scaled_grid_gap"] = hit }
            if a.margins.y < (hits["min_shifted_jensen"]?.margins[1] ?? 0) { hits["min_shifted_jensen"] = hit }
            let q = a.diagnostics.z
            let key = q < 0 ? "negative_closest_rate" : "positive_closest_rate"
            if abs(q) < abs(hits[key]?.diagnostics[2] ?? .infinity) { hits[key] = hit }
            if a.margins.w > (hits["max_grid_sign_changes"]?.margins[3] ?? -1) { hits["max_grid_sign_changes"] = hit }
            if abs(q) < 0.005 && a.margins.x < (hits["near_rate_min_grid_gap"]?.margins[0] ?? .infinity) {
                hits["near_rate_min_grid_gap"] = hit
            }
        }
        cycle += 1
        if cycle % 4 == 0 || cycle < 3 || gpu >= seconds {
            print("sample \(utc.string(from: Date())) batches=\(cycle) submitted=\(total) valid=\(valid) gpu=\(gpu)")
            fflush(stdout)
        }
    } while gpu < seconds
    let obj: [String:Any] = ["device": device.name, "submitted":total, "valid":valid, "counts":counts,
      "cycles":cycle, "panels":Int(panels), "batch":batch, "seed":Int(seed), "gpuSeconds":gpu,
      "wallSeconds":ProcessInfo.processInfo.systemUptime-start, "startUTC":startUTC,
      "endUTC":utc.string(from: Date()), "hits":try JSONSerialization.jsonObject(with:JSONEncoder().encode(hits))]
    try JSONSerialization.data(withJSONObject:obj,options:[.prettyPrinted,.sortedKeys])
       .write(to:URL(fileURLWithPath:args[5]),options:.atomic)
}
try run()
