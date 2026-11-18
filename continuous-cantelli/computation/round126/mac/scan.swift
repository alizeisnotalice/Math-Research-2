import Foundation
import Metal
import Darwin

struct Settings { var seed, panels, mode, count: UInt32 }
struct Candidate { var geometry, shape, overlap, values: SIMD4<Float> }
struct Hit: Codable {
    var geometry, shape, overlap, values: [Float]
    var cycle, mode, index: Int
}
func four(_ a: SIMD4<Float>) -> [Float] { [a.x,a.y,a.z,a.w] }
func run() throws {
    let args=CommandLine.arguments
    guard args.count==6,let panels=UInt32(args[1]),let seconds=Double(args[2]),let batch=Int(args[3]),
          let seed=UInt32(args[4]),batch>0,panels>0,let device=MTLCreateSystemDefaultDevice(),
          let queue=device.makeCommandQueue() else {throw NSError(domain:"panels gpu_seconds batch seed output",code:1)}
    let source=try String(contentsOfFile:"overlap.metal",encoding:.utf8)
    let options=MTLCompileOptions()
    options.mathMode = .safe
    let library=try device.makeLibrary(source:source,options:options)
    let pipe=try device.makeComputePipelineState(function:library.makeFunction(name:"overlaps")!)
    let buffer=device.makeBuffer(length:batch*MemoryLayout<Candidate>.stride,options:.storageModeShared)!
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let clock=ISO8601DateFormatter()
    clock.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    var hits=[String:Hit](), counts=[String:Int]()
    var submitted=0, accepted=0, cycle=0
    var gpu=0.0
    let start=ProcessInfo.processInfo.systemUptime
    let startUTC=clock.string(from:Date())
    print("start \(startUTC) pid=\(getpid()) device=\(device.name) panels=\(panels) batch=\(batch) seed=\(seed)")
    fflush(stdout)
    repeat {
        let mode=cycle%3
        var cfg=Settings(seed:seed &+ UInt32(cycle) &* 0x9e3779b9,panels:panels,mode:UInt32(mode),count:UInt32(batch))
        let command=queue.makeCommandBuffer()!,encoder=command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipe)
        encoder.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:0)
        encoder.setBuffer(buffer,offset:0,index:1)
        encoder.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        encoder.endEncoding();command.commit();command.waitUntilCompleted()
        guard command.status == .completed, command.gpuEndTime>command.gpuStartTime
        else {throw command.error ?? NSError(domain:"Metal timing or execution",code:2)}
        gpu+=command.gpuEndTime-command.gpuStartTime
        submitted+=batch
        let output=buffer.contents().bindMemory(to:Candidate.self,capacity:batch)
        for i in 0..<batch {
            let a=output[i]
            if a.values.x<0 {continue}
            accepted+=1
            if a.overlap.x==0 {counts["no_safe_shift",default:0]+=1;continue}
            if !a.values.z.isFinite {counts["nonfinite",default:0]+=1;continue}
            if a.values.z>0.001 {counts["improvement_positive",default:0]+=1}
            if a.values.w>0 {counts["separates_at_energy_1",default:0]+=1}
            let stratum=Int(log2(a.geometry.x))/2
            let key="best_gain_stratum_\(stratum)"
            var keys=[String]()
            if a.values.z>(hits[key]?.values[2] ?? -.infinity) {keys.append(key)}
            if a.values.x<1 && a.values.y<1 && a.overlap.w>(hits["best_strict_complement"]?.overlap[3] ?? 0) {
                keys.append("best_strict_complement")
            }
            if max(a.overlap.y,a.overlap.z)<0.90 && a.values.x<1 && a.values.y<1 && a.overlap.w>(hits["best_strict_q_below_point9"]?.overlap[3] ?? 0) {
                keys.append("best_strict_q_below_point9")
            }
            if !keys.isEmpty {
                let hit=Hit(geometry:four(a.geometry),shape:four(a.shape),overlap:four(a.overlap),
                    values:four(a.values),cycle:cycle,mode:mode,index:i)
                for k in keys {hits[k]=hit}
            }
        }
        cycle+=1
        if cycle%4==0 || cycle<3 || gpu>=seconds {
            print("sample \(clock.string(from:Date())) batches=\(cycle) submitted=\(submitted) accepted=\(accepted) gpu=\(gpu)")
            fflush(stdout)
        }
    } while gpu<seconds
    let obj:[String:Any]=["device":device.name,"submitted":submitted,"accepted":accepted,"counts":counts,
        "cycles":cycle,"seed":Int(seed),"panels":Int(panels),"batch":batch,"gpuSeconds":gpu,
        "wallSeconds":ProcessInfo.processInfo.systemUptime-start,"startUTC":startUTC,
        "endUTC":clock.string(from:Date()),"hits":try JSONSerialization.jsonObject(with:JSONEncoder().encode(hits))]
    try JSONSerialization.data(withJSONObject:obj,options:[.prettyPrinted,.sortedKeys])
        .write(to:URL(fileURLWithPath:args[5]),options:.atomic)
}
try run()
