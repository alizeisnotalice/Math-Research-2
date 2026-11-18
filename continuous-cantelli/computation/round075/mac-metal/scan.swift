import Foundation
import Metal
import Darwin
struct Settings {var panels,base,seed,density:UInt32;var power:Float}
struct Candidate {var scores:SIMD4<Float>;var config:SIMD4<UInt32>}
struct Hit: Codable {var scores:[Float];var mask,depth:UInt32;var height:Float;var cycle,index:Int}
struct Row: Codable {var base:Int;var exponent:Float;var densityIndex:Int;var examined:Int;var best6,best3:Hit?}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==6,let panels=Int(a[1]),let seconds=Double(a[2]),let batch=Int(a[3]),let seed=UInt32(a[4]),
          panels>0,batch>0 else{throw NSError(domain:"panels seconds batch seed out",code:1)}
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else{throw NSError(domain:"Metal device",code:2)}
    let source=try String(contentsOfFile:"cantor.metal",encoding:.utf8)
    let options=MTLCompileOptions();options.mathMode = .safe
    let lib=try device.makeLibrary(source:source,options:options)
    let pipe=try device.makeComputePipelineState(function:lib.makeFunction(name:"cantor_scan")!)
    let buffer=device.makeBuffer(length:batch*MemoryLayout<Candidate>.stride,options:.storageModeShared)!
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    var rows=[Row]()
    let powers:[Float]=[0.25,0.5,1,2,3]
    for b in [3,4,5,7,9,13,19] {for p in powers {for d in 0..<3 {rows.append(Row(base:b,exponent:p,densityIndex:d,examined:0))}}}
    var cycle=0,finite=0
    var gpu=0.0
    let start=ProcessInfo.processInfo.systemUptime
    print("device=\(device.name) strata=\(rows.count) panels=\(panels) batch=\(batch)");fflush(stdout)
    repeat {
        for j in 0..<rows.count {
            var row=rows[j]
            var cfg=Settings(panels:UInt32(panels),base:UInt32(row.base),seed:seed &+ UInt32(cycle*rows.count+j),
               density:UInt32(row.densityIndex),power:row.exponent)
            let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
            e.setComputePipelineState(pipe);e.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:0)
            e.setBuffer(buffer,offset:0,index:1)
            e.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
            e.endEncoding();cb.commit();cb.waitUntilCompleted()
            guard cb.status == .completed else{throw cb.error ?? NSError(domain:"GPU failed",code:3)}
            let elapsed=cb.gpuEndTime-cb.gpuStartTime
            guard elapsed>0 else{throw NSError(domain:"GPU timing unavailable",code:4)}
            gpu+=elapsed
            let output=buffer.contents().bindMemory(to:Candidate.self,capacity:batch)
            for i in 0..<batch {
                let x=output[i]
                if !x.scores.x.isFinite || !x.scores.y.isFinite {continue}
                finite+=1
                let win6=x.scores.x<(row.best6?.scores[0] ?? .infinity)
                let win3=x.scores.y<(row.best3?.scores[1] ?? .infinity)
                if win6 || win3 {
                    let hit=Hit(scores:[x.scores.x,x.scores.y,x.scores.z,x.scores.w],mask:x.config.x,depth:x.config.y,
                        height:Float(bitPattern:x.config.z),cycle:cycle,index:i)
                    if win6 {row.best6=hit};if win3 {row.best3=hit}
                }
            }
            row.examined+=batch;rows[j]=row
        }
        cycle+=1
        print("cycle=\(cycle) examined=\(cycle*rows.count*batch) gpu_seconds=\(gpu) wall=\(ProcessInfo.processInfo.systemUptime-start)");fflush(stdout)
    } while gpu<seconds
    let o:[String:Any]=["device":device.name,"panels":panels,"batch":batch,"seed":seed,"cycles":cycle,"gpuSeconds":gpu,
      "finite":finite,"evaluated":cycle*rows.count*batch,"wallSeconds":ProcessInfo.processInfo.systemUptime-start,
      "rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:o,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[5]),options:.atomic)
}
try run()
