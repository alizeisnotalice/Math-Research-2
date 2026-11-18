import Darwin
import Foundation
import Metal

struct Settings {var mode,seed,steps,unused:UInt32}
struct Candidate {var shape,origin,endpoint,sweep:SIMD4<Float>}
struct Hit:Codable {var shape,origin,endpoint,sweep:[Float];var index,cycle:Int}
struct Row:Codable {
    var mode:Int
    var count=0, endPass=0, interiorOnly=0, onlyOnePhase=0, nonfinite=0
    var winner:Hit?=nil, selectedPhase:Hit?=nil, fixed:Hit?=nil
}
func a4(_ v:SIMD4<Float>)->[Float] {[v.x,v.y,v.z,v.w]}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==6,let seconds=Double(a[1]),let batch=Int(a[2]),let seed=UInt32(a[3]),let steps=UInt32(a[4]),
          batch>0,steps>2,let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue()
    else {throw NSError(domain:"seconds batch seed steps out",code:1)}
    let source=try String(contentsOfFile:"sine.metal",encoding:.utf8)
    let options=MTLCompileOptions();options.mathMode = .safe
    let library=try device.makeLibrary(source:source,options:options)
    let pipe=try device.makeComputePipelineState(function:library.makeFunction(name:"positive_sine")!)
    let buffer=device.makeBuffer(length:batch*MemoryLayout<Candidate>.stride,options:.storageModeShared)!
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let utc=ISO8601DateFormatter();utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    var rows=(0..<4).map {Row(mode:$0)}
    var cycle=0,gpu=0.0
    let start=ProcessInfo.processInfo.systemUptime,startUTC=utc.string(from:Date())
    print("start \(startUTC) pid=\(getpid()) device=\(device.name) batch=\(batch) steps=\(steps)")
    fflush(stdout)
    repeat {
        var row=rows[cycle%4],cfg=Settings(mode:UInt32(cycle%4),seed:seed &+ UInt32(cycle) &* 0x9e3779b9,steps:steps,unused:0)
        let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
        e.setComputePipelineState(pipe);e.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:0)
        e.setBuffer(buffer,offset:0,index:1)
        e.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        e.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed,cb.gpuEndTime>cb.gpuStartTime else {
            throw cb.error ?? NSError(domain:"Metal timing or execution",code:2)
        }
        gpu += cb.gpuEndTime-cb.gpuStartTime
        let p=buffer.contents().bindMemory(to:Candidate.self,capacity:batch)
        var best=row.winner.map {max($0.sweep[0],$0.sweep[1])} ?? 0
        var bestPhase=row.selectedPhase.map {max($0.sweep[0],$0.sweep[1])} ?? 0
        for i in 0..<batch {
            let v=p[i],top=max(v.sweep.x,v.sweep.y)
            if !top.isFinite || !v.endpoint.w.isFinite {row.nonfinite+=1;continue}
            row.count+=1
            if i==0 && row.mode==3 && row.fixed==nil {
                row.fixed=Hit(shape:a4(v.shape),origin:a4(v.origin),endpoint:a4(v.endpoint),sweep:a4(v.sweep),index:i,cycle:cycle)
            }
            let limit=max(max(v.endpoint.x,v.endpoint.y),max(v.endpoint.z,v.endpoint.w))
            if limit>=0.999 {continue}
            row.endPass+=1
            if top<=1.001 {continue}
            row.interiorOnly+=1
            let phaseOnly=min(v.sweep.x,v.sweep.y)<0.999
            if phaseOnly {row.onlyOnePhase+=1}
            if top>best || (phaseOnly && top>bestPhase){
                let hit=Hit(shape:a4(v.shape),origin:a4(v.origin),endpoint:a4(v.endpoint),sweep:a4(v.sweep),index:i,cycle:cycle)
                if top>best {best=top;row.winner=hit}
                if phaseOnly && top>bestPhase {bestPhase=top;row.selectedPhase=hit}
            }
        }
        rows[row.mode]=row;cycle+=1
        if cycle%10==0 || cycle<=2 || ProcessInfo.processInfo.systemUptime-start>=seconds {
            print("sample \(utc.string(from:Date())) commands=\(cycle) gpu=\(gpu)")
            fflush(stdout)
        }
    } while ProcessInfo.processInfo.systemUptime-start<seconds
    let summary:[String:Any]=["device":device.name,"submitted":cycle*batch,"cycles":cycle,"seed":Int(seed),
        "steps":Int(steps),"batch":batch,"gpuSeconds":gpu,"wallSeconds":ProcessInfo.processInfo.systemUptime-start,
        "startUTC":startUTC,"endUTC":utc.string(from:Date()),"rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:summary,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[5]),options:.atomic)
}
try run()
