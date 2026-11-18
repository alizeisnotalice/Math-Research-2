import Foundation
import Metal
import Darwin

struct Settings {var seed,mode,count,unused:UInt32}
struct Candidate {var point,stable,raw,coeff:SIMD4<Float>}
struct Hit:Codable {var point,stable,raw,coeff:[Float];var index,cycle:Int}
struct Row:Codable {
    var mode:Int
    var count=0, zeros=0, nonfinite=0
    var small:Hit?=nil, large:Hit?=nil, error:Hit?=nil
}
func arr(_ v:SIMD4<Float>)->[Float] {[v.x,v.y,v.z,v.w]}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==5,let seconds=Double(a[1]),let batch=Int(a[2]),let seed=UInt32(a[3]),batch>0,
        let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue()
    else {throw NSError(domain:"seconds batch seed out",code:1)}
    let source=try String(contentsOfFile:"wells.metal",encoding:.utf8)
    let opts=MTLCompileOptions();opts.mathMode = .safe
    let lib=try device.makeLibrary(source:source,options:opts)
    let pipe=try device.makeComputePipelineState(function:lib.makeFunction(name:"two_wells")!)
    let buffer=device.makeBuffer(length:batch*MemoryLayout<Candidate>.stride,options:.storageModeShared)!
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let utc=ISO8601DateFormatter();utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    var rows=(0..<30).map {Row(mode:$0)}
    var cycle=0,gpu=0.0
    let start=ProcessInfo.processInfo.systemUptime
    let startUTC=utc.string(from:Date())
    print("start \(startUTC) pid=\(getpid()) device=\(device.name) batch=\(batch) seed=\(seed)")
    fflush(stdout)
    repeat {
        var row=rows[cycle%rows.count]
        var cfg=Settings(seed:seed &+ UInt32(cycle) &* 0x9e3779b9,mode:UInt32(row.mode),count:UInt32(batch),unused:0)
        let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
        e.setComputePipelineState(pipe);e.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:0)
        e.setBuffer(buffer,offset:0,index:1)
        e.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        e.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed,cb.gpuEndTime>cb.gpuStartTime else {
            throw cb.error ?? NSError(domain:"no Metal timing",code:2)
        }
        gpu+=cb.gpuEndTime-cb.gpuStartTime
        let ptr=buffer.contents().bindMemory(to:Candidate.self,capacity:batch)
        var small=row.small?.point[0] ?? .infinity,large=row.large?.point[0] ?? -1,err=row.error?.raw[3] ?? -1
        for i in 0..<batch {
            let x=ptr[i]
            if !x.stable.z.isFinite || !x.raw.w.isFinite {row.nonfinite+=1;continue}
            row.count+=1
            if x.raw.z==0 {row.zeros+=1}
            if x.point.x<small || x.point.x>large || x.raw.w>err {
                let h=Hit(point:arr(x.point),stable:arr(x.stable),raw:arr(x.raw),coeff:arr(x.coeff),index:i,cycle:cycle)
                if x.point.x<small {small=x.point.x;row.small=h}
                if x.point.x>large {large=x.point.x;row.large=h}
                if x.raw.w>err {err=x.raw.w;row.error=h}
            }
        }
        rows[row.mode]=row
        cycle+=1
        if cycle<3 || cycle%10==0 || ProcessInfo.processInfo.systemUptime-start>=seconds {
            print("sample \(utc.string(from:Date())) commands=\(cycle) gpu=\(gpu)")
            fflush(stdout)
        }
    } while ProcessInfo.processInfo.systemUptime-start<seconds
    let summary:[String:Any]=["device":device.name,"submitted":cycle*batch,"cycles":cycle,"seed":Int(seed),"batch":batch,
       "gpuSeconds":gpu,"wallSeconds":ProcessInfo.processInfo.systemUptime-start,"startUTC":startUTC,
       "endUTC":utc.string(from:Date()),"rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:summary,options:[.sortedKeys,.prettyPrinted])
        .write(to:URL(fileURLWithPath:a[4]),options:.atomic)
}
try run()
