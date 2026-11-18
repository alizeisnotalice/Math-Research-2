import Foundation
import Metal
import Darwin
struct TablePoint {var c1,c2,s1,s2:SIMD4<Float>}
struct Settings {var panels,family,seed:UInt32;var amplitude:Float}
struct Candidate {var metrics,c1,c2,s1,s2:SIMD4<Float>}
struct Hit: Codable {var metrics,cosine,sine:[Float];var cycle,index:Int}
struct Row: Codable {var family:String;var amplitude:Float;var examined:Int;var below1e3,below1e4,below1e5:Int;var best:Hit?}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==6,let panels=Int(a[1]),let seconds=Double(a[2]),let batch=Int(a[3]),let seed=UInt32(a[4]),panels>=32,batch>0 else {
        throw NSError(domain:"panels seconds batch seed out",code:1)
    }
    guard let d=MTLCreateSystemDefaultDevice(),let queue=d.makeCommandQueue() else{throw NSError(domain:"Metal device",code:2)}
    let source=try String(contentsOfFile:"harmonics.metal",encoding:.utf8)
    let o=MTLCompileOptions();o.mathMode = .safe
    let lib=try d.makeLibrary(source:source,options:o)
    let pipe=try d.makeComputePipelineState(function:lib.makeFunction(name:"scan_harmonics")!)
    var table=[TablePoint]()
    for j in 0..<panels {
        let t=2*Double.pi*Double(j)/Double(panels)
        let c=(1...8).map{Float(cos(Double($0)*t))}
        let s=(1...8).map{Float(sin(Double($0)*t))}
        table.append(TablePoint(c1:SIMD4(c[0],c[1],c[2],c[3]),c2:SIMD4(c[4],c[5],c[6],c[7]),
                                s1:SIMD4(s[0],s[1],s[2],s[3]),s2:SIMD4(s[4],s[5],s[6],s[7])))
    }
    let buf=d.makeBuffer(bytes:table,length:table.count*MemoryLayout<TablePoint>.stride,options:.storageModeShared)!
    let out=d.makeBuffer(length:batch*MemoryLayout<Candidate>.stride,options:.storageModeShared)!
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let families=["low-only","high-only","corrected","subperiod","mixed"]
    let amplitudes:[Float]=[0.015625,0.03125,0.0625,0.125,0.25,0.5,0.75,1]
    var rows=[Row]()
    for f in families {for x in amplitudes {rows.append(Row(family:f,amplitude:x,examined:0,below1e3:0,below1e4:0,below1e5:0))}}
    var cycle=0
    var gpu=0.0
    var finite=0
    let started=ProcessInfo.processInfo.systemUptime
    print("device=\(d.name) panels=\(panels) batch=\(batch) seed=\(seed)");fflush(stdout)
    repeat {
        for j in 0..<rows.count {
            var settings=Settings(panels:UInt32(panels),family:UInt32(j/amplitudes.count),seed:seed &+ UInt32(cycle*rows.count+j),amplitude:rows[j].amplitude)
            let cb=queue.makeCommandBuffer()!, e=cb.makeComputeCommandEncoder()!
            e.setComputePipelineState(pipe);e.setBuffer(buf,offset:0,index:0)
            e.setBytes(&settings,length:MemoryLayout<Settings>.stride,index:1);e.setBuffer(out,offset:0,index:2)
            e.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
            e.endEncoding();cb.commit();cb.waitUntilCompleted()
            guard cb.status == .completed else{throw cb.error ?? NSError(domain:"GPU failed",code:3)}
            let elapsed=cb.gpuEndTime-cb.gpuStartTime
            guard elapsed>0 else{throw NSError(domain:"GPU timestamp unavailable",code:4)}
            gpu+=elapsed
            let v=out.contents().bindMemory(to:Candidate.self,capacity:batch)
            var row=rows[j]
            for i in 0..<batch {
                let value=v[i].metrics.w
                if !value.isFinite {continue}
                finite+=1
                if value<1e-3 {row.below1e3+=1}
                if value<1e-4 {row.below1e4+=1}
                if value<1e-5 {row.below1e5+=1}
                if value<(row.best?.metrics[3] ?? .infinity) {
                    let h=v[i]
                    row.best=Hit(metrics:[h.metrics.x,h.metrics.y,h.metrics.z,h.metrics.w],
                        cosine:[h.c1.x,h.c1.y,h.c1.z,h.c1.w,h.c2.x,h.c2.y,h.c2.z,h.c2.w],
                        sine:[h.s1.x,h.s1.y,h.s1.z,h.s1.w,h.s2.x,h.s2.y,h.s2.z,h.s2.w],cycle:cycle,index:i)
                }
            }
            row.examined+=batch;rows[j]=row
        }
        cycle+=1
        print("cycles=\(cycle) evaluated=\(cycle*rows.count*batch) gpu_seconds=\(gpu) wall=\(ProcessInfo.processInfo.systemUptime-started)");fflush(stdout)
    } while gpu<seconds
    let e=JSONEncoder()
    let payload:[String:Any]=["device":d.name,"panels":panels,"batch":batch,"seed":seed,"cycles":cycle,
        "gpuSeconds":gpu,"evaluated":cycle*rows.count*batch,"finite":finite,
        "wallSeconds":ProcessInfo.processInfo.systemUptime-started,"rows":try JSONSerialization.jsonObject(with:e.encode(rows))]
    try JSONSerialization.data(withJSONObject:payload,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[5]),options:.atomic)
}
try run()
