import Darwin
import Foundation
import Metal

struct Probe {var logK,logE:Float;var nearOne,pad:UInt32}
struct Pick:Codable {
    var group,index:Int
    var logK,logE,rootY,rootX,scaledOneGap,residual:Float
    var nearOne:UInt32
}
struct Row:Codable {
    var group,count,nonfinite,monotonicReversals:Int
    var minimumY,maximumY,maxAbsResidual:Float
    var selected:[Pick]
}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==4,let count=Int(a[1]),let tile=Int(a[2]),count>=16,tile>0,
      let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue()
    else {throw NSError(domain:"count_per_stratum tile output",code:1)}
    let opts=MTLCompileOptions();opts.mathMode = .safe
    let source=try String(contentsOfFile:"marchaud.metal",encoding:.utf8)
    let lib=try device.makeLibrary(source:source,options:opts)
    let pipe=try device.makeComputePipelineState(function:lib.makeFunction(name:"marchaud")!)
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let utc=ISO8601DateFormatter();utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    let before=ProcessInfo.processInfo.systemUptime,start=utc.string(from:Date())
    var rows=[Row](),gpu=0.0,commands=0
    var sourceSets=[[Probe]]()
    let anchor:[Double]=[0.000000000001,0.00001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.9999,1-pow(2.0,-20),1-pow(2.0,-40)]
    sourceSets.append(anchor.map{Probe(logK:Float(log($0)),logE:Float(log1p(-$0)),nearOne:$0>0.7 ? 1:0,pad:0)})
    sourceSets.append((0..<count).map{i in
        let k=0.001+0.998*Double(i)/Double(count-1)
        return Probe(logK:Float(log(k)),logE:Float(log1p(-k)),nearOne:k>0.7 ? 1:0,pad:0)
    })
    sourceSets.append((0..<count).map{i in
        let ell=0.7+599.3*Double(i)/Double(count-1)
        return Probe(logK:Float(-ell),logE:Float(log1p(-exp(-ell))),nearOne:0,pad:0)
    })
    sourceSets.append((0..<count).map{i in
        let e=0.7+79.3*Double(i)/Double(count-1)
        return Probe(logK:Float(log1p(-exp(-e))),logE:Float(-e),nearOne:1,pad:0)
    })
    print("start \(start) pid=\(getpid()) device=\(device.name) per_stratum=\(count) tile=\(tile)")
    fflush(stdout)
    var total=0
    for (group,probes) in sourceSets.enumerated(){
        var picks=[Pick](),minY=Float.infinity,maxY:Float=0,resid:Float=0,previous:Float?=nil,reversed=0,bad=0
        let selected=Set(group==0 ? Array(probes.indices) : (0...16).map {$0*(count-1)/16})
        for begin in stride(from:0,to:probes.count,by:tile){
            let end=min(probes.count,begin+tile),n=end-begin
            let input=Array(probes[begin..<end])
            let data=device.makeBuffer(bytes:input,length:n*MemoryLayout<Probe>.stride,options:.storageModeShared)!
            let output=device.makeBuffer(length:n*MemoryLayout<SIMD4<Float>>.stride,options:.storageModeShared)!
            let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
            e.setComputePipelineState(pipe);e.setBuffer(data,offset:0,index:0);e.setBuffer(output,offset:0,index:1)
            e.dispatchThreads(MTLSize(width:n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
            e.endEncoding();cb.commit();cb.waitUntilCompleted()
            guard cb.status == .completed,cb.gpuEndTime>cb.gpuStartTime else{throw cb.error ?? NSError(domain:"Metal timing/execution",code:2)}
            gpu+=cb.gpuEndTime-cb.gpuStartTime;commands+=1
            let ptr=output.contents().bindMemory(to:SIMD4<Float>.self,capacity:n)
            for j in 0..<n {
                let p=input[j],r=ptr[j],i=begin+j
                if !r.x.isFinite || !r.w.isFinite{bad+=1;continue}
                minY=min(minY,r.x);maxY=max(maxY,r.x);resid=max(resid,abs(r.w))
                if let prev=previous,group>0, (group==2 ? r.x < prev-2e-6 : r.x>prev+2e-6){reversed+=1}
                previous=r.x
                if selected.contains(i){
                    picks.append(Pick(group:group,index:i,logK:p.logK,logE:p.logE,rootY:r.x,rootX:r.y,
                        scaledOneGap:r.z,residual:r.w,nearOne:p.nearOne))
                }
            }
        }
        total+=probes.count
        rows.append(Row(group:group,count:probes.count,nonfinite:bad,monotonicReversals:reversed,
                        minimumY:minY,maximumY:maxY,maxAbsResidual:resid,selected:picks))
        print("sample \(utc.string(from:Date())) group=\(group) commands=\(commands) tested=\(total) gpu=\(gpu)");fflush(stdout)
    }
    let obj:[String:Any]=["device":device.name,"evaluated":total,"commands":commands,"perStratum":count,"tile":tile,
      "gpuSeconds":gpu,"wallSeconds":ProcessInfo.processInfo.systemUptime-before,"startUTC":start,
      "endUTC":utc.string(from:Date()),"rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:obj,options:[.sortedKeys,.prettyPrinted])
       .write(to:URL(fileURLWithPath:a[3]),options:.atomic)
}
try run()
