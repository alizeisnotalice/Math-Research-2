import Foundation
import Metal
import Darwin
struct Core: Codable {var a,w,r,length,b,moment,bandMass,valleyMass,logGain:Float}
struct Settings {var shifts,panels:UInt32}
struct Hit: Codable {var profile:Core;var lowerNu,upperNu,lowerResidual,upperResidual:Float}
func grid(_ low:Double,_ high:Double,_ i:Int,_ n:Int)->Double {exp(log(low)+(log(high)-log(low))*Double(i)/Double(n-1))}
func phi(_ x:Double)->Double {exp(-0.5*x*x)/sqrt(2*Double.pi)}
func run() throws {
    let args=CommandLine.arguments
    guard args.count==10 else {throw NSError(domain:"a w r L b shifts panels chunk out",code:1)}
    let n=args[1...8].compactMap(Int.init)
    guard n.count==8,n[0...5].allSatisfy({$0>=2}),n[6]>0,n[7]>0 else {throw NSError(domain:"Dimensions",code:2)}
    guard let d=MTLCreateSystemDefaultDevice(),let queue=d.makeCommandQueue() else {throw NSError(domain:"Metal device",code:3)}
    let source=try String(contentsOfFile:"window.metal",encoding:.utf8)
    let options=MTLCompileOptions();options.mathMode = .safe
    let library=try d.makeLibrary(source:source,options:options)
    let pipeline=try d.makeComputePipelineState(function:library.makeFunction(name:"window_scan")!)
    let threads=min(pipeline.maxTotalThreadsPerThreadgroup,4*pipeline.threadExecutionWidth)
    var settings=Settings(shifts:UInt32(n[5]),panels:UInt32(n[6]))
    var cores=[Core]()
    for ai in 0..<n[0] {
        let a=grid(32,1024,ai,n[0])
        for wi in 0..<n[1] {
            let w=grid(0.0625,1,wi,n[1])
            let pv=erf(w/(a*sqrt(2)))
            for ri in 0..<n[2] {
                let r=grid(0.5,2,ri,n[2])
                for li in 0..<n[3] {
                    let l=grid(0.25,2,li,n[3])
                    let outer=r+l
                    let pb=erf(outer/sqrt(2))-erf(r/sqrt(2))
                    let moment=1+2*(r*phi(r)-outer*phi(outer))/pb
                    for bi in 0..<n[4] {
                        let b=0.125+0.75*Double(bi)/Double(n[4]-1)
                        cores.append(Core(a:Float(a),w:Float(w),r:Float(r),length:Float(l),b:Float(b),moment:Float(moment),
                            bandMass:Float(pb),valleyMass:Float(pv),logGain:Float(log1p(b*pb/pv))))
                    }
                }
            }
        }
    }
    print("device=\(d.name) cores=\(cores.count) shifts=\(n[5]) panels=\(n[6])")
    fflush(stdout)
    var totalTime=0.0;var valid=0;var crosses=0;var negativeFirst=0;var positiveLast=0;var hits=[String:Hit]()
    for start in stride(from:0,to:cores.count,by:n[7]) {
        let end=min(start+n[7],cores.count)
        let input=Array(cores[start..<end])
        let count=input.count*n[5]
        let bu=d.makeBuffer(bytes:input,length:input.count*MemoryLayout<Core>.stride,options:.storageModeShared)!
        let out=d.makeBuffer(length:count*MemoryLayout<Float>.stride,options:.storageModeShared)!
        let cmd=queue.makeCommandBuffer()!
        let enc=cmd.makeComputeCommandEncoder()!
        enc.setComputePipelineState(pipeline);enc.setBuffer(bu,offset:0,index:0)
        enc.setBytes(&settings,length:MemoryLayout<Settings>.stride,index:1);enc.setBuffer(out,offset:0,index:2)
        enc.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        enc.endEncoding();cmd.commit();cmd.waitUntilCompleted()
        guard cmd.status == .completed else {throw cmd.error ?? NSError(domain:"Metal failed",code:4)}
        let elapsed=cmd.gpuEndTime-cmd.gpuStartTime;totalTime+=elapsed
        let values=out.contents().bindMemory(to:Float.self,capacity:count)
        let step:Float=10/Float(n[5]-1)
        for j in 0..<input.count {
            let base=j*n[5]
            if values[base]<0 {negativeFirst+=1}
            if values[base+n[5]-1]>0 {positiveLast+=1}
            var crossing:Int?=nil
            for i in 0..<n[5] {
                if values[base+i].isFinite {valid+=1}
                if i>0 && -2+step*Float(i)>0 && values[base+i-1]>0 && values[base+i]<=0 {crossing=i}
            }
            guard let i=crossing else {continue}
            crosses+=1
            let hit=Hit(profile:input[j],lowerNu:-2+step*Float(i-1),upperNu:-2+step*Float(i),
                lowerResidual:values[base+i-1],upperResidual:values[base+i])
            if hit.upperNu>(hits["latestCrossing"]?.upperNu ?? -.infinity) {hits["latestCrossing"]=hit}
            if hit.upperNu<(hits["earliestCrossing"]?.upperNu ?? .infinity) {hits["earliestCrossing"]=hit}
            if input[j].a>=1024 && hit.upperNu>(hits["latestA1024"]?.upperNu ?? -.infinity) {hits["latestA1024"]=hit}
            if input[j].a>=1024 && hit.upperNu<(hits["earliestA1024"]?.upperNu ?? .infinity) {hits["earliestA1024"]=hit}
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(elapsed) crosses=\(crosses)");fflush(stdout)
    }
    let o:[String:Any]=["device":d.name,"profiles":cores.count,"dimensions":Array(n[0...6]),"evaluated":cores.count*n[5],
        "finite":valid,"crossed":crosses,"negativeAtFirst":negativeFirst,"positiveAtLast":positiveLast,"gpuSeconds":totalTime,
        "extremalCrossings":try JSONSerialization.jsonObject(with:JSONEncoder().encode(hits))]
    let data=try JSONSerialization.data(withJSONObject:o,options:[.prettyPrinted,.sortedKeys])
    try data.write(to:URL(fileURLWithPath:args[9]),options:.atomic)
    print("finite=\(valid) crosses=\(crosses) positive_last=\(positiveLast) gpu_sum_seconds=\(totalTime)")
}
try run()
