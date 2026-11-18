import Foundation
import Metal
import Darwin
struct Core: Codable {var k,l,coefficient,lambda,center,mean,std:Float}
struct Settings {var panels,strikes:UInt32}
struct Hit:Codable {var core:Core;var signChanges:Int;var maxPositive,maxNegative,positiveStrike,negativeStrike:Float}
func profile(k:Double,l:Double,c:Double,margin:Double)->Core?{
    let rk=exp(-0.5*k*k),rl=exp(-0.5*l*l)
    let v1=(1-rk*rk)*(1-rk*rk)/2
    let v2=(1-rl*rl)*(1-rl*rl)/2
    let cov=(exp(-0.5*(k-l)*(k-l))+exp(-0.5*(k+l)*(k+l)))/2-rk*rl
    let variance=v1+c*c*v2+2*c*cov
    guard variance>1e-8 else{return nil}
    let slope=2*(k*k*rk+c*l*l*rl)/variance
    let sigma=abs(slope)*sqrt(variance)
    let center=rk+c*rl
    let bound=abs(slope)*(1+abs(c))
    guard sigma>=0.004,bound-slope*center<=1 else{return nil}
    let mean=bound+slope*center+margin
    return Core(k:Float(k),l:Float(l),coefficient:Float(c),lambda:Float(slope),center:Float(center),mean:Float(mean),std:Float(sigma))
}
func run()throws{
    let a=CommandLine.arguments
    guard a.count==9 else{throw NSError(domain:"k count, l count, eta count, margins, strikes, panels, batch, output",code:1)}
    let n=a[1...7].compactMap(Int.init)
    guard n.count==7,n[0...4].allSatisfy({$0>=2}),n[5]>0,n[6]>0 else{throw NSError(domain:"dimensions",code:2)}
    var grid=[Core]()
    for i in 0..<n[0]{
        let k=2.5+4.5*Double(i)/Double(n[0]-1)
        for j in 0..<n[1]{
            let l=k+0.15+(9-k-0.15)*Double(j)/Double(n[1]-1)
            for h in 0..<n[2]{
                let c = -2+4*Double(h)/Double(n[2]-1)
                for m in 0..<n[3]{
                    let margin=0.35*pow(6/0.35,Double(m)/Double(n[3]-1))
                    if let p=profile(k:k,l:l,c:c,margin:margin){grid.append(p)}
                }
            }
        }
    }
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else{throw NSError(domain:"Metal device",code:3)}
    let src=try String(contentsOfFile:"call.metal",encoding:.utf8)
    let option=MTLCompileOptions();option.mathMode = .safe
    let lib=try device.makeLibrary(source:src,options:option)
    let state=try device.makeComputePipelineState(function:lib.makeFunction(name:"scan_calls")!)
    let threads=min(state.maxTotalThreadsPerThreadgroup,state.threadExecutionWidth*4)
    var cfg=Settings(panels:UInt32(n[5]),strikes:UInt32(n[4]))
    var elapsed=0.0,finite=0
    var signCounts=[String:Int](), hits=[String:Hit]()
    print("device=\(device.name) cores=\(grid.count) strikes=\(n[4]) panels=\(n[5])");fflush(stdout)
    for start in stride(from:0,to:grid.count,by:n[6]){
        let end=min(grid.count,start+n[6]),tile=Array(grid[start..<end]),count=tile.count*n[4]
        let input=device.makeBuffer(bytes:tile,length:tile.count*MemoryLayout<Core>.stride,options:.storageModeShared)!
        let output=device.makeBuffer(length:count*MemoryLayout<Float>.stride,options:.storageModeShared)!
        let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
        e.setComputePipelineState(state);e.setBuffer(input,offset:0,index:0);e.setBuffer(output,offset:0,index:2)
        e.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:1)
        e.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        e.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed else{throw cb.error ?? NSError(domain:"GPU failed",code:4)}
        let gpu=cb.gpuEndTime-cb.gpuStartTime;elapsed+=gpu
        let values=output.contents().bindMemory(to:Float.self,capacity:count)
        for j in 0..<tile.count{
            let start=j*n[4]
            var hi:Float=0,lo:Float=0
            var ihi=0,ilo=0
            for i in 0..<n[4]{
                let v=values[start+i]
                if !v.isFinite {continue}
                finite+=1
                if v>hi {hi=v;ihi=i}
                if v<lo {lo=v;ilo=i}
            }
            let floor=max(1e-6*max(hi,-lo),1e-12)
            var sign=0,changes=0
            for i in 0..<n[4]{
                let v=values[start+i]
                if abs(v)<floor {continue}
                let s=v>0 ? 1 : -1
                if sign != 0 && s != sign {changes+=1}
                sign=s
            }
            signCounts[String(changes),default:0]+=1
            let f:Float=7/Float(n[4]-1)
            let hit=Hit(core:tile[j],signChanges:changes,maxPositive:hi,maxNegative:lo,positiveStrike:Float(ihi)*f,negativeStrike:Float(ilo)*f)
            if hi>(hits["maxPositive"]?.maxPositive ?? -.infinity){hits["maxPositive"]=hit}
            if lo<(hits["maxNegative"]?.maxNegative ?? .infinity){hits["maxNegative"]=hit}
            if changes>(hits["mostSignChanges"]?.signChanges ?? -1){hits["mostSignChanges"]=hit}
            let strength=max(hi,-lo)/tile[j].std
            let old=hits["leastScaledResidual"]
            if old==nil || strength<max(old!.maxPositive,-old!.maxNegative)/old!.core.std {hits["leastScaledResidual"]=hit}
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(gpu)");fflush(stdout)
    }
    let values:[String:Any]=["device":device.name,"finite":finite,"profiles":grid.count,"evaluated":grid.count*n[4],"gpuSeconds":elapsed,
        "dimensions":Array(n[0...5]),"signChangeCounts":signCounts,"representatives":try JSONSerialization.jsonObject(with:JSONEncoder().encode(hits))]
    try JSONSerialization.data(withJSONObject:values,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[8]),options:.atomic)
    print("evaluated=\(grid.count*n[4]) finite=\(finite) gpu_seconds=\(elapsed)")
}
try run()
