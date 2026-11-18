import Darwin
import Foundation
import Metal

struct Hermite:Codable {var u,re,im:Double}
struct Node {var u,re,im,pad:Float}
struct Config {var step,width,offset,scaleM:Float;var m,mid,order,pad:UInt32}
struct Result:Codable {
    var beta,m,j,nodes:Int
    var c,A,L,b,D2,s,xi,offset,offsetOverW,step,W,I,I0,ratio,centerPower,zeroPower,logM,exponent,asymptotic:Double
}
func run() throws {
    let a=CommandLine.arguments
    guard a.count==4,let maxLevel=Int(a[1]),let baseStep=Double(a[2]),maxLevel>=1,maxLevel<=12,baseStep>0,
          let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue()
    else {throw NSError(domain:"maxLevel step output",code:1)}
    let hs=try JSONDecoder().decode([Hermite].self,from:Data(contentsOf:URL(fileURLWithPath:"hermite.json")))
    let weights=hs.map {Node(u:Float($0.u),re:Float($0.re),im:Float($0.im),pad:0)}
    var shrinks=[SIMD2<Float>](),rot=[SIMD2<Float>]()
    for j in 1...13 {
        let r=pow(3.0,-Double(j)),hi=Float(r)
        shrinks.append(SIMD2<Float>(hi,Float(r-Double(hi))))
    }
    for h in hs{
        for j in 1...13{
            let v=h.u*pow(3.0,-Double(j))
            rot.append(SIMD2<Float>(Float(cos(v)),Float(sin(v))))
        }
    }
    let dWeights=device.makeBuffer(bytes:weights,length:weights.count*MemoryLayout<Node>.stride,options:.storageModeShared)!
    let dShrinks=device.makeBuffer(bytes:shrinks,length:shrinks.count*MemoryLayout<SIMD2<Float>>.stride,options:.storageModeShared)!
    let dRot=device.makeBuffer(bytes:rot,length:rot.count*MemoryLayout<SIMD2<Float>>.stride,options:.storageModeShared)!
    let opts=MTLCompileOptions();opts.mathMode = .safe
    let source=try String(contentsOfFile:"window.metal",encoding:.utf8)
    let lib=try device.makeLibrary(source:source,options:opts)
    let pipe=try device.makeComputePipelineState(function:lib.makeFunction(name:"cantor_window")!)
    let threads=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let utc=ISO8601DateFormatter();utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    let start=ProcessInfo.processInfo.systemUptime,startUTC=utc.string(from:Date())
    let d=log(2.0)/log(3.0),critical=1-d
    var buffers=0,gpu=0.0,evaluated=0
    var rows=[Result]()
    func perform(_ beta:Int,_ m:Int,_ c:Double,_ A:Double,_ rounded:Bool,_ step:Double)throws {
        let B=Double(beta),sc=pow(3.0,Double(m))
        let xi=2*Double.pi*sc
        let exact=exp((2*log(xi)+log(B*A/(2*c)))/(B+1))
        let j=rounded ? max(1,Int(round(exact/log(3.0)))):0
        let L=rounded ? Double(j)*log(3.0):exact
        let off=rounded ? xi*expm1((B+1)/2*log(L/exact)):0
        let b=B*A/pow(L,B),W=1/sqrt(b)
        let mid=Int(ceil((8*W+abs(off))/step))
        let size=2*mid+1
        guard let output=device.makeBuffer(length:size*MemoryLayout<SIMD4<Float>>.stride,options:.storageModeShared)
        else {throw NSError(domain:"Metal output allocation",code:3)}
        var cfg=Config(step:Float(step),width:Float(W),offset:Float(off),scaleM:Float(1/sc),
                       m:UInt32(m),mid:UInt32(mid),order:UInt32(hs.count),pad:0)
        let cb=queue.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
        e.setComputePipelineState(pipe)
        e.setBuffer(dWeights,offset:0,index:0);e.setBuffer(dShrinks,offset:0,index:1)
        e.setBuffer(dRot,offset:0,index:2);e.setBytes(&cfg,length:MemoryLayout<Config>.stride,index:3)
        e.setBuffer(output,offset:0,index:4)
        e.dispatchThreads(MTLSize(width:size,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        e.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed,cb.gpuEndTime>cb.gpuStartTime
        else {throw cb.error ?? NSError(domain:"Metal clock or native execution",code:4)}
        gpu+=cb.gpuEndTime-cb.gpuStartTime;buffers+=1;evaluated+=size
        let p=output.contents().bindMemory(to:SIMD4<Float>.self,capacity:size)
        var total=0.0,total0=0.0
        for k in 0..<size {total+=Double(p[k].x);total0+=Double(p[k].y)}
        let I=total*step/(2*Double.pi),I0=total0*step/(2*Double.pi)
        guard I>0,I0>0,I.isFinite,I0.isFinite else {throw NSError(domain:"nonfinite energy",code:5)}
        let logM=c*L+log(I)
        rows.append(Result(beta:beta,m:m,j:j,nodes:size,c:c,A:A,L:L,b:b,D2:2*c*b*L,s:(B-1)*A/pow(L,B),xi:xi,
          offset:off,offsetOverW:off/W,step:step,W:W,I:I,I0:I0,ratio:I/I0,centerPower:Double(p[mid].z),
          zeroPower:Double(p[mid].w),logM:logM,exponent:logM/L-critical,asymptotic:c-critical))
    }
    print("start \(startUTC) pid=\(getpid()) device=\(device.name) maxLevel=\(maxLevel) step=\(baseStep)")
    fflush(stdout)
    for beta in 1...3{
        for m in 1...maxLevel{
            for c in [0.25,critical,0.55]{
                for A in [0.5,2.0]{
                    try perform(beta,m,c,A,false,baseStep)
                    try perform(beta,m,c,A,true,baseStep)
                }
            }
            if [1,4,7,10].contains(m) {try perform(beta,m,critical,0.5,false,baseStep/2)}
            print("sample \(utc.string(from:Date())) beta=\(beta) m=\(m) commands=\(buffers) frequencies=\(evaluated) gpu=\(gpu)")
            fflush(stdout)
        }
    }
    let obj:[String:Any]=["device":device.name,"commands":buffers,"frequencies":evaluated,"gpuSeconds":gpu,
        "wallSeconds":ProcessInfo.processInfo.systemUptime-start,"startUTC":startUTC,"endUTC":utc.string(from:Date()),
        "rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:obj,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[3]),options:.atomic)
}
try run()
