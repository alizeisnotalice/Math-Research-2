import Foundation
import Metal
import Darwin

struct Case: Codable {
    let m,j,n:Int
    let sigma,kappa:Double
    let beta,theoremC:Double?
    let energy,midpoint,gap:[Double]
}
struct Report: Codable {
    let device,os,startUTC,endUTC:String
    let cases:[Case]
    let order,rows,pid:Int
    let nativeSeconds,wallSeconds:Double
}
func gaussLegendre(_ order:Int)->[SIMD2<Float>] {
    var nodes=[SIMD2<Float>](repeating:SIMD2.zero,count:order)
    for i in 0..<(order/2) {
        var x=cos(Double.pi*(Double(i)+0.75)/(Double(order)+0.5)),d=0.0
        for _ in 0..<10 {
            var p=0.0,q=1.0
            for k in 1...order {
                let r=p
                p=q
                q=((2*Double(k)-1)*x*p-(Double(k)-1)*r)/Double(k)
            }
            d=Double(order)*(x*q-p)/(x*x-1)
            x-=q/d
        }
        let w=1/((1-x*x)*d*d)
        nodes[i]=SIMD2(Float(-x/2),Float(w))
        nodes[order-1-i]=SIMD2(Float(x/2),Float(w))
    }
    return nodes
}
let tag=CommandLine.arguments.dropFirst().first ?? "main"
let order=tag=="smoke" ? 8 : 16
let maxima=tag=="smoke" ? [0,2,5] : [0,1,2,3,5,7,9]
let kappas=tag=="smoke" ? [0.42,1.4] : [0.20,0.42,0.8,1.4,2.5,3.5]
var pts:[UInt32]=[0]
var layer:[UInt32]=[0]
for _ in 1...14 {
    layer=layer.flatMap {[$0*3,$0*3+2]}
    pts+=layer
}
var tasks=[SIMD4<UInt32>]()
var sigmas=[Float]()
var spans=[(m:Int,j:Int,n:Int,k:Double,s:Float,lo:Int,hi:Int,beta:Double?,c:Double?)]()
func appendCase(_ m:Int,_ j:Int,_ n:Int,_ k:Double,_ beta:Double?=nil,_ c:Double?=nil) {
    let s=Float(k*pow(3.0,-Double(m)))
    let start=tasks.count
    for i in 0..<(1<<n) {
        tasks.append(SIMD4(UInt32(n),UInt32(j),UInt32(m),UInt32(i)))
        sigmas.append(s)
    }
    spans.append((m,j,n,k,s,start,tasks.count,beta,c))
}
if tag=="theory" {
    let dimension=log(2.0)/log(3.0)
    for beta in [2.0,3.0,4.0] {
        for m in 0...3 {
            let omega=2*Double.pi*pow(3.0,Double(m))
            for j in max(2,m+1)...11 {
                let L=Double(j)*log(3.0)
                let b=(beta-1)*pow(L,-beta)+pow(L+log(6.0),-beta)
                let c=b*omega*omega/(2*L)
                if c>1-dimension && c<1 {
                    appendCase(m,j,max(m+6,j+2),sqrt(b)*pow(3.0,Double(m)),beta,c)
                }
            }
        }
    }
} else {
    for m in maxima {
        for j in Set([-1,0,1,2,3,4].map {min(12,max(0,m+$0))}).sorted() {
            for k in kappas {appendCase(m,j,min(14,m+6),k)}
        }
    }
}
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Metal GPU unavailable")
}
func make<T>(_ array:[T])->MTLBuffer {
    array.withUnsafeBytes {gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
}
let qbuf=make(tasks),sbuf=make(sigmas),pbuf=make(pts),nbuf=make(gaussLegendre(order))
let result=gpu.makeBuffer(length:tasks.count*6*4,options:.storageModeShared)!
let source=try String(contentsOfFile:"cantor.metal",encoding:.utf8)
let opts=MTLCompileOptions()
opts.mathMode = .safe
let library=try gpu.makeLibrary(source:source,options:opts)
let pipeline=try gpu.makeComputePipelineState(function:library.makeFunction(name:"cantor")!)
let command=queue.makeCommandBuffer()!
let encoder=command.makeComputeCommandEncoder()!
encoder.setComputePipelineState(pipeline)
encoder.setBuffer(qbuf,offset:0,index:0)
encoder.setBuffer(sbuf,offset:0,index:1)
encoder.setBuffer(pbuf,offset:0,index:2)
encoder.setBuffer(nbuf,offset:0,index:3)
var nOrder=UInt32(order)
encoder.setBytes(&nOrder,length:4,index:4)
encoder.setBuffer(result,offset:0,index:5)
encoder.dispatchThreads(MTLSize(width:tasks.count,height:1,depth:1),
                        threadsPerThreadgroup:MTLSize(width:128,height:1,depth:1))
encoder.endEncoding()
let utc=ISO8601DateFormatter()
utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
let start=utc.string(from:Date())
print("start=\(start) device=\(gpu.name) pid=\(getpid()) cases=\(spans.count) rows=\(tasks.count) order=\(order)")
fflush(stdout)
let begin=ProcessInfo.processInfo.systemUptime
command.commit()
command.waitUntilCompleted()
guard command.status == .completed else {fatalError("Metal failed: \(String(describing:command.error))")}
let seconds=command.gpuEndTime-command.gpuStartTime
let wall=ProcessInfo.processInfo.systemUptime-begin
let end=utc.string(from:Date())
let rows=result.contents().bindMemory(to:Float.self,capacity:tasks.count*6)
var cases=[Case]()
for c in spans {
    var s=Array(repeating:0.0,count:6)
    for i in c.lo..<c.hi {
        for k in 0..<6 {s[k]+=Double(rows[6*i+k])}
    }
    s=s.map {$0/Double(c.hi-c.lo)}
    cases.append(Case(m:c.m,j:c.j,n:c.n,sigma:Double(c.s),kappa:c.k,
                      beta:c.beta,theoremC:c.c,energy:Array(s[0..<2]),midpoint:Array(s[2..<4]),gap:Array(s[4..<6])))
}
let report=Report(device:gpu.name,os:ProcessInfo.processInfo.operatingSystemVersionString,
                  startUTC:start,endUTC:end,cases:cases,order:order,rows:tasks.count,
                  pid:Int(getpid()),nativeSeconds:seconds,wallSeconds:wall)
let json=JSONEncoder()
json.outputFormatting=[.prettyPrinted,.sortedKeys]
try json.encode(report).write(to:URL(fileURLWithPath:"report-\(tag).json"))
print("end=\(end) rows=\(tasks.count) gpu_seconds=\(seconds) wall_seconds=\(wall)")
