import Foundation
import Metal
import Darwin

func gaussLegendre(_ order:Int)->[SIMD2<Float>] {
    var nodes=[SIMD2<Float>](repeating: SIMD2.zero,count:order)
    for i in 0..<(order/2) {
        var x=cos(Double.pi*(Double(i)+0.75)/(Double(order)+0.5))
        var derivative=0.0
        for _ in 0..<12 {
            var previous=0.0
            var value=1.0
            for j in 1...order {
                let earlier=previous
                previous=value
                value=((2*Double(j)-1)*x*previous-(Double(j)-1)*earlier)/Double(j)
            }
            derivative=Double(order)*(x*value-previous)/(x*x-1)
            x-=value/derivative
        }
        let w=Float(1/((1-x*x)*derivative*derivative))
        nodes[i]=SIMD2(Float(-x/2),w)
        nodes[order-1-i]=SIMD2(Float(x/2),w)
    }
    return nodes
}
let smoke=CommandLine.arguments.dropFirst().first=="smoke"
let domainStages=[0,1,2,4,6,8,10,12]
let centerStages=domainStages+[14,16]
var centers=[0.0]
var gaps=[SIMD2<Float>]()
var width=1.0
for stage in 1...16 {
    let theta=pow(2.0,-Double(stage))/256
    let gap=width*theta
    gaps+=centers.map {SIMD2(Float($0),Float(gap))}
    let newWidth=(width-gap)/2
    let offset=(width-newWidth)/2
    centers=centers.flatMap {[$0-offset,$0+offset]}
    width=newWidth
}
var queries=[SIMD4<Float>]()
let coarseV=smoke ? 5 : 41
let coarseX=smoke ? 17 : 257
let fineV=smoke ? 9 : 513
func grid(_ i:Int,_ count:Int,_ lo:Double,_ hi:Double)->Float {
    Float(lo+(hi-lo)*Double(i)/Double(count-1))
}
for j in domainStages {
    for v in 0..<coarseV {
        for x in 0..<coarseX {
            queries.append(SIMD4(grid(x,coarseX,-0.64,0.64),grid(v,coarseV,0.125,1.0/6),Float(j),4))
        }
    }
}
for j in centerStages {
    for order in [1,4] {
        for v in 0..<fineV {
            queries.append(SIMD4(0,grid(v,fineV,0.14,0.16),Float(j),Float(order)))
        }
    }
}
let nodes=gaussLegendre(16)+gaussLegendre(4)
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Metal device unavailable")
}
func make<T>(_ array:[T])->MTLBuffer {
    array.withUnsafeBytes {gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
}
let qbuf=make(queries),gbuf=make(gaps),nbuf=make(nodes)
let result=gpu.makeBuffer(length:queries.count*6*MemoryLayout<Float>.size,options:.storageModeShared)!
let source=try String(contentsOfFile:"cantor.metal",encoding:.utf8)
let opts=MTLCompileOptions()
opts.mathMode = .safe
let library=try gpu.makeLibrary(source:source,options:opts)
let pipeline=try gpu.makeComputePipelineState(function:library.makeFunction(name:"cantor")!)
let command=queue.makeCommandBuffer()!
let encoder=command.makeComputeCommandEncoder()!
encoder.setComputePipelineState(pipeline)
encoder.setBuffer(qbuf,offset:0,index:0)
encoder.setBuffer(gbuf,offset:0,index:1)
encoder.setBuffer(nbuf,offset:0,index:2)
encoder.setBuffer(result,offset:0,index:3)
encoder.dispatchThreads(MTLSize(width:queries.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:128,height:1,depth:1))
encoder.endEncoding()
let tag=smoke ? "smoke" : "main"
let utc=ISO8601DateFormatter()
utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
let start=utc.string(from:Date())
print("start=\(start) device=\(gpu.name) pid=\(getpid()) rows=\(queries.count)")
fflush(stdout)
let before=ProcessInfo.processInfo.systemUptime
command.commit()
command.waitUntilCompleted()
guard command.status == .completed else {fatalError("Metal failed: \(String(describing:command.error))")}
let gpuSeconds=command.gpuEndTime-command.gpuStartTime
let wallSeconds=ProcessInfo.processInfo.systemUptime-before
let end=utc.string(from:Date())
let out=Data(bytes:result.contents(),count:result.length)
try out.write(to:URL(fileURLWithPath:"values-\(tag).f32"))
try queries.withUnsafeBytes { try Data($0).write(to:URL(fileURLWithPath:"queries-\(tag).f32")) }
let receipt:[String:Any]=["device":gpu.name,"os":ProcessInfo.processInfo.operatingSystemVersionString,
                         "startUTC":start,"endUTC":end,"rows":queries.count,"nativeSeconds":gpuSeconds,
                         "wallSeconds":wallSeconds,"pid":getpid(),"coarseV":coarseV,
                         "coarseX":coarseX,"fineV":fineV,"domainStages":domainStages,"centerStages":centerStages]
try JSONSerialization.data(withJSONObject:receipt,options:[.prettyPrinted,.sortedKeys])
    .write(to:URL(fileURLWithPath:"receipt-\(tag).json"))
print("end=\(end) tested=\(queries.count) gpu_seconds=\(gpuSeconds) wall_seconds=\(wallSeconds)")
