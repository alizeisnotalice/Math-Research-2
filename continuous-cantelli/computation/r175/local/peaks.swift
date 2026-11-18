import Foundation
import Metal
import Darwin

struct Shape: Codable {
    let c,L,s,w1,w2,h1,h2,l1,r1,l2,r2:Float
}
struct Row: Codable {
    let shape:Shape
    let D0,D1,D2,parab,det,sine,score,maxValue:Float
    let mode,index,tested,witness,bad,salt:UInt32
}
struct Ledger: Codable {
    var tested=0
    var bad=0
    var witness=0
    var best: Row?
    mutating func add(_ row:Row) {
        tested+=Int(row.tested)
        bad+=Int(row.bad)
        witness+=Int(row.witness)
        if row.score.isFinite && row.score > (best?.score ?? -.infinity) { best=row }
    }
}
struct Result:Codable {
    let startUTC,endUTC,device,os:String
    let batchSize,batches,order:Int
    let commandSeconds,wallSeconds:Double
    let strata:[Ledger]
}
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
            let change=value/derivative
            x-=change
            if abs(change)<2e-16 {break}
        }
        let w=Float(1/((1-x*x)*derivative*derivative))
        nodes[i]=SIMD2(Float((1-x)/2),w)
        nodes[order-1-i]=SIMD2(Float((1+x)/2),w)
    }
    return nodes
}
let args=CommandLine.arguments
let seconds=args.count>1 ? Double(args[1])! : 45
let order=args.count>2 ? Int(args[2])! : 72
let tag=args.count>3 ? args[3] : "main"
let batch=1<<17,groups=batch/32
precondition(MemoryLayout<Row>.stride==100 && order%2==0)
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Metal device unavailable")
}
let source=try String(contentsOfFile:"peaks.metal",encoding:.utf8)
let options=MTLCompileOptions()
options.mathMode = .safe
let lib=try gpu.makeLibrary(source:source,options:options)
let pipeline=try gpu.makeComputePipelineState(function:lib.makeFunction(name:"peaks_scan")!)
precondition(pipeline.threadExecutionWidth==32)
let nodes=gaussLegendre(order).withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
}
let resultBuffer=gpu.makeBuffer(length:100*groups,options:.storageModeShared)!
let rows=resultBuffer.contents().bindMemory(to:Row.self,capacity:groups)
let formatter=ISO8601DateFormatter()
formatter.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
let startUTC=formatter.string(from:Date())
print("start=\(startUTC) device=\(gpu.name) pid=\(getpid()) order=\(order) batch=\(batch)")
fflush(stdout)
let start=ProcessInfo.processInfo.systemUptime
var ledgers=Array(repeating:Ledger(),count:6)
var count=0
var timed=0.0
repeat {
    let base=UInt64(count)*UInt64(batch)
    var p=SIMD4(UInt32(batch),UInt32(order),UInt32(truncatingIfNeeded:base),UInt32(base>>32))
    let cb=queue.makeCommandBuffer()!
    let enc=cb.makeComputeCommandEncoder()!
    enc.setComputePipelineState(pipeline)
    enc.setBuffer(nodes,offset:0,index:0)
    enc.setBuffer(resultBuffer,offset:0,index:1)
    enc.setBytes(&p,length:MemoryLayout<SIMD4<UInt32>>.stride,index:2)
    enc.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1))
    enc.endEncoding()
    cb.commit()
    cb.waitUntilCompleted()
    guard cb.status == .completed else {fatalError("Metal command failed: \(String(describing:cb.error))")}
    timed+=cb.gpuEndTime-cb.gpuStartTime
    for r in 0..<groups {
        let row=rows[r]
        precondition(row.mode<6)
        ledgers[Int(row.mode)].add(row)
    }
    count+=1
    if count==1 || count%100==0 {
        print("utc=\(formatter.string(from:Date())) batches=\(count) tested=\(count*batch) gpu_seconds=\(timed)")
        fflush(stdout)
    }
} while ProcessInfo.processInfo.systemUptime-start<seconds
let elapsed=ProcessInfo.processInfo.systemUptime-start
let endUTC=formatter.string(from:Date())
precondition(ledgers.reduce(0){$0+$1.tested}==count*batch)
let output=Result(startUTC:startUTC,endUTC:endUTC,device:gpu.name,
                  os:ProcessInfo.processInfo.operatingSystemVersionString,batchSize:batch,
                  batches:count,order:order,commandSeconds:timed,wallSeconds:elapsed,strata:ledgers)
let json=JSONEncoder()
json.outputFormatting=[.prettyPrinted,.sortedKeys]
try json.encode(output).write(to:URL(fileURLWithPath:"peaks-\(tag).json"),options:.atomic)
print("end=\(endUTC) tested=\(count*batch) gpu_seconds=\(timed) wall_seconds=\(elapsed)")
