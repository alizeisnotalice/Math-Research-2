import Foundation
import Metal
import Darwin

struct Row: Codable {
    let eps,turns,sine,rawSine,first,second,stableScaled,rawScaled: Float
    let secondReal,secondImag,mass2,quadFloor,score: Float
    let n,mode,band,tested,index,bad,salt: UInt32
}
struct Ledger: Codable {
    var tested=0
    var bad=0
    var best: Row?
    mutating func add(_ row:Row) {
        tested+=Int(row.tested)
        bad+=Int(row.bad)
        if row.score.isFinite && row.score > (best?.score ?? -.infinity) { best=row }
    }
}
struct Result:Codable {
    let startUTC,endUTC,device,os:String
    let batchSize,batches,order:Int
    let commandSeconds,wallSeconds:Double
    let metrics:[String]
    let strata:[[[Ledger]]]
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
let batch=1<<16,groups=batch/32
precondition(MemoryLayout<Row>.stride==80 && order%2==0)
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Metal device unavailable")
}
let source=try String(contentsOfFile:"wells.metal",encoding:.utf8)
let options=MTLCompileOptions()
options.mathMode = .safe
let lib=try gpu.makeLibrary(source:source,options:options)
let pipeline=try gpu.makeComputePipelineState(function:lib.makeFunction(name:"wells_scan")!)
precondition(pipeline.threadExecutionWidth==32)
let nodes=gaussLegendre(order).withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
}
let resultBuffer=gpu.makeBuffer(length:80*groups*3,options:.storageModeShared)!
let rows=resultBuffer.contents().bindMemory(to:Row.self,capacity:groups*3)
let formatter=ISO8601DateFormatter()
formatter.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
let startUTC=formatter.string(from:Date())
print("start=\(startUTC) device=\(gpu.name) pid=\(getpid()) order=\(order) batch=\(batch)")
fflush(stdout)
let start=ProcessInfo.processInfo.systemUptime
var ledgers=Array(repeating:Ledger(),count:4*5*3)
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
    for r in 0..<groups*3 {
        let row=rows[r]
        precondition(row.n>=16 && row.mode<5 && row.band<4)
        ledgers[(Int(row.band)*5+Int(row.mode))*3+r%3].add(row)
    }
    count+=1
    if count==1 || count%100==0 {
        print("utc=\(formatter.string(from:Date())) batches=\(count) tested=\(count*batch) gpu_seconds=\(timed)")
        fflush(stdout)
    }
} while ProcessInfo.processInfo.systemUptime-start<seconds
let elapsed=ProcessInfo.processInfo.systemUptime-start
let endUTC=formatter.string(from:Date())
precondition(ledgers.enumerated().filter {$0.offset%3==0}.reduce(0){$0+$1.element.tested}==count*batch)
let table=(0..<4).map { band in
    (0..<5).map { mode in Array(ledgers[(band*5+mode)*3..<(band*5+mode)*3+3]) }
}
let output=Result(startUTC:startUTC,endUTC:endUTC,device:gpu.name,
                  os:ProcessInfo.processInfo.operatingSystemVersionString,batchSize:batch,
                  batches:count,order:order,commandSeconds:timed,wallSeconds:elapsed,
                  metrics:["raw_phase_error","raw_scaled_error","second_remainder_error"],strata:table)
let json=JSONEncoder()
json.outputFormatting=[.prettyPrinted,.sortedKeys]
try json.encode(output).write(to:URL(fileURLWithPath:"wells-\(tag).json"),options:.atomic)
print("end=\(endUTC) tested=\(count*batch) gpu_seconds=\(timed) wall_seconds=\(elapsed)")
