import Foundation
import Metal
import Darwin

struct Hit: Codable {
    let a,b,p,q,score,aux: Float
    let ap,am,pp,pm: Float
    let caseId,mode,index,tested,negative,nonfinite: UInt32
}
struct Row: Codable {
    var tested=0
    var negative=0
    var nonfinite=0
    var best: Hit?
    mutating func add(_ row: Hit) {
        tested += Int(row.tested)
        negative += Int(row.negative)
        nonfinite += Int(row.nonfinite)
        if row.score.isFinite && row.score > (best?.score ?? -.infinity) { best=row }
    }
}
struct Report: Codable {
    let device,os,startUTC,endUTC: String
    let profiles,batches,batchSize: Int
    let gridSeconds,commandSeconds: Double
    let byProfile: [[Row]]
    let modes: [String]
}
let modes=["real_pos","real_neg","imag","complex","sub_pos","sub_neg","critical_pos","critical_neg","super"]
let args=CommandLine.arguments
let budget=args.count>1 ? Double(args[1])! : 45
let runs=args.count>2 ? Int(args[2])! : 32
let tag=args.count>3 ? args[3] : "main"
let data=try Data(contentsOf:URL(fileURLWithPath:"profiles.bin"))
let sizes=data.withUnsafeBytes { raw in
    (raw.loadUnaligned(fromByteOffset:0,as:UInt32.self),
     raw.loadUnaligned(fromByteOffset:4,as:UInt32.self))
}
let profiles=Int(sizes.0),stride=Int(sizes.1)+5
precondition(data.count==8+4*profiles*stride)
precondition(MemoryLayout<Hit>.stride==64)
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Apple Metal device unavailable")
}
let input=data.withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!.advanced(by:8),length:data.count-8,options:.storageModeShared)!
}
let source=try String(contentsOfFile:"scan.metal",encoding:.utf8)
let options=MTLCompileOptions()
options.mathMode = .safe
let lib=try gpu.makeLibrary(source:source,options:options)
let pipeline=try gpu.makeComputePipelineState(function:lib.makeFunction(name:"validation")!)
precondition(pipeline.threadExecutionWidth==32)
let nodes: [SIMD2<Float>] = [
    SIMD2(0.01985507175,0.05061426815),SIMD2(0.10166676129,0.11119051723),
    SIMD2(0.23723379504,0.15685332294),SIMD2(0.40828267875,0.18134189169),
    SIMD2(0.59171732125,0.18134189169),SIMD2(0.76276620496,0.15685332294),
    SIMD2(0.89833323871,0.11119051723),SIMD2(0.98014492825,0.05061426815),
]
let nodeBuf=nodes.withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!, length:$0.count,options:.storageModeShared)!
}
let groups=profiles*9*runs
let batchSize=groups*32
let output=gpu.makeBuffer(length:64*groups,options:.storageModeShared)!
let out=output.contents().bindMemory(to:Hit.self,capacity:groups)
let formatter=ISO8601DateFormatter()
formatter.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
if args.count>4, let epoch=Double(args[4]), epoch>0 {
    let remaining=epoch-Date().timeIntervalSince1970
    if remaining>0 { Thread.sleep(forTimeInterval:remaining) }
}
let begin=formatter.string(from:Date())
print("start=\(begin) device=\(gpu.name) pid=\(getpid()) profiles=\(profiles) rows=\(batchSize)")
fflush(stdout)
var stat=Array(repeating:Row(),count:profiles*9)
let start=ProcessInfo.processInfo.systemUptime
var native=0.0
var batches=0
repeat {
    var p=SIMD4(sizes.0,UInt32(truncatingIfNeeded:batches * batchSize),UInt32(stride),UInt32(runs))
    let cb=queue.makeCommandBuffer()!
    let encoder=cb.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(input,offset:0,index:0)
    encoder.setBuffer(nodeBuf,offset:0,index:1)
    encoder.setBuffer(output,offset:0,index:2)
    encoder.setBytes(&p,length:MemoryLayout<SIMD4<UInt32>>.stride,index:3)
    encoder.dispatchThreads(MTLSize(width:batchSize,height:1,depth:1),
                            threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1))
    encoder.endEncoding()
    cb.commit()
    cb.waitUntilCompleted()
    guard cb.status == .completed else { fatalError("Metal command failed: \(String(describing:cb.error))") }
    native += cb.gpuEndTime - cb.gpuStartTime
    for j in 0..<groups {
        stat[j/runs].add(out[j])
    }
    batches+=1
    if batches==1 || batches%25==0 {
        print("utc=\(formatter.string(from:Date())) batches=\(batches) tested=\(batches*batchSize) gpu_seconds=\(native)")
        fflush(stdout)
    }
} while ProcessInfo.processInfo.systemUptime-start<budget
let elapsed=ProcessInfo.processInfo.systemUptime-start
let end=formatter.string(from:Date())
precondition(stat.allSatisfy {$0.tested==batches*runs*32})
let byProfile=(0..<profiles).map {Array(stat[$0*9..<$0*9+9])}
let report=Report(device:gpu.name,os:ProcessInfo.processInfo.operatingSystemVersionString,
                  startUTC:begin,endUTC:end,profiles:profiles,batches:batches,batchSize:batchSize,
                  gridSeconds:elapsed,commandSeconds:native,byProfile:byProfile,modes:modes)
let json=JSONEncoder()
json.outputFormatting=[.prettyPrinted,.sortedKeys]
try json.encode(report).write(to:URL(fileURLWithPath:"scan-\(tag).json"),options:.atomic)
print("end=\(end) batches=\(batches) tested=\(batches*batchSize) gpu_seconds=\(native) wall_seconds=\(elapsed)")
