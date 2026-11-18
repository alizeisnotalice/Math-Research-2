import Foundation
import Metal
import Darwin
let tag=CommandLine.arguments.dropFirst().first ?? "main"
let ms=tag=="smoke" ? [9,24,128] : [9,12,16,24,32,48,64,96,128,192,256,384,512]
let steps=tag=="smoke" ? 2 : 240
var inputs=[SIMD4<Float>]()
for m in ms {
    for k in 0...1 {
        for i in 0...steps {
            let d=Float(5+10*Double(i)/Double(steps))
            inputs.append(SIMD4(Float(m),d,Float(k),1))
        }
        if m<=128 {
            for d:Float in [5,7.5,10,12.5,15] {
                inputs.append(SIMD4(Float(m),d,Float(k),2))
            }
        }
    }
}
guard let gpu=MTLCreateSystemDefaultDevice(),let queue=gpu.makeCommandQueue() else {
    fatalError("Metal unavailable")
}
let nodes=try Data(contentsOf:URL(fileURLWithPath:"nodes.bin"))
precondition(nodes.count==2*64*8)
let nodeBuf=nodes.withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
}
let spec=inputs.withUnsafeBytes {
    gpu.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
}
let output=gpu.makeBuffer(length:inputs.count*6*4,options:.storageModeShared)!
let source=try String(contentsOfFile:"cluster.metal",encoding:.utf8)
let opts=MTLCompileOptions()
opts.mathMode = .safe
let library=try gpu.makeLibrary(source:source,options:opts)
let pipeline=try gpu.makeComputePipelineState(function:library.makeFunction(name:"cluster")!)
let cb=queue.makeCommandBuffer()!
let enc=cb.makeComputeCommandEncoder()!
enc.setComputePipelineState(pipeline)
enc.setBuffer(spec,offset:0,index:0)
enc.setBuffer(nodeBuf,offset:0,index:1)
enc.setBuffer(output,offset:0,index:2)
enc.dispatchThreads(MTLSize(width:inputs.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:32,height:1,depth:1))
enc.endEncoding()
let utc=ISO8601DateFormatter()
utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
let start=utc.string(from:Date())
print("start=\(start) device=\(gpu.name) pid=\(getpid()) rows=\(inputs.count)")
fflush(stdout)
let tic=ProcessInfo.processInfo.systemUptime
cb.commit()
cb.waitUntilCompleted()
guard cb.status == .completed else {fatalError("Metal failed: \(String(describing:cb.error))")}
let wall=ProcessInfo.processInfo.systemUptime-tic,native=cb.gpuEndTime-cb.gpuStartTime
let end=utc.string(from:Date())
let p=output.contents().bindMemory(to:Float.self,capacity:inputs.count*6)
var rows=[[String:Any]]()
for (i,c) in inputs.enumerated() {
    rows.append(["m":Int(c.x),"D":Double(c.y),"kappa":c.z==0 ? 0.25:0.5,
                 "depth":Int(c.w),"scales":(0..<5).map {Double(p[i*6+$0])},"at0":Double(p[6*i+5])])
}
let result:[String:Any]=["rows":rows,"n":[0,1,2,4,8],"device":gpu.name,
                         "os":ProcessInfo.processInfo.operatingSystemVersionString,
                         "start":start,"end":end,"pid":Int(getpid()),"nativeSeconds":native,
                         "wallSeconds":wall]
try JSONSerialization.data(withJSONObject:result,options:[.prettyPrinted,.sortedKeys])
   .write(to:URL(fileURLWithPath:"report-\(tag).json"))
print("end=\(end) rows=\(inputs.count) gpu_seconds=\(native) wall_seconds=\(wall)")
