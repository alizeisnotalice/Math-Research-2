import Foundation
import Metal
import Darwin

struct Core: Codable { var center, width, compression, drift: Float }
struct Settings { var offsets, halfWidths, panels: UInt32 }
struct Hit: Codable { var core: Core; var halfWidth, offset, ratio: Float }

func logGrid(_ low: Double, _ high: Double, _ i: Int, _ count: Int) -> Double {
    exp(log(low) + (log(high) - log(low)) * Double(i) / Double(count - 1))
}

func run() throws {
    let args=CommandLine.arguments
    guard args.count==10 else { throw NSError(domain:"Usage: a w compression drift offset halfwidth panels chunk out",code:1) }
    let n=args[1...8].compactMap(Int.init)
    guard n.count==8,n[0...5].allSatisfy({$0>=2}),n[6]>0,n[7]>0 else {throw NSError(domain:"Dimensions",code:2)}
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else {throw NSError(domain:"Metal device",code:3)}
    let source=try String(contentsOfFile:"dual.metal",encoding:.utf8)
    let compiler=MTLCompileOptions()
    compiler.mathMode = .safe
    let library=try device.makeLibrary(source:source,options:compiler)
    let pipeline=try device.makeComputePipelineState(function:library.makeFunction(name:"dual_scan")!)
    let threadCount=min(4*pipeline.threadExecutionWidth,pipeline.maxTotalThreadsPerThreadgroup)
    var settings=Settings(offsets:UInt32(n[4]),halfWidths:UInt32(n[5]),panels:UInt32(n[6]))
    var grid=[Core]()
    for a in 0..<n[0] {
        for w in 0..<n[1] {
            for p in 0..<n[2] {
                for h in 0..<n[3] {
                    grid.append(Core(center:Float(logGrid(16,512,a,n[0])),width:Float(logGrid(0.125,8,w,n[1])),
                        compression:Float(p)/Float(n[2]-1),drift:-8+16*Float(h)/Float(n[3]-1)))
                }
            }
        }
    }

    print("device=\(device.name) cores=\(grid.count) offsets=\(n[4]) halfwidths=\(n[5]) panels=\(n[6])")
    fflush(stdout)
    let points=n[4]*n[5]
    var bestByRegime=[String:Hit]()
    var minimumByRegime=[String:Hit]()
    var counts=[String:Int]()
    var finite=0
    var above=0
    var gpuSeconds=0.0
    for start in stride(from:0,to:grid.count,by:n[7]) {
        let end=min(start+n[7],grid.count)
        let cores=Array(grid[start..<end])
        let count=cores.count*points
        let input=device.makeBuffer(bytes:cores,length:cores.count*MemoryLayout<Core>.stride,options:.storageModeShared)!
        let output=device.makeBuffer(length:count*MemoryLayout<Float>.stride,options:.storageModeShared)!
        let command=queue.makeCommandBuffer()!
        let encoder=command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(input,offset:0,index:0)
        encoder.setBytes(&settings,length:MemoryLayout<Settings>.stride,index:1)
        encoder.setBuffer(output,offset:0,index:2)
        encoder.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threadCount,height:1,depth:1))
        encoder.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {throw command.error ?? NSError(domain:"Metal failed",code:4)}
        let elapsed=command.gpuEndTime-command.gpuStartTime
        gpuSeconds += elapsed
        let values=output.contents().bindMemory(to:Float.self,capacity:count)
        var tileBest:Float=0
        for j in 0..<cores.count {
            let p=cores[j]
            var candidateIndex=0
            var candidate:Float=0
            for i in 0..<points {
                let v=values[j*points+i]
                if !v.isFinite {continue}
                finite += 1
                if v>1 {above+=1}
                if v>candidate {candidate=v; candidateIndex=i}
            }
            tileBest=max(tileBest,candidate)
            let offset = -2+4*Float(candidateIndex/n[5])/Float(n[4]-1)
            let halfWidth=Float(logGrid(0.0625,2,candidateIndex%n[5],n[5]))
            let best=Hit(core:p,halfWidth:halfWidth,offset:offset,ratio:candidate)
            var keys=["overall"]
            if p.width<=0.5 {keys.append("widthAtMostHalf")}
            if p.width<=1.0 {keys.append("widthAtMostOne")}
            if p.width<=2.0 {keys.append("widthAtMostTwo")}
            if p.center>=64 {keys.append("centerAtLeast64")}
            if p.compression==0 {keys.append("completeCollapse")}
            if p.compression==1 {keys.append("noCompression")}
            if p.width<=1 && p.compression==0 {keys.append("completeCollapseWidthAtMostOne")}
            if p.width<=1 && p.compression==1 {keys.append("noCompressionWidthAtMostOne")}
            for key in keys {
                if candidate>(bestByRegime[key]?.ratio ?? 0) {bestByRegime[key]=best}
                if candidate>1 {
                    counts[key,default:0]+=1
                    let old=minimumByRegime[key]
                    if old == nil || p.drift<old!.core.drift || (p.drift==old!.core.drift && candidate>old!.ratio) {
                        minimumByRegime[key]=best
                    }
                }
            }
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(elapsed) best=\(tileBest)")
        fflush(stdout)
    }
    let e=JSONEncoder()
    let payload:[String:Any]=["device":device.name,"profiles":grid.count,"dimensions":Array(n[0...6]),"evaluated":grid.count*points,
        "finite":finite,"aboveUnity":above,"gpuSeconds":gpuSeconds,
        "minimumDriftByRegime":try JSONSerialization.jsonObject(with:e.encode(minimumByRegime)),
        "bestByRegime":try JSONSerialization.jsonObject(with:e.encode(bestByRegime)),"witnessedCoresByRegime":counts]
    let data=try JSONSerialization.data(withJSONObject:payload,options:[.prettyPrinted,.sortedKeys])
    try data.write(to:URL(fileURLWithPath:args[9]),options:.atomic)
    print("finite=\(finite) above=\(above) gpu_sum_seconds=\(gpuSeconds)")
}
try run()
