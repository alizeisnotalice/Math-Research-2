import Foundation
import Metal
import Darwin
struct Core:Codable {var b,t,beta,Ccrit:Float}
struct Settings {var panels,lambdas,offsets:UInt32}
struct Hit:Codable {var profile:Core;var lambda,offset:Double;var result:[Float]}
func run()throws{
    let a=CommandLine.arguments
    guard a.count==9 else{throw NSError(domain:"b t beta offsets lambda panels batch out",code:1)}
    let n=a[1...7].compactMap(Int.init)
    guard n.count==7,n[0...4].allSatisfy({$0>=2}),n[5]>0,n[6]>0 else{throw NSError(domain:"dimensions",code:2)}
    var grid=[Core]()
    for i in 0..<n[0]{
        let b=0.25*pow(16,Double(i)/Double(n[0]-1))
        for j in 0..<n[1]{
            let t=0.25*pow(16,Double(j)/Double(n[1]-1))
            let logCap=log(2*b/t)-0.5*log(2*Double.pi*t)-b*b/(2*t)
            for k in 0..<n[2]{
                grid.append(Core(b:Float(b),t:Float(t),beta:3.25+1.5*Float(k)/Float(n[2]-1),Ccrit:Float(b*b-2*logCap)))
            }
        }
    }
    guard let d=MTLCreateSystemDefaultDevice(),let q=d.makeCommandQueue() else{throw NSError(domain:"Metal device",code:3)}
    let opts=MTLCompileOptions();opts.mathMode = .safe
    let source=try String(contentsOfFile:"laplace.metal",encoding:.utf8)
    let lib=try d.makeLibrary(source:source,options:opts)
    let kernel=try d.makeComputePipelineState(function:lib.makeFunction(name:"scan_laplace")!)
    let threads=min(kernel.maxTotalThreadsPerThreadgroup,4*kernel.threadExecutionWidth)
    var cfg=Settings(panels:UInt32(n[5]),lambdas:UInt32(n[4]),offsets:UInt32(n[3]))
    let points=n[3]*n[4]
    var gpu=0.0,finite=0
    var counts=[String:Int](),positive=[String:Int](),hits=[String:Hit]()
    print("device=\(d.name) profiles=\(grid.count) points_per_profile=\(points) panels=\(n[5])");fflush(stdout)
    for lo in stride(from:0,to:grid.count,by:n[6]){
        let hi=min(grid.count,lo+n[6]),tile=Array(grid[lo..<hi]),count=tile.count*points
        let input=d.makeBuffer(bytes:tile,length:tile.count*MemoryLayout<Core>.stride,options:.storageModeShared)!
        let output=d.makeBuffer(length:count*MemoryLayout<SIMD4<Float>>.stride,options:.storageModeShared)!
        let cb=q.makeCommandBuffer()!,e=cb.makeComputeCommandEncoder()!
        e.setComputePipelineState(kernel);e.setBuffer(input,offset:0,index:0)
        e.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:1);e.setBuffer(output,offset:0,index:2)
        e.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:threads,height:1,depth:1))
        e.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed else{throw cb.error ?? NSError(domain:"GPU failed",code:4)}
        let time=cb.gpuEndTime-cb.gpuStartTime;gpu+=time
        let v=output.contents().bindMemory(to:SIMD4<Float>.self,capacity:count)
        for j in 0..<tile.count{
            let p=tile[j]
            let tag=p.beta==4 ? "critical": p.beta<4 ? "subcritical":"supercritical"
            var localCounts=0,localPositive=0
            for i in 0..<points{
                let z=v[j*points+i]
                if !z.x.isFinite {continue}
                finite+=1;localCounts+=1
                if z.x>0 {localPositive+=1}
                let selectCritical=p.beta==4 && i%n[3]==(n[3]-1)/2
                let high=i/n[3]==n[4]-1
                var keys:[String]=[]
                if z.x>(hits[tag]?.result[0] ?? -.infinity){keys.append(tag)}
                if selectCritical && z.x>(hits["criticalEqualLeading"]?.result[0] ?? -.infinity){keys.append("criticalEqualLeading")}
                if high && selectCritical && abs(z.w)>(hits["maxFirstCorrectionError"]?.result[3].magnitude ?? 0){
                    keys.append("maxFirstCorrectionError")
                }
                if !keys.isEmpty{
                    let lam=32*pow(2048,Double(i/n[3])/Double(n[4]-1))
                    let off=2*Double(i%n[3])/Double(n[3]-1)-1
                    let hit=Hit(profile:p,lambda:lam,offset:off,result:[z.x,z.y,z.z,z.w])
                    for k in keys {hits[k]=hit}
                }
            }
            counts[tag,default:0]+=localCounts;positive[tag,default:0]+=localPositive
        }
        print("batch=\(lo)..<\(hi) gpu_seconds=\(time)");fflush(stdout)
    }
    let o:[String:Any]=["device":d.name,"evaluated":grid.count*points,"finite":finite,"gpuSeconds":gpu,
        "dimensions":Array(n[0...5]),"counts":counts,"positiveCounts":positive,
        "hits":try JSONSerialization.jsonObject(with:JSONEncoder().encode(hits))]
    try JSONSerialization.data(withJSONObject:o,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[8]),options:.atomic)
    print("finite=\(finite) gpu_sum_seconds=\(gpu)")
}
try run()
