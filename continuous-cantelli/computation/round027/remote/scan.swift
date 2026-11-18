import Foundation
import Metal
import Darwin

struct Profile: Codable {
    var lambda, h, m, downWidth, plateauWidth, logBase: Float
}
struct Settings { var offsets, panels: UInt32 }
struct Crossing: Codable {
    var profile: Profile
    var curvature: Double
    var lowerKappa, upperKappa, lowerLogF, upperLogF: Float
}

func logGrid(_ lo: Double, _ hi: Double, _ i: Int, _ n: Int) -> Double {
    exp(log(lo)+(log(hi)-log(lo))*Double(i)/Double(n-1))
}
func smoothEdge(_ s: Double) -> Double {
    if s <= 0 { return 0 }
    if s >= 1 { return 1 }
    let e = exp(-1/s)
    return e/(e+exp(-1/(1-s)))
}
func baseline(_ l: Double, _ h: Double, _ c: Double) -> Double {
    let normal = sqrt(2 * Double.pi)
    let radius = 1/sqrt(c)
    let innerRadius = radius/2
    let precision = 1+l*l*h*c
    let gaussianMean = l/precision
    let root = sqrt(precision/2)
    let mass = (erf((innerRadius-gaussianMean)*root)-erf((-innerRadius-gaussianMean)*root))/2
    let core = mass*exp(l*l/(2*precision))/sqrt(precision)
    let restore = exp(-(h-1)*l*l/2)
    let backgroundInside = restore*(erf((radius-l)/sqrt(2))-erf((-radius-l)/sqrt(2)))/2
    let n = 256
    let width = innerRadius/Double(n)
    var shell = 0.0
    for i in 0...n {
        let x = innerRadius + Double(i)*width
        let transition = smoothEdge(2*sqrt(c)*x-1)
        let cutoff = c*x*x*(1-transition)+transition
        let e = -x*x/2 - l*l*h*cutoff/2
        let weight = i == 0 || i == n ? 1.0 : (i.isMultiple(of: 2) ? 2.0 : 4.0)
        shell += weight*(exp(e+l*x)+exp(e-l*x))
    }
    return log(core + shell*width/(3*normal) + restore - backgroundInside)
}

func run() throws {
    let args = CommandLine.arguments
    guard args.count == 11 else { throw NSError(domain: "l m h c ramps plateaus offsets panels batch out", code: 1) }
    let n = args[1...9].compactMap(Int.init)
    guard n.count == 9, n[0...6].allSatisfy({$0>=2}), n[7]>0, n[8]>0 else { throw NSError(domain: "Dimensions", code: 2) }
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else { throw NSError(domain:"Metal device",code:3) }
    let source=try String(contentsOfFile:"tail.metal",encoding:.utf8)
    let compile=MTLCompileOptions()
    compile.mathMode = .safe
    let library=try device.makeLibrary(source:source,options:compile)
    let pipeline=try device.makeComputePipelineState(function:library.makeFunction(name:"scan_tail")!)
    let threadCount=min(4*pipeline.threadExecutionWidth,pipeline.maxTotalThreadsPerThreadgroup)
    var settings=Settings(offsets:UInt32(n[6]),panels:UInt32(n[7]))
    var profiles=[Profile]()
    var curvatures=[Double]()
    var skipped=0
    for li in 0..<n[0] {
        let l=logGrid(16,256,li,n[0])
        for hi in 0..<n[2] {
            let h=logGrid(1.25,4,hi,n[2])
            for ci in 0..<n[3] {
                let c=logGrid(0.25,4,ci,n[3])
                let b=baseline(l,h,c)
                for mi in 0..<n[1] {
                    let m=logGrid(0.1,0.96,mi,n[1])
                    for di in 0..<n[4] {
                        let d=logGrid(0.01,96,di,n[4])
                        // Keep the entire support right of x=lambda. This
                        // guarantees disjointness and monotonic translation.
                        if d+4*log(l)>=l*l*(1-m) {skipped += n[5]; continue}
                        for wi in 0..<n[5] {
                            let w=logGrid(0.01,48,wi,n[5])
                            profiles.append(Profile(lambda:Float(l),h:Float(h),m:Float(m),downWidth:Float(d),plateauWidth:Float(w),logBase:Float(b)))
                            curvatures.append(c)
                        }
                    }
                }
            }
        }
    }
    print("device=\(device.name) cores=\(profiles.count) offsets=\(n[6]) panels=\(n[7])")
    fflush(stdout)
    var bestByRegime=[String:Crossing]()
    var seconds=0.0
    var finite=0
    var crossed=0
    var positiveAtLast=0
    var negativeAtFirst=0
    var monotoneViolations=0
    var maxAbsScaledBaseRelative=0.0
    for start in stride(from:0,to:profiles.count,by:n[8]) {
        let end=min(start+n[8],profiles.count)
        let tile=Array(profiles[start..<end])
        let count=tile.count*n[6]
        let input=device.makeBuffer(bytes:tile,length:tile.count*MemoryLayout<Profile>.stride,options:.storageModeShared)!
        let output=device.makeBuffer(length:count*MemoryLayout<SIMD2<Float>>.stride,options:.storageModeShared)!
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
        guard command.status == .completed else { throw command.error ?? NSError(domain:"Metal failed",code:4) }
        let elapsed=command.gpuEndTime-command.gpuStartTime
        seconds += elapsed
        let values=output.contents().bindMemory(to:SIMD2<Float>.self,capacity:count)
        for j in 0..<tile.count {
            let p=profiles[start+j]
            let c=curvatures[start+j]
            let asymptotic=exp(1/(2*c*Double(p.h)))/sqrt(c*Double(p.h))
            maxAbsScaledBaseRelative=max(maxAbsScaledBaseRelative,abs(Double(p.lambda)*exp(Double(p.logBase))/asymptotic-1))
            let base=j*n[6]
            if values[base].x<0 { negativeAtFirst += 1 }
            if values[base+n[6]-1].x>0 { positiveAtLast += 1 }
            var crossingIndex: Int?=nil
            for i in 0..<n[6] {
                let v=values[base+i].x
                finite += v.isFinite ? 1 : 0
                if i>0 {
                    let old=values[base+i-1].x
                    if v>old+1e-6 { monotoneViolations += 1 }
                    if old>=0 && v<0 { crossingIndex=i }
                }
            }
            guard let i=crossingIndex else {continue}
            crossed += 1
            let step: Float=16/Float(n[6]-1)
            let hit=Crossing(profile:p,curvature:c,lowerKappa:-4+step*Float(i-1),upperKappa:-4+step*Float(i),lowerLogF:values[base+i-1].x,upperLogF:values[base+i].x)
            let keyList=["largestTransition"]+(p.downWidth<=0.05 ? ["sharpRamp"] : [])+(p.lambda>=128 ? ["lambdaAtLeast128"] : [])+(p.plateauWidth<=0.1 ? ["shortPlateau"] : [])
            for key in keyList where hit.lowerKappa>(bestByRegime[key]?.lowerKappa ?? -.infinity) {bestByRegime[key]=hit}
        }
        print("batch=\(start)..<\(end) gpu_seconds=\(elapsed) crossing_profiles=\(crossed)")
        fflush(stdout)
    }
    let e=JSONEncoder()
    let o:[String:Any]=["device":device.name,"profileCount":profiles.count,"finite":finite,"evaluated":profiles.count*n[6],"crossed":crossed,
        "skippedProfiles":skipped,"positiveAtLast":positiveAtLast,"negativeAtFirst":negativeAtFirst,"monotoneViolations":monotoneViolations,
        "gpuSeconds":seconds,"dimensions":Array(n[0...7]),"maxOriginRelativeDiscrepancy":maxAbsScaledBaseRelative,
        "bestByRegime":try JSONSerialization.jsonObject(with:e.encode(bestByRegime))]
    let data=try JSONSerialization.data(withJSONObject:o,options:[.prettyPrinted,.sortedKeys])
    try data.write(to:URL(fileURLWithPath:args[10]),options:.atomic)
    print("valid=\(finite) gpu_sum_seconds=\(seconds) crosses=\(crossed) positive_last=\(positiveAtLast) negative_first=\(negativeAtFirst)")
}
try run()
