import Darwin
import Foundation
import Metal

struct Level {var length,shift,tail,pad:Float}
struct Profile {var variance,c:Float;var geometry,pad:UInt32}
struct Settings {var points,begin,uniformPoints,pad:UInt32}
struct Result:Codable {
    var c,v,center,maxValue,x,centerCurvature,maxSlope:Float
}
func run() throws {
    let args=CommandLine.arguments
    guard args.count==6,let n=Int(args[1]),let nV=Int(args[2]),let nX=Int(args[3]),let tile=Int(args[4]),
        n>=2,nV>=2,nX>=64,tile>0,let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue()
    else {throw NSError(domain:"c_variations v_variations x_points tile output",code:1)}
    var cList=(0..<n).map {0.001953125*pow(704.0,Float($0)/Float(n-1))}
    cList += [1.0/256, 1.0/64, 0.125, 0.25, 0.5, 0.75, 1.0, 1.25]
    cList=Array(Set(cList)).sorted()
    var geometries=[Level](),profiles=[Profile]()
    let nxUniform=min(257,nX/2)
    for (i,c) in cList.enumerated(){
        var lengths=[1.0]
        for j in 1...52 {lengths.append(lengths[j-1]*(1-Double(c)*pow(2.0,-Double(j)))*0.5)}
        var p=[Double](repeating:1,count:53)
        for j in stride(from:51,through:0,by:-1) {p[j]=p[j+1]*(1-Double(c)*pow(2.0,-Double(j+1)))}
        for j in 0...21 {
            geometries.append(Level(length:Float(lengths[j]),shift:j==0 ? 0:Float((lengths[j-1]-lengths[j])*0.5),tail:Float(p[j]),pad:0))
        }
        for j in 0..<nV{
            let v=pow(2.0,-32.0+31.0*Float(j)/Float(nV-1))
            profiles.append(Profile(variance:v,c:c,geometry:UInt32(i),pad:0))
        }
    }
    let options=MTLCompileOptions();options.mathMode = .safe
    let source=try String(contentsOfFile:"cantor.metal",encoding:.utf8)
    let lib=try device.makeLibrary(source:source,options:options)
    let pipe=try device.makeComputePipelineState(function:lib.makeFunction(name:"cantor_scan")!)
    let dGeom=device.makeBuffer(bytes:geometries,length:geometries.count*MemoryLayout<Level>.stride,options:.storageModeShared)!
    let dProfiles=device.makeBuffer(bytes:profiles,length:profiles.count*MemoryLayout<Profile>.stride,options:.storageModeShared)!
    let buf=device.makeBuffer(length:tile*nX*MemoryLayout<SIMD4<Float>>.stride,options:.storageModeShared)!
    let size=min(pipe.maxTotalThreadsPerThreadgroup,4*pipe.threadExecutionWidth)
    let utc=ISO8601DateFormatter();utc.formatOptions=[.withInternetDateTime,.withFractionalSeconds]
    let wall=ProcessInfo.processInfo.systemUptime,start=utc.string(from:Date())
    var rows=[Result](),gpu=0.0,cmd=0
    print("start \(start) pid=\(getpid()) device=\(device.name) c=\(cList.count) v=\(nV) x=\(nX)")
    fflush(stdout)
    for off in stride(from:0,to:profiles.count,by:tile){
        let count=min(tile,profiles.count-off)
        let cb=queue.makeCommandBuffer()!,enc=cb.makeComputeCommandEncoder()!
        var cfg=Settings(points:UInt32(nX),begin:UInt32(off),uniformPoints:UInt32(nxUniform),pad:0)
        enc.setComputePipelineState(pipe)
        enc.setBuffer(dGeom,offset:0,index:0);enc.setBuffer(dProfiles,offset:0,index:1)
        enc.setBytes(&cfg,length:MemoryLayout<Settings>.stride,index:2);enc.setBuffer(buf,offset:0,index:3)
        enc.dispatchThreads(MTLSize(width:count*nX,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:size,height:1,depth:1))
        enc.endEncoding();cb.commit();cb.waitUntilCompleted()
        guard cb.status == .completed,cb.gpuEndTime>cb.gpuStartTime else {
            throw cb.error ?? NSError(domain:"Metal command or clock",code:2)
        }
        gpu+=cb.gpuEndTime-cb.gpuStartTime;cmd+=1
        let ptr=buf.contents().bindMemory(to:SIMD4<Float>.self,capacity:count*nX)
        for i in 0..<count{
            let test=profiles[off+i],center=ptr[i*nX]
            var winner=center
            for k in 1..<nX{
                let candidate=ptr[i*nX+k]
                if candidate.x>winner.x {winner=candidate}
            }
            rows.append(Result(c:test.c,v:test.variance,center:center.x,maxValue:winner.x,x:winner.w,
                centerCurvature:center.z,maxSlope:winner.y))
        }
        if cmd%10==0 || cmd<3 || off+count==profiles.count{
            print("sample \(utc.string(from:Date())) commands=\(cmd) cases=\(off+count) gpu=\(gpu)")
            fflush(stdout)
        }
    }
    let summary:[String:Any]=["device":device.name,"cProfiles":cList.count,"variances":nV,"spatialSites":nX,"tile":tile,
       "sites":profiles.count*nX,"commands":cmd,"gpuSeconds":gpu,"wallSeconds":ProcessInfo.processInfo.systemUptime-wall,
       "startUTC":start,"endUTC":utc.string(from:Date()),"rows":try JSONSerialization.jsonObject(with:JSONEncoder().encode(rows))]
    try JSONSerialization.data(withJSONObject:summary,options:[.sortedKeys])
        .write(to:URL(fileURLWithPath:args[5]),options:.atomic)
}
try run()
