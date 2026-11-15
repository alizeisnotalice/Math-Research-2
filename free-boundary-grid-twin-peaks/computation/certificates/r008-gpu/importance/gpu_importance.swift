import Foundation
import Metal

// Importance sampler for legal connected completions on an n by n free grid.
// Proposal: choose uniformly from the current outer vertex boundary.  A path
// gets weight product_t |d S_t|, so E_q[weight]=h(S).  This is exploratory
// floating-point evidence only; the ratio estimator is biased at finite N.

let metalSource = """
#include <metal_stdlib>
using namespace metal;
inline uint nextRand(thread uint &x) { x ^= x << 13; x ^= x >> 17; x ^= x << 5; return x; }
struct Params { uint n,N,count; uint s0,s1,s2,s3,s4; uint c0,c1,c2,c3,c4; };
struct Result { float logW; uint event,w,v,b0,b1; };
inline bool has(thread uint *a,uint v) { return (a[v>>5] & (1u<<(v&31))) != 0; }
inline void put(thread uint *a,uint v) { a[v>>5] |= 1u<<(v&31); }
inline bool touches(thread uint *a,uint v,uint n,uint N) {
  uint x=v%n, y=v/n;
  return (x && has(a,v-1)) || (x+1<n && has(a,v+1)) ||
         (y && has(a,v-n)) || (y+1<n && has(a,v+n));
}
kernel void sample(constant Params &p [[buffer(0)]], device uint *seeds [[buffer(1)]],
                   device Result *out [[buffer(2)]], uint tid [[thread_position_in_grid]]) {
  if(tid>=p.count) return;
  thread uint marked[5]={p.s0,p.s1,p.s2,p.s3,p.s4};
  thread uint boundary[144]; uint state=seeds[tid]; float lw=0.0f;
  uint ev=0u, first=0xffffffffu, second=0xffffffffu, bfirst=0u, bsecond=0u;
  for(uint step=0; step<p.N; ++step) {
    uint nb=0u;
    for(uint z=0; z<p.N; ++z) if(!has(marked,z) && touches(marked,z,p.n,p.N)) boundary[nb++]=z;
    if(nb==0u) break;
    uint chosen=boundary[nextRand(state)%nb];
    if(step==0u) { first=chosen; bfirst=nb; }
    if(step==1u) { second=chosen; bsecond=nb; uint c=(chosen>>5); if(((c==0u?p.c0:c==1u?p.c1:c==2u?p.c2:c==3u?p.c3:p.c4) & (1u<<(chosen&31)))==0u) ev=1u; }
    lw += log(float(nb)); put(marked,chosen);
  }
  out[tid]=(Result){lw,ev,first,second,bfirst,bsecond}; seeds[tid]=state;
}
"""

struct ParamsSwift { var n,N,count: UInt32; var s0,s1,s2,s3,s4: UInt32; var c0,c1,c2,c3,c4: UInt32 }
struct GPUResult { var logW: Float; var event,w,v,b0,b1: UInt32 }
struct Shape { let name:String; let n:Int; let cells:[(Int,Int)]; let exactScaled:Double? }
struct Trace { let logW:Double; let event:Bool; let first:Int; let second:Int; let b0:Int; let b1:Int; let allBounds:[Int] }

func nowISO() -> String { ISO8601DateFormatter().string(from: Date()) }
func words(_ n:Int, _ cells:[(Int,Int)]) -> [UInt32] {
  var a=[UInt32](repeating:0,count:5)
  for (x,y) in cells { let v=y*n+x; a[v>>5] |= UInt32(1)<<UInt32(v&31) }
  return a
}
func closed(_ n:Int,_ seed:[UInt32]) -> [UInt32] {
  var r=seed; let N=n*n
  for v in 0..<N where (seed[v>>5] & (UInt32(1)<<UInt32(v&31))) != 0 {
    let x=v%n,y=v/n
    for q in [x>0 ? v-1 : -1, x+1<n ? v+1 : -1, y>0 ? v-n : -1, y+1<n ? v+n : -1] where q >= 0 { r[q>>5] |= UInt32(1)<<UInt32(q&31) }
  }
  return r
}
func marked(_ a:[UInt32],_ v:Int)->Bool { (a[v>>5] & (UInt32(1)<<UInt32(v&31))) != 0 }
func put(_ a: inout [UInt32],_ v:Int) { a[v>>5] |= UInt32(1)<<UInt32(v&31) }
func rng(_ x: inout UInt32)->UInt32 { x ^= x &<< 13; x ^= x >> 17; x ^= x &<< 5; return x }
func cpuTrace(n:Int, start:[UInt32], c:[UInt32], seed:UInt32)->Trace {
  let N=n*n; var a=start, s=seed, lw=0.0, bs=[Int](), first = -1, second = -1, ev=false
  for step in 0..<N {
    var b=[Int](); b.reserveCapacity(N)
    for v in 0..<N where !marked(a,v) { let x=v%n,y=v/n
      if (x>0 && marked(a,v-1)) || (x+1<n && marked(a,v+1)) || (y>0 && marked(a,v-n)) || (y+1<n && marked(a,v+n)) { b.append(v) }
    }
    if b.isEmpty { break }; bs.append(b.count); let z=b[Int(rng(&s)%UInt32(b.count))]
    if step==0 { first=z }; if step==1 { second=z; ev = !marked(c,z) }
    lw += log(Double(b.count)); put(&a,z)
  }
  return Trace(logW:lw,event:ev,first:first,second:second,b0:bs.first ?? 0,b1:bs.dropFirst().first ?? 0,allBounds:bs)
}
func logAdd(_ x:Double,_ y:Double)->Double { if x == -.infinity { return y }; if y == -.infinity { return x }; let m=max(x,y); return m + log(exp(x-m)+exp(y-m)) }
func cellText(_ n:Int,_ v:Int)->String { v < 0 ? "null" : "[\(v%n),\(v/n)]" }
func jsonEscape(_ s:String)->String { s.replacingOccurrences(of:"\\",with:"\\\\").replacingOccurrences(of:"\"",with:"\\\"") }

let outDir=URL(fileURLWithPath:CommandLine.arguments.dropFirst().first ?? ".")
try FileManager.default.createDirectory(at:outDir, withIntermediateDirectories:true)
let output=outDir.appendingPathComponent("results.jsonl")
FileManager.default.createFile(atPath:output.path, contents:nil)
let fh=try FileHandle(forWritingTo:output)
func write(_ s:String) { try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!) }

let dev=MTLCreateSystemDefaultDevice()!
let library=try dev.makeLibrary(source:metalSource, options:nil)
let pipeline=try dev.makeComputePipelineState(function:library.makeFunction(name:"sample")!)
let queue=dev.makeCommandQueue()!
let samplesPerBatch=16384
let batchPlan:[Int:Int]=[4:30,5:60,6:40,8:100,12:100]
// The streams are deliberately not adjacent ranges: base+i would make the
// three nominal seeds share almost every initial PRNG state in a batch.
let seeds:[UInt32]=[0xA0220801,0xB14E1D53,0xC35A29A7]
let shapes:[Shape]=[
  Shape(name:"n4_singleton_exactmax",n:4,cells:[(1,0)],exactScaled:Double(276996)/Double(137767)),
  Shape(name:"n5_four_point_exactmax",n:5,cells:[(1,0),(2,0),(3,0),(1,1)],exactScaled:Double(6716680961734.0)/Double(2581972569155.0)),
  Shape(name:"n6_corner_cluster",n:6,cells:[(0,0),(1,0),(0,1),(1,1)],exactScaled:nil),
  Shape(name:"n6_four_point_peak",n:6,cells:[(1,0),(2,0),(3,0),(1,1)],exactScaled:nil),
  Shape(name:"n6_concave_interface",n:6,cells:[(1,1),(2,1),(3,1),(1,2),(3,2),(1,3),(2,3),(3,3)],exactScaled:nil),
  Shape(name:"n8_corner_cluster",n:8,cells:[(0,0),(1,0),(0,1),(1,1),(2,0)],exactScaled:nil),
  Shape(name:"n8_four_point_peak",n:8,cells:[(2,0),(3,0),(4,0),(2,1)],exactScaled:nil),
  Shape(name:"n8_concave_interface",n:8,cells:[(2,2),(3,2),(4,2),(2,3),(4,3),(2,4),(3,4),(4,4)],exactScaled:nil),
  Shape(name:"n12_corner_cluster",n:12,cells:[(0,0),(1,0),(0,1),(1,1),(2,0),(0,2)],exactScaled:nil),
  Shape(name:"n12_four_point_peak",n:12,cells:[(4,0),(5,0),(6,0),(4,1)],exactScaled:nil),
  Shape(name:"n12_concave_interface",n:12,cells:[(4,4),(5,4),(6,4),(7,4),(4,5),(7,5),(4,6),(7,6),(4,7),(5,7),(6,7),(7,7)],exactScaled:nil)
]
write("{\"kind\":\"metadata\",\"time\":\"\(nowISO())\",\"device\":\"\(jsonEscape(dev.name))\",\"swift\":\"6.4.2\",\"proposal\":\"uniform-current-boundary\",\"samples_per_batch\":\(samplesPerBatch),\"seeds\":[\(seeds.map(String.init).joined(separator:","))]}")

for sh in shapes {
  let start=words(sh.n,sh.cells), c=closed(sh.n,start), N=sh.n*sh.n
  var p=ParamsSwift(n:UInt32(sh.n),N:UInt32(N),count:UInt32(samplesPerBatch),s0:start[0],s1:start[1],s2:start[2],s3:start[3],s4:start[4],c0:c[0],c1:c[1],c2:c[2],c3:c[3],c4:c[4])
  let pb=dev.makeBuffer(bytes:&p,length:MemoryLayout<ParamsSwift>.stride,options:.storageModeShared)!
  let sb=dev.makeBuffer(length:samplesPerBatch*MemoryLayout<UInt32>.stride,options:.storageModeShared)!
  let rb=dev.makeBuffer(length:samplesPerBatch*MemoryLayout<GPUResult>.stride,options:.storageModeShared)!
  let tr=cpuTrace(n:sh.n,start:start,c:c,seed:seeds[0])
  // Same seed, proposal and linear boundary order: a one-thread GPU/CPU trajectory check.
  p.count=1; pb.contents().bindMemory(to:ParamsSwift.self,capacity:1).pointee=p
  sb.contents().bindMemory(to:UInt32.self,capacity:1).pointee=seeds[0]
  do { let cb=queue.makeCommandBuffer()!; let enc=cb.makeComputeCommandEncoder()!; enc.setComputePipelineState(pipeline); enc.setBuffer(pb,offset:0,index:0); enc.setBuffer(sb,offset:0,index:1); enc.setBuffer(rb,offset:0,index:2); enc.dispatchThreads(MTLSize(width:1,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:1,height:1,depth:1)); enc.endEncoding(); cb.commit(); cb.waitUntilCompleted() }
  let gr=rb.contents().bindMemory(to:GPUResult.self,capacity:1).pointee
  let traceOK=abs(Double(gr.logW)-tr.logW)<1e-3 && gr.event == (tr.event ? 1:0) && Int(gr.w)==tr.first && Int(gr.v)==tr.second && Int(gr.b0)==tr.b0 && Int(gr.b1)==tr.b1
  write("{\"kind\":\"trajectory_check\",\"shape\":\"\(sh.name)\",\"ok\":\(traceOK),\"cpu\":{\"logw\":\(tr.logW),\"event\":\(tr.event),\"w\":\(cellText(sh.n,tr.first)),\"v\":\(cellText(sh.n,tr.second)),\"b0\":\(tr.b0),\"b1\":\(tr.b1),\"boundaries\":[\(tr.allBounds.map(String.init).joined(separator:","))]},\"gpu\":{\"logw\":\(gr.logW),\"event\":\(gr.event),\"w\":\(cellText(sh.n,Int(gr.w))),\"v\":\(cellText(sh.n,Int(gr.v))),\"b0\":\(gr.b0),\"b1\":\(gr.b1)}}")
  precondition(traceOK,"GPU/CPU trajectory mismatch for \(sh.name)")
  p.count=UInt32(samplesPerBatch); pb.contents().bindMemory(to:ParamsSwift.self,capacity:1).pointee=p
  for (_,base) in seeds.enumerated() {
    var logAll:Double = -.infinity, logEvt:Double = -.infinity, logSq:Double = -.infinity, maxW:Double = -.infinity
    var rawEvent=0, total=0, gpuSeconds=0.0, wallSeconds=0.0, hwAvailable=true
    for batch in 0..<(batchPlan[sh.n] ?? 8) {
      let sp=sb.contents().bindMemory(to:UInt32.self,capacity:samplesPerBatch)
      for i in 0..<samplesPerBatch { sp[i] = base &+ (UInt32(i) &* 0x9e3779b9) &+ (UInt32(batch) &* 0x7f4a7c15) }
      let cb=queue.makeCommandBuffer()!; let enc=cb.makeComputeCommandEncoder()!
      enc.setComputePipelineState(pipeline); enc.setBuffer(pb,offset:0,index:0); enc.setBuffer(sb,offset:0,index:1); enc.setBuffer(rb,offset:0,index:2)
      let tpg=min(pipeline.maxTotalThreadsPerThreadgroup,256)
      enc.dispatchThreads(MTLSize(width:samplesPerBatch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:tpg,height:1,depth:1)); enc.endEncoding()
      let began=Date(); cb.commit(); cb.waitUntilCompleted(); let wall=Date().timeIntervalSince(began), hw=cb.gpuEndTime-cb.gpuStartTime
      wallSeconds += wall; if hw > 0 { gpuSeconds += hw } else { hwAvailable=false }
      let rp=rb.contents().bindMemory(to:GPUResult.self,capacity:samplesPerBatch)
      for i in 0..<samplesPerBatch { let z=Double(rp[i].logW); logAll=logAdd(logAll,z); logSq=logAdd(logSq,2*z); if rp[i].event != 0 { logEvt=logAdd(logEvt,z); rawEvent += 1 }; maxW=max(maxW,z); total += 1 }
    }
    let pHat=exp(logEvt-logAll), ess=exp(2*logAll-logSq), maxShare=exp(maxW-logAll), B=tr.b0 // initial boundary is deterministic
    let scaled=Double(B)*pHat, rate=Double(total)/max(wallSeconds,1e-12)
    let exact = sh.exactScaled.map { ",\"exact_scaled\":\($0),\"error_scaled\":\(scaled-$0)" } ?? ""
    let gpuField = hwAvailable ? String(gpuSeconds) : "null"
    write("{\"kind\":\"seed_result\",\"shape\":\"\(sh.name)\",\"n\":\(sh.n),\"seed\":\(base),\"samples\":\(total),\"raw_event_samples\":\(rawEvent),\"boundary\":\(B),\"p_hat\":\(pHat),\"scaled_boundary_p\":\(scaled),\"ess\":\(ess),\"ess_fraction\":\(ess/Double(total)),\"max_weight_share\":\(maxShare),\"metal_gpu_seconds\":\(gpuField),\"submit_wait_wall_seconds\":\(wallSeconds),\"samples_per_wall_second\":\(rate)\(exact)}")
  }
}
try fh.close()
print("wrote \(output.path)")
