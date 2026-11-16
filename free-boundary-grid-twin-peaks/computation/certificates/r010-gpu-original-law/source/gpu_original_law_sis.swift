import Foundation
import Metal

// One GPU thread = one full descending label order.  mode 0 is the exact
// <=2-birth proposal; mode 1 is the independent one-peak proposal.
let msl = """
#include <metal_stdlib>
using namespace metal;
struct P { uint n,N,count,mode,seed; };
struct O { float logw; uint births,p1,p2,e2,rho2; };
inline uint mix(uint x) { x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16; return x | 1u; }
inline uint rnd(thread uint& x) { x ^= x << 13; x ^= x >> 17; x ^= x << 5; return x; }
inline bool has(const thread uint* a,uint v) { return (a[v>>5] & (1u << (v&31))) != 0; }
inline void add(thread uint* a,uint v) { a[v>>5] |= (1u << (v&31)); }
inline bool touch(const thread uint* a,uint v,uint n) { uint x=v%n,y=v/n; return (x && has(a,v-1)) || (x+1<n && has(a,v+1)) || (y && has(a,v-n)) || (y+1<n && has(a,v+n)); }
kernel void sample(constant P& p [[buffer(0)]], device O* out [[buffer(1)]], uint id [[thread_position_in_grid]]) {
  if (id >= p.count) return;
  thread uint a[5] = {0,0,0,0,0};
  thread uint choices[144]; uint rng=mix(p.seed ^ (id*0x9e3779b9u));
  uint births=0,p1=0xffffffffu,p2=0xffffffffu,e2=0,rho2=0; float lw=0.0f;
  for (uint rank=0; rank<p.N; ++rank) {
    uint c=0;
    for (uint v=0; v<p.N; ++v) if (!has(a,v)) {
      bool legal = (p.mode==0) ? (births<2 || touch(a,v,p.n)) : (rank==0 || touch(a,v,p.n));
      if (legal) choices[c++]=v;
    }
    if (!c) { births=0xffffffffu; break; }
    uint v=choices[rnd(rng)%c]; lw += log(float(c)); bool born=!touch(a,v,p.n);
    if (born) { if (births==0) p1=v; else if (births==1) p2=v; births++; }
    add(a,v);
  }
  if (births==2) { int dx=int(p1%p.n)-int(p2%p.n), dy=int(p1/p.n)-int(p2/p.n); uint ax=uint(abs(dx)), ay=uint(abs(dy)); e2=ax*ax+ay*ay; rho2=(ax+ay)*(ax+ay); }
  out[id]=(O){lw,births,p1,p2,e2,rho2};
}
"""

struct P { var n,N,count,mode,seed: UInt32 }
struct O { var logw: Float; var births,p1,p2,e2,rho2: UInt32 }
struct Target { let n: Int; let a1: String?; let a2: String?; let q: Double? }

func mix(_ x: UInt32) -> UInt32 { var y=x; y ^= y >> 16; y &*= 0x7feb352d; y ^= y >> 15; y &*= 0x846ca68b; y ^= y >> 16; return y | 1 }
func rnd(_ x: inout UInt32) -> UInt32 { x ^= x &<< 13; x ^= x >> 17; x ^= x &<< 5; return x }
func has(_ a: [UInt32], _ v: Int) -> Bool { (a[v >> 5] & (1 << UInt32(v & 31))) != 0 }
func touch(_ a: [UInt32], _ v: Int, _ n: Int) -> Bool { let x=v%n,y=v/n; return (x>0 && has(a,v-1)) || (x+1<n && has(a,v+1)) || (y>0 && has(a,v-n)) || (y+1<n && has(a,v+n)) }
func logadd(_ a: Double, _ b: Double) -> Double { if a == -.infinity { return b }; if b == -.infinity { return a }; let m=max(a,b); return m + log(exp(a-m)+exp(b-m)) }
func logFactorial(_ n: Int) -> Double { (2...max(2,n)).reduce(0.0) { $0 + log(Double($1)) } }
func jsonNumber(_ x: Double) -> String { x.isFinite ? String(x) : "null" }

// Exact single-thread CPU mirror used before each GPU batch family.
func trace(n: Int, mode: UInt32, base: UInt32) -> O {
  let N=n*n; var a=[UInt32](repeating:0,count:5), rng=mix(base), births:UInt32=0, p1=UInt32.max,p2=UInt32.max,e2:UInt32=0,rho2:UInt32=0; var lw=0.0
  for rank in 0..<N {
    var choices=[Int](); choices.reserveCapacity(N)
    for v in 0..<N where !has(a,v) { let legal = mode==0 ? (births<2 || touch(a,v,n)) : (rank==0 || touch(a,v,n)); if legal { choices.append(v) } }
    precondition(!choices.isEmpty)
    let v=choices[Int(rnd(&rng)%UInt32(choices.count))]; lw += log(Double(choices.count)); let born = !touch(a,v,n)
    if born { if births==0 { p1=UInt32(v) } else if births==1 { p2=UInt32(v) }; births += 1 }
    a[v>>5] |= 1 << UInt32(v&31)
  }
  if births==2 { let dx=abs(Int(p1 % UInt32(n))-Int(p2 % UInt32(n))),dy=abs(Int(p1 / UInt32(n))-Int(p2 / UInt32(n))); e2=UInt32(dx*dx+dy*dy);rho2=UInt32((dx+dy)*(dx+dy)) }
  return O(logw:Float(lw),births:births,p1:p1,p2:p2,e2:e2,rho2:rho2)
}

let arguments=Array(CommandLine.arguments.dropFirst())
let directory=URL(fileURLWithPath:arguments.first ?? ".")
let heavy=arguments.contains("--heavy")
try FileManager.default.createDirectory(at:directory,withIntermediateDirectories:true)
let path=directory.appendingPathComponent("results.jsonl")
guard !FileManager.default.fileExists(atPath:path.path) else { fatalError("output directory already has results.jsonl; choose a new directory") }
FileManager.default.createFile(atPath:path.path,contents:nil)
let fh=try FileHandle(forWritingTo:path)
func emit(_ s: String) { try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!) }
let dev=MTLCreateSystemDefaultDevice()!
let lib=try dev.makeLibrary(source:msl,options:nil)
let pipe=try dev.makeComputePipelineState(function:lib.makeFunction(name:"sample")!)
let queue=dev.makeCommandQueue()!
let batch=131072
let targets=[
 Target(n:2,a1:"16",a2:"8",q:0.5), Target(n:3,a1:"29568",a2:"143808",q:1100.0/2247.0),
 Target(n:4,a1:"35704394880",a2:"779299174080",q:784569917.0/1855474224.0),
 Target(n:5,a1:"88366736882654697600",a2:"7063261481362958280000",q:67631454448823894520000.0/(25.0*7063261481362958280000.0)),
 Target(n:6,a1:"1080582891996915763252912268789760",a2:"304351937651819767235871772704107520",q:768470066630672011431926133632875.0/2139974561614357738377223401825756.0),
 Target(n:7,a1:"134267527713152677457983522684235911142268448963584",a2:"137077928808481869759712347602605837601908319657172480",q:30737076327432172193126879326262773610108769852594.0/89243443234688717291479392970446508855409062276805.0),
 Target(n:8,a1:nil,a2:nil,q:nil), Target(n:12,a1:nil,a2:nil,q:nil)
]
emit("{\"kind\":\"metadata\",\"device\":\"\(dev.name)\",\"swift\":\"6.4.2\",\"proposal\":\"descending original-law SIS; mode0 K<=2, mode1 K=1\",\"batch\":\(batch),\"rng\":\"mix32 then xorshift32; deterministic streams only\"}")

for target in targets where !heavy || target.n >= 7 {
  let n=target.n, N=n*n, batches=(n <= 7 ? 4 : 8)
  for mode: UInt32 in [0,1] {
    var p=P(n:UInt32(n),N:UInt32(N),count:UInt32(batch),mode:mode,seed:0xa022c001)
    let pb=dev.makeBuffer(bytes:&p,length:MemoryLayout<P>.stride,options:.storageModeShared)!
    let ob=dev.makeBuffer(length:batch*MemoryLayout<O>.stride,options:.storageModeShared)!
    // CPU/GPU equality is checked for id=0, whose seed is mix(base).
    p.count=1; pb.contents().bindMemory(to:P.self,capacity:1).pointee=p
    let cpu=trace(n:n,mode:mode,base:p.seed)
    let cb0=queue.makeCommandBuffer()!, enc0=cb0.makeComputeCommandEncoder()!; enc0.setComputePipelineState(pipe); enc0.setBuffer(pb,offset:0,index:0); enc0.setBuffer(ob,offset:0,index:1); enc0.dispatchThreads(MTLSize(width:1,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:1,height:1,depth:1)); enc0.endEncoding(); cb0.commit();cb0.waitUntilCompleted()
    guard cb0.status == .completed else { fatalError("trace command failed: \(cb0.error?.localizedDescription ?? "unknown")") }
    let gpu=ob.contents().bindMemory(to:O.self,capacity:1).pointee
    precondition(abs(Double(cpu.logw)-Double(gpu.logw)) < 1e-4 && cpu.births==gpu.births && cpu.p1==gpu.p1 && cpu.p2==gpu.p2 && cpu.e2==gpu.e2 && cpu.rho2==gpu.rho2,"CPU/GPU trace mismatch")
    emit("{\"kind\":\"trajectory_check\",\"n\":\(n),\"mode\":\(mode),\"ok\":true,\"logw\":\(gpu.logw),\"births\":\(gpu.births),\"e2\":\(gpu.e2),\"rho2\":\(gpu.rho2)}")
    p.count=UInt32(batch);pb.contents().bindMemory(to:P.self,capacity:1).pointee=p
    var lw1 = -Double.infinity, lw2 = -Double.infinity, lw2e = -Double.infinity, lw2r = -Double.infinity, lw1sq = -Double.infinity, lw2sq = -Double.infinity
    var max1 = -Double.infinity, max2 = -Double.infinity, gpuSeconds=0.0, wallSeconds=0.0, k1=0,k2=0,bad=0
    for b in 0..<(heavy ? batches * 64 : batches) {
      p.seed = 0xa022c001 &+ UInt32(n*131+Int(mode)*17+b*7919); pb.contents().bindMemory(to:P.self,capacity:1).pointee=p
      let cb=queue.makeCommandBuffer()!, enc=cb.makeComputeCommandEncoder()!;enc.setComputePipelineState(pipe);enc.setBuffer(pb,offset:0,index:0);enc.setBuffer(ob,offset:0,index:1);enc.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1));enc.endEncoding();let st=Date();cb.commit();cb.waitUntilCompleted();wallSeconds += Date().timeIntervalSince(st);guard cb.status == .completed else { fatalError("batch \(b) failed: \(cb.error?.localizedDescription ?? "unknown")") };gpuSeconds += cb.gpuEndTime-cb.gpuStartTime
      let o=ob.contents().bindMemory(to:O.self,capacity:batch)
      for i in 0..<batch { let x=Double(o[i].logw); if o[i].births==1 { k1 += 1; lw1=logadd(lw1,x);lw1sq=logadd(lw1sq,2*x);max1=max(max1,x) } else if o[i].births==2 { k2 += 1;lw2=logadd(lw2,x);lw2sq=logadd(lw2sq,2*x);max2=max(max2,x);lw2e=logadd(lw2e,x+log(Double(o[i].e2)));lw2r=logadd(lw2r,x+log(Double(o[i].rho2))) } else { bad += 1 } }
      let rec="{\"kind\":\"batch\",\"n\":\(n),\"mode\":\(mode),\"batch_index\":\(b),\"gpu_seconds\":\(cb.gpuEndTime-cb.gpuStartTime),\"wall_seconds\":\(Date().timeIntervalSince(st))}";emit(rec)
    }
    let totalBatches = heavy ? batches * 64 : batches
    let M=Double(batch*totalBatches), logA1=lw1-log(M),logA2=lw2-log(M),qE=exp(lw2e-lw2)/Double(N),qR=exp(lw2r-lw2)/Double(N)
    let ess1=exp(2*lw1-lw1sq),ess2=exp(2*lw2-lw2sq), top1=exp(max1-lw1),top2=exp(max2-lw2)
    let exactA = mode==1 ? target.a1 : target.a2
    let exactQ = mode==0 ? target.q : nil
    let ratio = exactA.flatMap { Double($0) }.map { exp((mode==1 ? logA1 : logA2)-log($0)) }
    let qerr = exactQ.map { qE-$0 }
    let ratioText = ratio.map(jsonNumber) ?? "null"
    let qerrText = qerr.map(jsonNumber) ?? "null"
    let summary = "{\"kind\":\"summary\",\"n\":\(n),\"mode\":\(mode),\"samples\":\(batch*totalBatches),\"K1\":\(k1),\"K2\":\(k2),\"invalid\":\(bad),\"log_A1_hat\":\(jsonNumber(logA1)),\"log_A2_hat\":\(jsonNumber(logA2)),\"Q_euclidean_hat\":\(jsonNumber(qE)),\"Q_manhattan_hat\":\(jsonNumber(qR)),\"ESS_A1\":\(jsonNumber(ess1)),\"ESS_A2\":\(jsonNumber(ess2)),\"top_weight_A1\":\(jsonNumber(top1)),\"top_weight_A2\":\(jsonNumber(top2)),\"gpu_seconds\":\(jsonNumber(gpuSeconds)),\"submit_wait_seconds\":\(jsonNumber(wallSeconds)),\"log_factorial_N\":\(jsonNumber(logFactorial(N))),\"A_ratio_to_exact\":\(ratioText),\"Q_error_to_exact\":\(qerrText)}"
    emit(summary)
  }
}
try fh.close()
