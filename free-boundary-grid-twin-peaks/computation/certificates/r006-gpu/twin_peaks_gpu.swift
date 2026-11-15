// Metal batch heuristic for free-boundary P_n square P_n twin peaks.
// It is a construction finder only: CPU independently verifies every saved hit.
import Foundation
import Metal

struct Config { let n: Int; let chains: Int; let rounds: Int; let steps: Int; let seed: UInt64 }
struct GPUParams { var n: UInt32; var count: UInt32; var steps: UInt32; var temperature: UInt32 }

let shader = """
#include <metal_stdlib>
using namespace metal;
struct Params { uint n, count, steps, temperature; };
inline uint nextRand(thread uint &s) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }
inline int energy(device const uint *a, uint n, thread uint &p1, thread uint &p2) {
  uint N=n*n, k=0; p1=N; p2=N;
  for (uint v=0; v<N; ++v) {
    uint x=v % n, y=v / n, z=a[v]; bool peak=true;
    if (x && a[v-1]>z) peak=false;
    if (x+1<n && a[v+1]>z) peak=false;
    if (y && a[v-n]>z) peak=false;
    if (y+1<n && a[v+n]>z) peak=false;
    if (peak) { if(k==0)p1=v; else if(k==1)p2=v; ++k; }
  }
  if (k==2) { int dx=int(p1%n)-int(p2%n), dy=int(p1/n)-int(p2/n); return -abs(dx)-abs(dy); }
  return 10000*abs(int(k)-2);
}
kernel void search(device uint *labels [[buffer(0)]], device uint *seeds [[buffer(1)]],
                   device uint *outK [[buffer(2)]], device uint2 *outP [[buffer(3)]],
                   constant Params &q [[buffer(4)]], device uint *outAccepted [[buffer(5)]], uint gid [[thread_position_in_grid]]) {
  if (gid>=q.count) return;
  uint N=q.n*q.n, base=gid*N; device uint *a=labels+base; uint s=seeds[gid] ^ (gid*747796405u+2891336453u);
  uint p1,p2,acceptedMoves=0; int e=energy(a,q.n,p1,p2);
  for(uint it=0;it<q.steps;++it){
    uint i=nextRand(s)%N, j=nextRand(s)%N; if(i==j)continue;
    uint tmp=a[i];a[i]=a[j];a[j]=tmp;
    uint u,v; int ne=energy(a,q.n,u,v); bool accept=ne<=e;
    if(!accept && q.temperature) { uint r=nextRand(s)%10000u; uint slack=uint(min(ne-e,9999)); accept=(r < q.temperature/(1u+slack/200u)); }
    if(accept){e=ne;p1=u;p2=v;++acceptedMoves;} else {tmp=a[i];a[i]=a[j];a[j]=tmp;}
  }
  uint fk=0; p1=N;p2=N; (void)energy(a,q.n,p1,p2);
  for(uint v=0;v<N;++v){uint x=v%q.n,y=v/q.n,z=a[v];bool peak=true;
    if(x&&a[v-1]>z)peak=false;if(x+1<q.n&&a[v+1]>z)peak=false;if(y&&a[v-q.n]>z)peak=false;if(y+1<q.n&&a[v+q.n]>z)peak=false;
    if(peak)++fk;}
  outK[gid]=fk; outP[gid]=uint2(p1,p2); outAccepted[gid]=acceptedMoves; seeds[gid]=s;
}
"""

@inline(__always) func next(_ s: inout UInt64) -> UInt64 { s ^= s << 7; s ^= s >> 9; return s }
func shuffle<T>(_ a: inout [T], _ s: inout UInt64) { if a.count>1 { for i in stride(from:a.count-1,through:1,by:-1) { let j=Int(next(&s)%UInt64(i+1)); a.swapAt(i,j) } } }
func peakInfo(_ a:[UInt32], _ n:Int)->(Int,[Int]) { var p:[Int]=[]; for v in 0..<a.count { let x=v%n,y=v/n,z=a[v]; var ok=true; if x>0 && a[v-1]>z {ok=false}; if x+1<n && a[v+1]>z {ok=false}; if y>0 && a[v-n]>z {ok=false}; if y+1<n && a[v+n]>z {ok=false}; if ok {p.append(v)} }; return(p.count,p) }
func distance(_ u:Int,_ v:Int,_ n:Int)->Int { abs(u%n-v%n)+abs(u/n-v/n) }

// Two complete rectangles, each ranked strictly downwards from an outer root.
func cutSeed(n:Int, state:inout UInt64)->[UInt32] {
  let N=n*n, cut=n/2; var left:[(Int,Int)]=[], right:[(Int,Int)]=[]
  let rl=Int(next(&state)%UInt64(n)), rr=Int(next(&state)%UInt64(n))
  for y in 0..<n { for x in 0..<n { let v=y*n+x; if x<cut { left.append((abs(y-rl)+x,v)) } else { right.append((abs(y-rr)+(n-1-x),v)) } } }
  left.sort { $0.0 == $1.0 ? $0.1 < $1.1 : $0.0 < $1.0 }; right.sort { $0.0 == $1.0 ? $0.1 < $1.1 : $0.0 < $1.0 }
  var side=Array(repeating:false,count:N); for i in 0..<left.count {side[i]=true}; shuffle(&side,&state)
  var a=Array(repeating:UInt32(0),count:N), i=0,j=0
  for rank in stride(from:N,through:1,by:-1) { if side[N-rank] { a[left[i].1]=UInt32(rank);i+=1 } else { a[right[j].1]=UInt32(rank);j+=1 } }
  return a
}
func randomSeed(n:Int,state:inout UInt64)->[UInt32] { var a=(1...n*n).map(UInt32.init);shuffle(&a,&state);return a }
func jsonEscape(_ a:[UInt32])->String { a.map(String.init).joined(separator:",") }

let args=CommandLine.arguments
let cutOnly=args.contains("--cut-only"), longRun=args.contains("--long")
let requested: [Config] = longRun ? [
  Config(n:24,chains:8192, rounds:64, steps:4096, seed:0xA022_1024),
  Config(n:32,chains:6144, rounds:64, steps:4096, seed:0xA022_1032)
] : [
  Config(n:8, chains:8192, rounds:20, steps:2048, seed:0xA022_0008),
  Config(n:12,chains:8192, rounds:24, steps:2048, seed:0xA022_0012),
  Config(n:16,chains:6144, rounds:28, steps:2048, seed:0xA022_0016),
  Config(n:24,chains:4096, rounds:32, steps:2048, seed:0xA022_0024),
  Config(n:32,chains:3072, rounds:36, steps:2048, seed:0xA022_0032)
]
let only = args.dropFirst().compactMap { Int($0) }; let configs=requested.filter { only.isEmpty || only.contains($0.n) }
guard let device=MTLCreateSystemDefaultDevice(), let queue=device.makeCommandQueue() else { fatalError("Metal device unavailable") }
let library:MTLLibrary; do { library=try device.makeLibrary(source:shader,options:nil) } catch { fatalError("Metal compiler: \(error)") }
guard let fun=library.makeFunction(name:"search"), let pipe=try? device.makeComputePipelineState(function:fun) else { fatalError("pipeline unavailable") }
print("device=\(device.name) maxThreads=\(pipe.maxTotalThreadsPerThreadgroup)")
let base=URL(fileURLWithPath:FileManager.default.currentDirectoryPath); let outURL=base.appendingPathComponent("gpu-search-results.jsonl")
FileManager.default.createFile(atPath:outURL.path,contents:nil); let handle=try FileHandle(forWritingTo:outURL)
for cfg in configs {
  let N=cfg.n*cfg.n; var host=[UInt32]();host.reserveCapacity(cfg.chains*N);var hostState=cfg.seed
  for c in 0..<cfg.chains { host += (cutOnly || c % 4 == 0 ? cutSeed(n:cfg.n,state:&hostState) : randomSeed(n:cfg.n,state:&hostState)) }
  // CPU/GPU agree before any heuristic move, on all seeds in this batch.
  var cpuK=[UInt32]();cpuK.reserveCapacity(cfg.chains);for c in 0..<cfg.chains {cpuK.append(UInt32(peakInfo(Array(host[c*N..<(c+1)*N]),cfg.n).0))}
  let labels=device.makeBuffer(bytes:host,length:host.count*MemoryLayout<UInt32>.stride,options:.storageModeShared)!
  var seeds=(0..<cfg.chains).map { UInt32(truncatingIfNeeded: cfg.seed &+ UInt64($0)*0x9e3779b9) }
  let seedBuf=device.makeBuffer(bytes:&seeds,length:seeds.count*4,options:.storageModeShared)!
  let kBuf=device.makeBuffer(length:cfg.chains*4,options:.storageModeShared)!, pBuf=device.makeBuffer(length:cfg.chains*MemoryLayout<SIMD2<UInt32>>.stride,options:.storageModeShared)!, acceptBuf=device.makeBuffer(length:cfg.chains*4,options:.storageModeShared)!
  var zero=GPUParams(n:UInt32(cfg.n),count:UInt32(cfg.chains),steps:0,temperature:0)
  let par=device.makeBuffer(bytes:&zero,length:MemoryLayout<GPUParams>.stride,options:.storageModeShared)!
  func dispatch(_ params:GPUParams)->Double { par.contents().bindMemory(to:GPUParams.self,capacity:1).pointee=params; let cb=queue.makeCommandBuffer()!;let enc=cb.makeComputeCommandEncoder()!;enc.setComputePipelineState(pipe);enc.setBuffer(labels,offset:0,index:0);enc.setBuffer(seedBuf,offset:0,index:1);enc.setBuffer(kBuf,offset:0,index:2);enc.setBuffer(pBuf,offset:0,index:3);enc.setBuffer(par,offset:0,index:4);enc.setBuffer(acceptBuf,offset:0,index:5);let tg=min(256,pipe.maxTotalThreadsPerThreadgroup);enc.dispatchThreads(MTLSize(width:cfg.chains,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:tg,height:1,depth:1));enc.endEncoding();let t=Date();cb.commit();cb.waitUntilCompleted();if cb.status != .completed { fatalError("GPU command failed \(String(describing:cb.error))") };return Date().timeIntervalSince(t) }
  _=dispatch(zero);let got=kBuf.contents().bindMemory(to:UInt32.self,capacity:cfg.chains);for c in 0..<cfg.chains { if got[c] != cpuK[c] {fatalError("CPU/GPU peak disagreement n=\(cfg.n) chain=\(c): \(cpuK[c]) vs \(got[c])")} }
  var gpuSeconds=0.0, verified=0, far=0, bestDist=0, bestK=Int.max; var acceptedTotal:UInt64=0
  for r in 0..<cfg.rounds {
    gpuSeconds += dispatch(GPUParams(n:UInt32(cfg.n),count:UInt32(cfg.chains),steps:UInt32(cfg.steps),temperature: cutOnly ? 0 : (r<4 ? 2500 : 150)))
    let lab=labels.contents().bindMemory(to:UInt32.self,capacity:host.count), acc=acceptBuf.contents().bindMemory(to:UInt32.self,capacity:cfg.chains)
    acceptedTotal += (0..<cfg.chains).reduce(0) { $0 + UInt64(acc[$1]) }
    for c in 0..<cfg.chains { let k=Int(got[c]);bestK=min(bestK,k);if k==2 {let a=Array(UnsafeBufferPointer(start:lab+c*N,count:N));let info=peakInfo(a,cfg.n);precondition(info.0==2,"CPU mismatch after GPU search");let d=distance(info.1[0],info.1[1],cfg.n);verified+=1;bestDist=max(bestDist,d);if d>=cfg.n/2 {far+=1;if far<=64 {let rec="{\"n\":\(cfg.n),\"round\":\(r),\"chain\":\(c),\"k\":2,\"distance\":\(d),\"labels\":[\(jsonEscape(a))]}\n";try! handle.write(contentsOf:rec.data(using:.utf8)!)} } } }
    let elapsed=String(format:"%.3f",gpuSeconds)
    print("n=\(cfg.n) round=\(r+1)/\(cfg.rounds) gpu_s=\(elapsed) K2_checked=\(verified) far=\(far) bestDist=\(bestDist)")
  }
  let moves=Double(cfg.chains*cfg.rounds*cfg.steps);let rec="{\"summary\":true,\"n\":\(cfg.n),\"chains\":\(cfg.chains),\"rounds\":\(cfg.rounds),\"steps\":\(cfg.steps),\"seed\":\(cfg.seed),\"cut_only\":\(cutOnly),\"gpu_seconds\":\(gpuSeconds),\"moves\":\(moves),\"moves_per_second\":\(moves/gpuSeconds),\"accepted_moves\":\(acceptedTotal),\"cpu_gpu_initial_check\":true,\"cpu_verified_k2_visits\":\(verified),\"cpu_verified_far_visits\":\(far),\"best_distance\":\(bestDist),\"best_k_seen\":\(bestK),\"buffer_bytes\":\(host.count*4+cfg.chains*16)}\n";try! handle.write(contentsOf:rec.data(using:.utf8)!);print(rec)
}
try handle.close()
