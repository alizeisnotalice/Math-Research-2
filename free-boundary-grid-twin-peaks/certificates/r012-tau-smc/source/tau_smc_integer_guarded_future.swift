import Foundation
import Metal

// Original-law tau-stratified particle filter.
// Reference law: a uniformly random descending vertex permutation.  A fixed
// tau forces boundary choices except for its second birth.  All resampling
// weights are integer allowed-set counts.  This is numerical code; PRNG use
// does not constitute an independence proof.

let metal = """
#include <metal_stdlib>
using namespace metal;
struct Par { uint n,N,tau,rank,seedlo,seedhi; };
struct State { uint m0,m1,m2,m3,m4,root1,root2; };
inline ulong seed64(constant Par&p){ return (ulong(p.seedhi)<<32)|ulong(p.seedlo); }
inline ulong next64(thread ulong&s){s+=0x9e3779b97f4a7c15ul;ulong z=s;z=(z^(z>>30))*0xbf58476d1ce4e5b9ul;z=(z^(z>>27))*0x94d049bb133111ebul;return z^(z>>31);}
inline ulong childkey(constant Par&p,uint child,uint domain){return seed64(p) ^ (ulong(p.rank)<<32) ^ (ulong(child)*0x9e3779b97f4a7c15ul) ^ (ulong(domain)*0xd1b54a32d192ed03ul);}
inline ulong below(thread ulong&s,ulong b){ulong threshold=(0ul-b)%b;while(true){ulong x=next64(s);if(x>=threshold)return x%b;}}
inline bool has(thread State&s,uint v){uint q=v>>5,b=v&31,a=q==0?s.m0:q==1?s.m1:q==2?s.m2:q==3?s.m3:s.m4;return (a&(1u<<b))!=0;}
inline void add(thread State&s,uint v){uint q=v>>5,b=v&31;if(q==0)s.m0|=1u<<b;else if(q==1)s.m1|=1u<<b;else if(q==2)s.m2|=1u<<b;else if(q==3)s.m3|=1u<<b;else s.m4|=1u<<b;}
inline bool touch(thread State&s,uint v,uint n){uint x=v%n,y=v/n;return(x&&has(s,v-1))||(x+1<n&&has(s,v+1))||(y&&has(s,v-n))||(y+1<n&&has(s,v+n));}
inline bool ok(thread State&s,uint v,constant Par&p){if(has(s,v))return false;if(p.rank==1)return true;bool t=touch(s,v,p.n);return p.rank==p.tau ? !t : t;}
kernel void counts(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],constant Par&p[[buffer(2)]],uint i[[thread_position_in_grid]]){thread State s=src[i];uint z=0;for(uint v=0;v<p.N;v++)if(ok(s,v,p))z++;a[i]=z;}
kernel void groupSums(device uint*a[[buffer(0)]],device ulong*out[[buffer(1)]],constant uint&count[[buffer(2)]],uint g[[thread_position_in_grid]]){uint lo=g*256,hi=min(lo+256,count);ulong s=0;for(uint i=lo;i<hi;i++)s+=ulong(a[i]);out[g]=s;}
kernel void resprop(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device uint*group[[buffer(2)]],device ulong*resid[[buffer(3)]],device State*dst[[buffer(4)]],device uint*parentOut[[buffer(5)]],constant Par&p[[buffer(6)]],constant uint&count[[buffer(7)]],uint child[[thread_position_in_grid]]){if(child>=count)return;uint lo=group[child]*256,hi=min(lo+256,count);ulong r=resid[child],c=0;uint par=lo;for(uint k=lo;k<hi;k++){uint w=a[k];if(r<c+ulong(w)){par=k;break;}c+=ulong(w);}thread State s=src[par];parentOut[child]=par;ulong rng=childkey(p,child,2);uint num=a[par];ulong pick=below(rng,ulong(num));uint seen=0,v=0;for(v=0;v<p.N;v++)if(ok(s,v,p)){if(seen==pick)break;seen++;}if(p.rank==1)s.root1=v;if(p.rank==p.tau)s.root2=v;add(s,v);dst[child]=s;}
kernel void terminalGroups(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*moment[[buffer(2)]],constant Par&p[[buffer(3)]],constant uint&count[[buffer(4)]],uint g[[thread_position_in_grid]]){uint lo=g*256,hi=min(lo+256,count);ulong s=0;for(uint i=lo;i<hi;i++){thread State x=src[i];uint r2=x.root2;if(p.rank==p.tau){for(uint v=0;v<p.N;v++)if(ok(x,v,p)){r2=v;break;}}int dx=int(x.root1%p.n)-int(r2%p.n),dy=int(x.root1/p.n)-int(r2/p.n);s+=ulong(a[i])*ulong(dx*dx+dy*dy);}moment[g]=s;}
"""

struct Par { var n, N, tau, rank, seedlo, seedhi: UInt32 }
struct State { var m0,m1,m2,m3,m4,root1,root2: UInt32 }
// This tracks only direct parents in one resampling step; it is not a
// multi-generation genealogical-ancestor count.
struct Rep { let logZ: Double; let z: Double; let u: Double; let minOneStepParents: Int; let dead: Bool; let gpuSeconds: Double; let dispatchSeconds: Double }
struct DispatchTiming { let gpu: Double; let synchronous: Double }

let cache = "/Volumes/data/CodexTemp/cache/swift-modules"
let device = MTLCreateSystemDefaultDevice()!
let library = try device.makeLibrary(source: metal, options: nil)
let countPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "counts")!)
let groupPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "groupSums")!)
let propPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "resprop")!)
let termPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "terminalGroups")!)
let queue = device.makeCommandQueue()!

func mix(_ x: UInt64) -> UInt64 { var z = x &+ 0x9e3779b97f4a7c15; z = (z ^ (z >> 30)) &* 0xbf58476d1ce4e5b9; z = (z ^ (z >> 27)) &* 0x94d049bb133111eb; return z ^ (z >> 31) }
func below(_ state: inout UInt64, _ b: UInt64) -> UInt64 { let threshold = (0 &- b) % b; while true { state &+= 0x9e3779b97f4a7c15; let x = mix(state &- 0x9e3779b97f4a7c15); if x >= threshold { return x % b } } }
func childKey(_ seed: UInt64, _ rank: Int, _ child: Int, _ domain: UInt64) -> UInt64 { seed ^ (UInt64(rank) << 32) ^ (UInt64(child) &* 0x9e3779b97f4a7c15) ^ (domain &* 0xd1b54a32d192ed03) }
func groups(_ p: Int) -> Int { (p + 255) / 256 }
func dispatch(_ pipe: MTLComputePipelineState, _ buffers: [(MTLBuffer, Int)], _ threads: Int) -> DispatchTiming { let begin=Date(); let cb = queue.makeCommandBuffer()!; let e = cb.makeComputeCommandEncoder()!; e.setComputePipelineState(pipe); for (b,i) in buffers { e.setBuffer(b, offset: 0, index: i) }; e.dispatchThreads(MTLSize(width: threads, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: 256, height: 1, depth: 1)); e.endEncoding(); cb.commit(); cb.waitUntilCompleted(); precondition(cb.status == .completed, cb.error?.localizedDescription ?? "Metal failure"); return DispatchTiming(gpu: max(0,cb.gpuEndTime-cb.gpuStartTime), synchronous: -begin.timeIntervalSinceNow) }
func summarize(_ values: [UInt64]) -> UInt64 { values.reduce(0, &+) }
func uniqueCount(_ x: UnsafePointer<UInt32>, _ n: Int) -> Int { Set(UnsafeBufferPointer(start: x, count: n)).count }

func runRep(n: Int, tau: Int, P: Int, seed: UInt64) -> Rep {
    let N = n*n, G = groups(P)
    var par = Par(n: UInt32(n), N: UInt32(N), tau: UInt32(tau), rank: 1, seedlo: UInt32(truncatingIfNeeded: seed), seedhi: UInt32(truncatingIfNeeded: seed >> 32))
    let parBuf = device.makeBuffer(bytes: &par, length: MemoryLayout<Par>.stride, options: .storageModeShared)!
    let countBuf = device.makeBuffer(length: P*MemoryLayout<UInt32>.stride, options: .storageModeShared)!
    let groupBuf = device.makeBuffer(length: G*MemoryLayout<UInt64>.stride, options: .storageModeShared)!
    let momentBuf = device.makeBuffer(length: G*MemoryLayout<UInt64>.stride, options: .storageModeShared)!
    let parentBuf = device.makeBuffer(length: P*MemoryLayout<UInt32>.stride, options: .storageModeShared)!
    let targetGroup = device.makeBuffer(length: P*MemoryLayout<UInt32>.stride, options: .storageModeShared)!
    let targetResid = device.makeBuffer(length: P*MemoryLayout<UInt64>.stride, options: .storageModeShared)!
    let stateBytes = P*MemoryLayout<State>.stride
    let x = device.makeBuffer(length: stateBytes, options: .storageModeShared)!
    let y = device.makeBuffer(length: stateBytes, options: .storageModeShared)!
    x.contents().bindMemory(to: State.self, capacity: P).initialize(repeating: State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0), count: P)
    var src=x, dst=y, logZ=0.0, minOneStepParents=P, gpuSeconds=0.0, dispatchSeconds=0.0
    var pCount=UInt32(P)
    let pCountBuf=device.makeBuffer(bytes:&pCount,length:4,options:.storageModeShared)!
    for rank in 1...N {
        par.rank=UInt32(rank); parBuf.contents().bindMemory(to:Par.self,capacity:1).pointee=par
        var timing=dispatch(countPipe, [(src,0),(countBuf,1),(parBuf,2)], P); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        timing=dispatch(groupPipe, [(countBuf,0),(groupBuf,1),(pCountBuf,2)], G); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        let sums=Array(UnsafeBufferPointer(start:groupBuf.contents().bindMemory(to:UInt64.self,capacity:G),count:G))
        let total=summarize(sums), remaining=UInt64(N-rank+1)
        if total == 0 { return Rep(logZ: -.infinity,z:0,u:0,minOneStepParents:0,dead:true,gpuSeconds:gpuSeconds,dispatchSeconds:dispatchSeconds) }
        logZ += log(Double(total)/(Double(P)*Double(remaining)))
        if rank == N {
            timing=dispatch(termPipe, [(src,0),(countBuf,1),(momentBuf,2),(parBuf,3),(pCountBuf,4)], G); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
            let raw=summarize(Array(UnsafeBufferPointer(start:momentBuf.contents().bindMemory(to:UInt64.self,capacity:G),count:G)))
            let z=exp(logZ), u=z*Double(raw)/Double(total)/Double(N)
            return Rep(logZ:logZ,z:z,u:u,minOneStepParents:minOneStepParents,dead:false,gpuSeconds:gpuSeconds,dispatchSeconds:dispatchSeconds)
        }
        var cumulative=[UInt64](); cumulative.reserveCapacity(G); var c:UInt64=0; for v in sums { c &+= v; cumulative.append(c) }
        let gp=targetGroup.contents().bindMemory(to:UInt32.self,capacity:P)
        let rp=targetResid.contents().bindMemory(to:UInt64.self,capacity:P)
        for child in 0..<P { var r=childKey(seed,rank,child,1); let u=below(&r,total); var lo=0,hi=G-1; while lo<hi { let mid=(lo+hi)/2; if u<cumulative[mid] {hi=mid} else {lo=mid+1} }; let prev=lo==0 ? 0:cumulative[lo-1]; gp[child]=UInt32(lo); rp[child]=u-prev }
        timing=dispatch(propPipe, [(src,0),(countBuf,1),(targetGroup,2),(targetResid,3),(dst,4),(parentBuf,5),(parBuf,6),(pCountBuf,7)], P); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        minOneStepParents=min(minOneStepParents,uniqueCount(parentBuf.contents().bindMemory(to:UInt32.self,capacity:P),P)); swap(&src,&dst)
    }
    fatalError("unreachable")
}

let args=CommandLine.arguments
guard args.count==7 else { fatalError("usage: tau_smc_integer n particles replicates tauLo tauHi output.jsonl") }
let n=Int(args[1])!, P=Int(args[2])!, R=Int(args[3])!, lo=Int(args[4])!, hi=Int(args[5])!
guard (2...12).contains(n), P > 0, R > 0, lo >= 2, hi >= lo, hi <= n*n else {
    fatalError("require 2<=n<=12, particles>0, replicates>0, and 2<=tauLo<=tauHi<=n*n")
}
let out=URL(fileURLWithPath:args[6]); FileManager.default.createFile(atPath:out.path, contents:nil); let fh=try FileHandle(forWritingTo:out); let wallStart=Date()
func emit(_ s:String){try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!)}
emit("{\"kind\":\"metadata\",\"device\":\"\(device.name)\",\"n\":\(n),\"particles\":\(P),\"replicates\":\(R),\"tau_lo\":\(lo),\"tau_hi\":\(hi),\"method\":\"integer group resampling; terminal no-resample\",\"timing\":\"gpu_hardware_seconds=sum Metal GPU timestamps; synchronous_dispatch_seconds=submit+wait on host\",\"genealogy_diagnostic\":\"min_one_step_distinct_parents is not a cross-generation ancestor count\"}")
for tau in lo...hi { for rep in 0..<R { let seed=0xA022000000000000 &+ UInt64(n*1000000+tau*10000+rep); let q=runRep(n:n,tau:tau,P:P,seed:seed); let logZJSON = q.dead ? "null" : "\(q.logZ)"; emit("{\"kind\":\"replicate\",\"tau\":\(tau),\"replicate\":\(rep),\"logZ\":\(logZJSON),\"Z\":\(q.z),\"U\":\(q.u),\"min_one_step_distinct_parents\":\(q.minOneStepParents),\"dead\":\(q.dead),\"gpu_hardware_seconds\":\(q.gpuSeconds),\"synchronous_dispatch_seconds\":\(q.dispatchSeconds)}") } }
emit("{\"kind\":\"run_end\",\"wall_seconds\":\(-wallStart.timeIntervalSinceNow)}")
try fh.close()
