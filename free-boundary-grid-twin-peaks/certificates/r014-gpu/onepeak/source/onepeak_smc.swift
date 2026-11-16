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
// GPU parent selection.  cdf is supplied by the host from exact group sums;
// domain 1 and domain 2 calls deliberately match the prior GPU kernel.
kernel void respropGPU(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*cdf[[buffer(2)]],constant ulong&total[[buffer(3)]],device State*dst[[buffer(4)]],device uint*parentOut[[buffer(5)]],constant Par&p[[buffer(6)]],constant uint&count[[buffer(7)]],uint child[[thread_position_in_grid]]){if(child>=count)return;ulong rng0=childkey(p,child,1);ulong u=below(rng0,total);uint G=(count+255)/256,loG=0,hiG=G-1;while(loG<hiG){uint mid=(loG+hiG)/2;if(u<cdf[mid])hiG=mid;else loG=mid+1;}ulong prev=loG==0?0:cdf[loG-1],r=u-prev;uint lo=loG*256,hi=min(lo+256,count);ulong c=0;uint par=lo;for(uint k=lo;k<hi;k++){uint w=a[k];if(r<c+ulong(w)){par=k;break;}c+=ulong(w);}thread State s=src[par];parentOut[child]=par;ulong rng=childkey(p,child,2);uint num=a[par];ulong pick=below(rng,ulong(num));uint seen=0,v=0;for(v=0;v<p.N;v++)if(ok(s,v,p)){if(seen==pick)break;seen++;}if(p.rank==1)s.root1=v;if(p.rank==p.tau)s.root2=v;add(s,v);dst[child]=s;}
kernel void terminalGroups(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*moment[[buffer(2)]],constant Par&p[[buffer(3)]],constant uint&count[[buffer(4)]],uint g[[thread_position_in_grid]]){uint lo=g*256,hi=min(lo+256,count);ulong s=0;for(uint i=lo;i<hi;i++){thread State x=src[i];uint r2=x.root2;if(p.rank==p.tau){for(uint v=0;v<p.N;v++)if(ok(x,v,p)){r2=v;break;}}int dx=int(x.root1%p.n)-int(r2%p.n),dy=int(x.root1/p.n)-int(r2/p.n);s+=ulong(a[i])*ulong(dx*dx+dy*dy);}moment[g]=s;}
// r013: all rank-to-rank dependencies stay on GPU in one command buffer.
kernel void countsFull(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],constant Par&p[[buffer(2)]],device uint*dead[[buffer(3)]],uint i[[thread_position_in_grid]]){if(dead[0]){a[i]=0;return;}thread State s=src[i];uint z=0;for(uint v=0;v<p.N;v++)if(ok(s,v,p))z++;a[i]=z;}
kernel void cdfFull(device ulong*groups[[buffer(0)]],device ulong*cdf[[buffer(1)]],device ulong*totals[[buffer(2)]],device uint*dead[[buffer(3)]],constant uint&G[[buffer(4)]],constant uint&rankIndex[[buffer(5)]],uint i[[thread_position_in_grid]]){if(i)return;ulong s=0;for(uint j=0;j<G;j++){s+=groups[j];cdf[j]=s;}totals[rankIndex]=s;if(s==0)dead[0]=1;}
kernel void respropFull(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*cdf[[buffer(2)]],device ulong*totals[[buffer(3)]],device State*dst[[buffer(4)]],constant Par&p[[buffer(5)]],constant uint&count[[buffer(6)]],constant uint&rankIndex[[buffer(7)]],device uint*dead[[buffer(8)]],constant uint&stratified[[buffer(9)]],uint child[[thread_position_in_grid]]){if(child>=count)return;if(dead[0]){dst[child]=src[child];return;}ulong total=totals[rankIndex];ulong rng0=childkey(p,child,1);ulong t=below(rng0,total);ulong u=stratified ? (ulong(child)*total+t)/ulong(count) : t;uint G=(count+255)/256,loG=0,hiG=G-1;while(loG<hiG){uint mid=(loG+hiG)/2;if(u<cdf[mid])hiG=mid;else loG=mid+1;}ulong prev=loG==0?0:cdf[loG-1],r=u-prev;uint lo=loG*256,hi=min(lo+256,count);ulong c=0;uint par=lo;for(uint k=lo;k<hi;k++){uint w=a[k];if(r<c+ulong(w)){par=k;break;}c+=ulong(w);}thread State s=src[par];ulong rng=childkey(p,child,2);uint num=a[par];ulong pick=below(rng,ulong(num));uint seen=0,v=0;for(v=0;v<p.N;v++)if(ok(s,v,p)){if(seen==pick)break;seen++;}if(p.rank==1)s.root1=v;if(p.rank==p.tau)s.root2=v;add(s,v);dst[child]=s;}
kernel void terminalFull(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*moment[[buffer(2)]],constant Par&p[[buffer(3)]],constant uint&count[[buffer(4)]],device uint*dead[[buffer(5)]],uint g[[thread_position_in_grid]]){if(dead[0]){moment[g]=0;return;}uint lo=g*256,hi=min(lo+256,count);ulong s=0;for(uint i=lo;i<hi;i++){thread State x=src[i];uint r2=x.root2;if(p.rank==p.tau){for(uint v=0;v<p.N;v++)if(ok(x,v,p)){r2=v;break;}}int dx=int(x.root1%p.n)-int(r2%p.n),dy=int(x.root1/p.n)-int(r2/p.n);s+=ulong(a[i])*ulong(dx*dx+dy*dy);}moment[g]=s;}
"""

struct Par { var n, N, tau, rank, seedlo, seedhi: UInt32 }
struct State { var m0,m1,m2,m3,m4,root1,root2: UInt32 }
struct Rep { let logZ: Double; let z: Double; let u: Double; let minAnc: Int; let dead: Bool; let gpuSeconds: Double; let dispatchSeconds: Double }
struct DispatchTiming { let gpu: Double; let synchronous: Double }

let cache = "/Volumes/data/CodexTemp/cache/swift-modules"
let device = MTLCreateSystemDefaultDevice()!
let library = try device.makeLibrary(source: metal, options: nil)
let countPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "counts")!)
let groupPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "groupSums")!)
let propPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "respropGPU")!)
let termPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "terminalGroups")!)
let fullCountPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "countsFull")!)
let fullCDFPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "cdfFull")!)
let fullPropPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "respropFull")!)
let fullTermPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "terminalFull")!)
let queue = device.makeCommandQueue()!

func mix(_ x: UInt64) -> UInt64 { var z = x &+ 0x9e3779b97f4a7c15; z = (z ^ (z >> 30)) &* 0xbf58476d1ce4e5b9; z = (z ^ (z >> 27)) &* 0x94d049bb133111eb; return z ^ (z >> 31) }
func below(_ state: inout UInt64, _ b: UInt64) -> UInt64 { let threshold = (0 &- b) % b; while true { state &+= 0x9e3779b97f4a7c15; let x = mix(state &- 0x9e3779b97f4a7c15); if x >= threshold { return x % b } } }
func childKey(_ seed: UInt64, _ rank: Int, _ child: Int, _ domain: UInt64) -> UInt64 { seed ^ (UInt64(rank) << 32) ^ (UInt64(child) &* 0x9e3779b97f4a7c15) ^ (domain &* 0xd1b54a32d192ed03) }
func groups(_ p: Int) -> Int { (p + 255) / 256 }
func dispatch(_ pipe: MTLComputePipelineState, _ buffers: [(MTLBuffer, Int)], _ threads: Int) -> DispatchTiming { let begin=Date(); let cb = queue.makeCommandBuffer()!; let e = cb.makeComputeCommandEncoder()!; e.setComputePipelineState(pipe); for (b,i) in buffers { e.setBuffer(b, offset: 0, index: i) }; e.dispatchThreads(MTLSize(width: threads, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: 256, height: 1, depth: 1)); e.endEncoding(); cb.commit(); cb.waitUntilCompleted(); precondition(cb.status == .completed, cb.error?.localizedDescription ?? "Metal failure"); return DispatchTiming(gpu: max(0,cb.gpuEndTime-cb.gpuStartTime), synchronous: -begin.timeIntervalSinceNow) }
func encode(_ cb: MTLCommandBuffer, _ pipe: MTLComputePipelineState, _ buffers: [(MTLBuffer, Int)], _ threads: Int) { let e=cb.makeComputeCommandEncoder()!;e.setComputePipelineState(pipe);for (b,i) in buffers{e.setBuffer(b,offset:0,index:i)};e.dispatchThreads(MTLSize(width:threads,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1));e.endEncoding() }
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
    let cdfBuf = device.makeBuffer(length: G*MemoryLayout<UInt64>.stride, options: .storageModeShared)!
    var totalStorage: UInt64 = 0
    let totalBuf = device.makeBuffer(bytes: &totalStorage, length: MemoryLayout<UInt64>.stride, options: .storageModeShared)!
    let stateBytes = P*MemoryLayout<State>.stride
    let x = device.makeBuffer(length: stateBytes, options: .storageModeShared)!
    let y = device.makeBuffer(length: stateBytes, options: .storageModeShared)!
    x.contents().bindMemory(to: State.self, capacity: P).initialize(repeating: State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0), count: P)
    var src=x, dst=y, logZ=0.0, minAnc = -1, gpuSeconds=0.0, dispatchSeconds=0.0
    var pCount=UInt32(P)
    let pCountBuf=device.makeBuffer(bytes:&pCount,length:4,options:.storageModeShared)!
    for rank in 1...N {
        par.rank=UInt32(rank); parBuf.contents().bindMemory(to:Par.self,capacity:1).pointee=par
        var timing=dispatch(countPipe, [(src,0),(countBuf,1),(parBuf,2)], P); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        timing=dispatch(groupPipe, [(countBuf,0),(groupBuf,1),(pCountBuf,2)], G); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        let sums=Array(UnsafeBufferPointer(start:groupBuf.contents().bindMemory(to:UInt64.self,capacity:G),count:G))
        let total=summarize(sums), remaining=UInt64(N-rank+1)
        if total == 0 { return Rep(logZ: -.infinity,z:0,u:0,minAnc:0,dead:true,gpuSeconds:gpuSeconds,dispatchSeconds:dispatchSeconds) }
        logZ += log(Double(total)/(Double(P)*Double(remaining)))
        if rank == N {
            timing=dispatch(termPipe, [(src,0),(countBuf,1),(momentBuf,2),(parBuf,3),(pCountBuf,4)], G); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
            let raw=summarize(Array(UnsafeBufferPointer(start:momentBuf.contents().bindMemory(to:UInt64.self,capacity:G),count:G)))
            let z=exp(logZ), u=z*Double(raw)/Double(total)/Double(N)
            return Rep(logZ:logZ,z:z,u:u,minAnc:minAnc,dead:false,gpuSeconds:gpuSeconds,dispatchSeconds:dispatchSeconds)
        }
        var cumulative=[UInt64](); cumulative.reserveCapacity(G); var c:UInt64=0; for v in sums { c &+= v; cumulative.append(c) }
        cumulative.withUnsafeBufferPointer { src in cdfBuf.contents().copyMemory(from: src.baseAddress!, byteCount: G*MemoryLayout<UInt64>.stride) }
        totalStorage=total; totalBuf.contents().bindMemory(to:UInt64.self,capacity:1).pointee=totalStorage
        timing=dispatch(propPipe, [(src,0),(countBuf,1),(cdfBuf,2),(totalBuf,3),(dst,4),(parentBuf,5),(parBuf,6),(pCountBuf,7)], P); gpuSeconds += timing.gpu; dispatchSeconds += timing.synchronous
        swap(&src,&dst)
    }
    fatalError("unreachable")
}

func runRepFull(n:Int,tau:Int,P:Int,seed:UInt64,stratified:Bool=false)->Rep {
    let N=n*n,G=groups(P),stateBytes=P*MemoryLayout<State>.stride
    var pc=UInt32(P),gg=UInt32(G),dead:UInt32=0,smode:UInt32=stratified ? 1:0
    let pcBuf=device.makeBuffer(bytes:&pc,length:4,options:.storageModeShared)!,gBuf=device.makeBuffer(bytes:&gg,length:4,options:.storageModeShared)!,deadBuf=device.makeBuffer(bytes:&dead,length:4,options:.storageModeShared)!,smodeBuf=device.makeBuffer(bytes:&smode,length:4,options:.storageModeShared)!
    let countBuf=device.makeBuffer(length:P*4,options:.storageModeShared)!,groupBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,cdfBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,totals=device.makeBuffer(length:N*8,options:.storageModeShared)!,moment=device.makeBuffer(length:G*8,options:.storageModeShared)!
    let x=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,y=device.makeBuffer(length:stateBytes,options:.storageModeShared)!;x.contents().bindMemory(to:State.self,capacity:P).initialize(repeating:State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0),count:P)
    var parBuffers=[MTLBuffer](),rankBuffers=[MTLBuffer]();for rank in 1...N{var par=Par(n:UInt32(n),N:UInt32(N),tau:UInt32(tau),rank:UInt32(rank),seedlo:UInt32(truncatingIfNeeded:seed),seedhi:UInt32(truncatingIfNeeded:seed>>32));var ri=UInt32(rank-1);parBuffers.append(device.makeBuffer(bytes:&par,length:MemoryLayout<Par>.stride,options:.storageModeShared)!);rankBuffers.append(device.makeBuffer(bytes:&ri,length:4,options:.storageModeShared)!)}
    let begin=Date(),cb=queue.makeCommandBuffer()!;var src=x,dst=y
    for r in 0..<N { encode(cb,fullCountPipe,[(src,0),(countBuf,1),(parBuffers[r],2),(deadBuf,3)],P);encode(cb,groupPipe,[(countBuf,0),(groupBuf,1),(pcBuf,2)],G);encode(cb,fullCDFPipe,[(groupBuf,0),(cdfBuf,1),(totals,2),(deadBuf,3),(gBuf,4),(rankBuffers[r],5)],1);if r+1==N { encode(cb,fullTermPipe,[(src,0),(countBuf,1),(moment,2),(parBuffers[r],3),(pcBuf,4),(deadBuf,5)],G) } else { encode(cb,fullPropPipe,[(src,0),(countBuf,1),(cdfBuf,2),(totals,3),(dst,4),(parBuffers[r],5),(pcBuf,6),(rankBuffers[r],7),(deadBuf,8),(smodeBuf,9)],P);swap(&src,&dst) } }
    cb.commit();cb.waitUntilCompleted();precondition(cb.status == .completed,cb.error?.localizedDescription ?? "full GPU failure");let wall = -begin.timeIntervalSinceNow
    let ts=Array(UnsafeBufferPointer(start:totals.contents().bindMemory(to:UInt64.self,capacity:N),count:N));if let z=ts.firstIndex(of:0){return Rep(logZ:-.infinity,z:0,u:0,minAnc:-1,dead:true,gpuSeconds:cb.gpuEndTime-cb.gpuStartTime,dispatchSeconds:wall)}
    var logZ=0.0;for (r,t) in ts.enumerated(){logZ += log(Double(t)/(Double(P)*Double(N-r)))};let raw=summarize(Array(UnsafeBufferPointer(start:moment.contents().bindMemory(to:UInt64.self,capacity:G),count:G))),z=exp(logZ),u=z*Double(raw)/Double(ts[N-1])/Double(N);return Rep(logZ:logZ,z:z,u:u,minAnc:-1,dead:false,gpuSeconds:cb.gpuEndTime-cb.gpuStartTime,dispatchSeconds:wall)
}

let args=CommandLine.arguments
guard args.count==6, let seedDomain=UInt64(args[4],radix:16) else {
    fatalError("usage: onepeak_smc n particles replicates seedDomainHex output.jsonl")
}
let n=Int(args[1])!, P=Int(args[2])!, R=Int(args[3])!, tau=n*n+1
precondition(n >= 2 && n <= 12 && P > 0 && R > 0, "invalid n/P/R")
precondition(UInt64(P) <= UInt64.max / UInt64(P) / UInt64(n*n), "P*T may overflow UInt64")
let out=URL(fileURLWithPath:args[5]); FileManager.default.createFile(atPath:out.path, contents:nil); let fh=try FileHandle(forWritingTo:out); let wallStart=Date()
func emit(_ s:String){try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!)}
emit("{\"kind\":\"metadata\",\"device\":\"\(device.name)\",\"n\":\(n),\"particles\":\(P),\"replicates\":\(R),\"tau\":\(tau),\"seed_domain_hex\":\"\(String(seedDomain,radix:16))\",\"method\":\"one-peak: tau=N+1, single-command-buffer integer-stratified; rank1=V and every later rank=boundary\",\"estimand\":\"Z=A1(n)/(n*n)!\",\"timing\":\"gpu_hardware_seconds=sum Metal GPU timestamps; synchronous_dispatch_seconds=submit+wait on host\"}")
for rep in 0..<R {
    let seed=seedDomain &+ UInt64(n*1_000_000+rep)
    let q: Rep = autoreleasepool { runRepFull(n:n,tau:tau,P:P,seed:seed,stratified:true) }
    let logZJSON = q.dead ? "null" : "\(q.logZ)"
    emit("{\"kind\":\"replicate\",\"tau\":\(tau),\"replicate\":\(rep),\"seed_hex\":\"\(String(seed,radix:16))\",\"logZ\":\(logZJSON),\"Z\":\(q.z),\"dead\":\(q.dead),\"gpu_hardware_seconds\":\(q.gpuSeconds),\"synchronous_dispatch_seconds\":\(q.dispatchSeconds)}")
}
emit("{\"kind\":\"run_end\",\"wall_seconds\":\(-wallStart.timeIntervalSinceNow)}")
try fh.close()
