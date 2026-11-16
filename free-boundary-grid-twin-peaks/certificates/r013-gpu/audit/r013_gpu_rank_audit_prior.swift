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
kernel void respropFull(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*cdf[[buffer(2)]],device ulong*totals[[buffer(3)]],device State*dst[[buffer(4)]],constant Par&p[[buffer(5)]],constant uint&count[[buffer(6)]],constant uint&rankIndex[[buffer(7)]],device uint*dead[[buffer(8)]],constant uint&stratified[[buffer(9)]],device uint*parentOut[[buffer(10)]],uint child[[thread_position_in_grid]]){if(child>=count)return;if(dead[0]){parentOut[child]=0xffffffffu;dst[child]=src[child];return;}ulong total=totals[rankIndex];ulong rng0=childkey(p,child,1);ulong t=below(rng0,total);ulong u=stratified ? (ulong(child)*total+t)/ulong(count) : t;uint G=(count+255)/256,loG=0,hiG=G-1;while(loG<hiG){uint mid=(loG+hiG)/2;if(u<cdf[mid])hiG=mid;else loG=mid+1;}ulong prev=loG==0?0:cdf[loG-1],r=u-prev;uint lo=loG*256,hi=min(lo+256,count);ulong c=0;uint par=lo;for(uint k=lo;k<hi;k++){uint w=a[k];if(r<c+ulong(w)){par=k;break;}c+=ulong(w);}parentOut[child]=par;thread State s=src[par];ulong rng=childkey(p,child,2);uint num=a[par];ulong pick=below(rng,ulong(num));uint seen=0,v=0;for(v=0;v<p.N;v++)if(ok(s,v,p)){if(seen==pick)break;seen++;}if(p.rank==1)s.root1=v;if(p.rank==p.tau)s.root2=v;add(s,v);dst[child]=s;}
kernel void terminalFull(device State*src[[buffer(0)]],device uint*a[[buffer(1)]],device ulong*moment[[buffer(2)]],constant Par&p[[buffer(3)]],constant uint&count[[buffer(4)]],device uint*dead[[buffer(5)]],uint g[[thread_position_in_grid]]){if(dead[0]){moment[g]=0;return;}uint lo=g*256,hi=min(lo+256,count);ulong s=0;for(uint i=lo;i<hi;i++){thread State x=src[i];uint r2=x.root2;if(p.rank==p.tau){for(uint v=0;v<p.N;v++)if(ok(x,v,p)){r2=v;break;}}int dx=int(x.root1%p.n)-int(r2%p.n),dy=int(x.root1/p.n)-int(r2/p.n);s+=ulong(a[i])*ulong(dx*dx+dy*dy);}moment[g]=s;}
// Audit-only observers: these copy kernel-produced data after each rank.  They
// do not affect the production count/CDF/selection/mutation kernels.
kernel void copyU32(device uint*src[[buffer(0)]],device uint*dst[[buffer(1)]],uint i[[thread_position_in_grid]]){dst[i]=src[i];}
kernel void copyU64(device ulong*src[[buffer(0)]],device ulong*dst[[buffer(1)]],uint i[[thread_position_in_grid]]){dst[i]=src[i];}
kernel void copyState(device State*src[[buffer(0)]],device State*dst[[buffer(1)]],uint i[[thread_position_in_grid]]){dst[i]=src[i];}
kernel void setDead(device uint*dead[[buffer(0)]],uint i[[thread_position_in_grid]]){if(i==0)dead[0]=1;}
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
let copyU32Pipe = try device.makeComputePipelineState(function: library.makeFunction(name: "copyU32")!)
let copyU64Pipe = try device.makeComputePipelineState(function: library.makeFunction(name: "copyU64")!)
let copyStatePipe = try device.makeComputePipelineState(function: library.makeFunction(name: "copyState")!)
let setDeadPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "setDead")!)
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
    let countBuf=device.makeBuffer(length:P*4,options:.storageModeShared)!,groupBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,cdfBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,totals=device.makeBuffer(length:N*8,options:.storageModeShared)!,moment=device.makeBuffer(length:G*8,options:.storageModeShared)!,parentObs=device.makeBuffer(length:P*4,options:.storageModeShared)!
    let x=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,y=device.makeBuffer(length:stateBytes,options:.storageModeShared)!;x.contents().bindMemory(to:State.self,capacity:P).initialize(repeating:State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0),count:P)
    var parBuffers=[MTLBuffer](),rankBuffers=[MTLBuffer]();for rank in 1...N{var par=Par(n:UInt32(n),N:UInt32(N),tau:UInt32(tau),rank:UInt32(rank),seedlo:UInt32(truncatingIfNeeded:seed),seedhi:UInt32(truncatingIfNeeded:seed>>32));var ri=UInt32(rank-1);parBuffers.append(device.makeBuffer(bytes:&par,length:MemoryLayout<Par>.stride,options:.storageModeShared)!);rankBuffers.append(device.makeBuffer(bytes:&ri,length:4,options:.storageModeShared)!)}
    let begin=Date(),cb=queue.makeCommandBuffer()!;var src=x,dst=y
    for r in 0..<N { encode(cb,fullCountPipe,[(src,0),(countBuf,1),(parBuffers[r],2),(deadBuf,3)],P);encode(cb,groupPipe,[(countBuf,0),(groupBuf,1),(pcBuf,2)],G);encode(cb,fullCDFPipe,[(groupBuf,0),(cdfBuf,1),(totals,2),(deadBuf,3),(gBuf,4),(rankBuffers[r],5)],1);if r+1==N { encode(cb,fullTermPipe,[(src,0),(countBuf,1),(moment,2),(parBuffers[r],3),(pcBuf,4),(deadBuf,5)],G) } else { encode(cb,fullPropPipe,[(src,0),(countBuf,1),(cdfBuf,2),(totals,3),(dst,4),(parBuffers[r],5),(pcBuf,6),(rankBuffers[r],7),(deadBuf,8),(smodeBuf,9),(parentObs,10)],P);swap(&src,&dst) } }
    cb.commit();cb.waitUntilCompleted();precondition(cb.status == .completed,cb.error?.localizedDescription ?? "full GPU failure");let wall = -begin.timeIntervalSinceNow
    let ts=Array(UnsafeBufferPointer(start:totals.contents().bindMemory(to:UInt64.self,capacity:N),count:N));if let z=ts.firstIndex(of:0){return Rep(logZ:-.infinity,z:0,u:0,minAnc:-1,dead:true,gpuSeconds:cb.gpuEndTime-cb.gpuStartTime,dispatchSeconds:wall)}
    var logZ=0.0;for (r,t) in ts.enumerated(){logZ += log(Double(t)/(Double(P)*Double(N-r)))};let raw=summarize(Array(UnsafeBufferPointer(start:moment.contents().bindMemory(to:UInt64.self,capacity:G),count:G))),z=exp(logZ),u=z*Double(raw)/Double(ts[N-1])/Double(N);return Rep(logZ:logZ,z:z,u:u,minAnc:-1,dead:false,gpuSeconds:cb.gpuEndTime-cb.gpuStartTime,dispatchSeconds:wall)
}

/* Original r013 CLI retained here for provenance; this audit has its own entry point.
let args=CommandLine.arguments
guard args.count==7 || (args.count==8 && args[7]=="full") || (args.count==9 && args[7]=="full" && args[8]=="strat") || (args.count==10 && args[7]=="full" && args[8]=="strat" && UInt64(args[9],radix:16) != nil) else { fatalError("usage: tau_smc_integer n particles replicates tauLo tauHi output.jsonl [full [strat SEEDDOMAINHEX]]") }
let n=Int(args[1])!, P=Int(args[2])!, R=Int(args[3])!, lo=Int(args[4])!, hi=Int(args[5])!
let fullGPU=args.count>=8, stratified=args.count>=9, seedDomain=args.count==10 ? UInt64(args[9],radix:16)! : UInt64(0)
precondition(n >= 2 && n <= 12 && P > 0 && R > 0 && lo >= 2 && hi <= n*n && lo <= hi, "invalid n/P/R/tau range")
precondition(UInt64(P) <= UInt64.max / UInt64(P) / UInt64(n*n), "P*T may overflow UInt64")
let out=URL(fileURLWithPath:args[6]); FileManager.default.createFile(atPath:out.path, contents:nil); let fh=try FileHandle(forWritingTo:out); let wallStart=Date()
func emit(_ s:String){try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!)}
emit("{\"kind\":\"metadata\",\"device\":\"\(device.name)\",\"n\":\(n),\"particles\":\(P),\"replicates\":\(R),\"tau_lo\":\(lo),\"tau_hi\":\(hi),\"seed_domain_hex\":\"\(String(seedDomain,radix:16))\",\"method\":\(stratified ? "\"single-command-buffer integer-stratified\"" : (fullGPU ? "\"single-command-buffer GPU CDF multinomial\"" : "\"integer group resampling; GPU parent CDF search\"")),\"timing\":\"gpu_hardware_seconds=sum Metal GPU timestamps; synchronous_dispatch_seconds=submit+wait on host\",\"distinct_parent_metric\":null}")
for tau in lo...hi { for rep in 0..<R { let seed=(0xA022000000000000 &+ seedDomain) &+ UInt64(n*1000000+tau*10000+rep); let q=fullGPU ? runRepFull(n:n,tau:tau,P:P,seed:seed,stratified:stratified) : runRep(n:n,tau:tau,P:P,seed:seed); let logZJSON = q.dead ? "null" : "\(q.logZ)"; emit("{\"kind\":\"replicate\",\"tau\":\(tau),\"replicate\":\(rep),\"seed_hex\":\"\(String(seed,radix:16))\",\"logZ\":\(logZJSON),\"Z\":\(q.z),\"U\":\(q.u),\"min_distinct_parents\":null,\"dead\":\(q.dead),\"gpu_hardware_seconds\":\(q.gpuSeconds),\"synchronous_dispatch_seconds\":\(q.dispatchSeconds)}") } }
emit("{\"kind\":\"run_end\",\"wall_seconds\":\(-wallStart.timeIntervalSinceNow)}")
try fh.close()
*/

// MARK: - Independent CPU reference and rank-by-rank audit
func stateEqual(_ a: State, _ b: State) -> Bool { a.m0 == b.m0 && a.m1 == b.m1 && a.m2 == b.m2 && a.m3 == b.m3 && a.m4 == b.m4 && a.root1 == b.root1 && a.root2 == b.root2 }
func hasCPU(_ s: State, _ v: Int) -> Bool { let q=v >> 5, b=UInt32(v & 31); let w:UInt32; switch q { case 0:w=s.m0; case 1:w=s.m1; case 2:w=s.m2; case 3:w=s.m3; default:w=s.m4 }; return (w & (UInt32(1) << b)) != 0 }
func addCPU(_ s: inout State, _ v: Int) { let q=v >> 5, b=UInt32(v & 31), bit=UInt32(1) << b; switch q { case 0:s.m0 |= bit; case 1:s.m1 |= bit; case 2:s.m2 |= bit; case 3:s.m3 |= bit; default:s.m4 |= bit } }
func touchCPU(_ s: State, _ v: Int, _ n: Int) -> Bool { let x=v % n,y=v / n; return (x > 0 && hasCPU(s,v-1)) || (x+1 < n && hasCPU(s,v+1)) || (y > 0 && hasCPU(s,v-n)) || (y+1 < n && hasCPU(s,v+n)) }
func allowedCPU(_ s: State, _ v: Int, _ n: Int, _ tau: Int, _ rank: Int) -> Bool { if hasCPU(s,v) { return false }; if rank == 1 { return true }; let t=touchCPU(s,v,n); return rank == tau ? !t : t }

struct RankAudit { let counts:[UInt32]; let groups:[UInt64]; let cdf:[UInt64]; let total:UInt64; let parents:[UInt32]; let states:[State]; let dead:Bool }
struct AuditRun { let ranks:[RankAudit]; let terminalTotal:UInt64; let z:Double; let u:Double; let dead:Bool }

func cpuReference(n:Int, tau:Int, P:Int, seed:UInt64, stratified:Bool, forceDeadAt:Int?=nil) -> AuditRun {
    let N=n*n, G=groups(P); var states=[State](repeating:State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0),count:P), dead=false, ranks=[RankAudit](), totals=[UInt64]()
    for rank in 1...N {
        if forceDeadAt == rank { dead=true }
        let counts:[UInt32] = dead ? [UInt32](repeating:0,count:P) : states.map { s in UInt32((0..<N).reduce(0) { $0 + (allowedCPU(s,$1,n,tau,rank) ? 1 : 0) }) }
        var gs=[UInt64](repeating:0,count:G); for i in 0..<P { gs[i/256] &+= UInt64(counts[i]) }
        var cdf=[UInt64](repeating:0,count:G), total:UInt64=0; for i in 0..<G { total &+= gs[i]; cdf[i]=total }
        if total == 0 { dead=true }
        var parents=[UInt32](repeating:UInt32.max,count:P)
        if rank < N && !dead {
            var dst=[State](repeating:states[0],count:P)
            for child in 0..<P {
                var r0=childKey(seed,rank,child,1); let t=below(&r0,total); let u=stratified ? (UInt64(child)*total+t)/UInt64(P) : t
                var lo=0,hi=G-1; while lo < hi { let mid=(lo+hi)/2; if u < cdf[mid] { hi=mid } else { lo=mid+1 } }
                let start=lo*256,end=min(start+256,P), previous=(lo == 0 ? 0 : cdf[lo-1]); var c:UInt64=0, parent=start
                for k in start..<end { if u-previous < c + UInt64(counts[k]) { parent=k; break }; c &+= UInt64(counts[k]) }
                parents[child]=UInt32(parent); var s=states[parent], rng=childKey(seed,rank,child,2); let pick=below(&rng,UInt64(counts[parent])); var seen:UInt64=0, v=0
                while v < N { if allowedCPU(s,v,n,tau,rank) { if seen == pick { break }; seen += 1 }; v += 1 }
                if rank == 1 { s.root1=UInt32(v) }; if rank == tau { s.root2=UInt32(v) }; addCPU(&s,v); dst[child]=s
            }
            states=dst
        }
        totals.append(total); ranks.append(RankAudit(counts:counts,groups:gs,cdf:cdf,total:total,parents:parents,states:states,dead:dead))
    }
    if totals.contains(0) { return AuditRun(ranks:ranks,terminalTotal:totals.last!,z:0,u:0,dead:true) }
    var logZ=0.0; for (r,t) in totals.enumerated() { logZ += log(Double(t)/(Double(P)*Double(N-r))) }
    let last=ranks[N-1], raw:UInt64 = (0..<P).reduce(0) { acc,i in var r2=last.states[i].root2; if tau == N { for v in 0..<N where allowedCPU(last.states[i],v,n,tau,N) { r2=UInt32(v); break } }; let dx=Int(last.states[i].root1 % UInt32(n))-Int(r2 % UInt32(n)),dy=Int(last.states[i].root1 / UInt32(n))-Int(r2 / UInt32(n)); return acc &+ UInt64(last.counts[i])*UInt64(dx*dx+dy*dy) }
    let z=exp(logZ); return AuditRun(ranks:ranks,terminalTotal:totals.last!,z:z,u:z*Double(raw)/Double(totals.last!)/Double(N),dead:false)
}

func read32(_ b:MTLBuffer,_ n:Int)->[UInt32] { Array(UnsafeBufferPointer(start:b.contents().bindMemory(to:UInt32.self,capacity:n),count:n)) }
func read64(_ b:MTLBuffer,_ n:Int)->[UInt64] { Array(UnsafeBufferPointer(start:b.contents().bindMemory(to:UInt64.self,capacity:n),count:n)) }
func readStates(_ b:MTLBuffer,_ n:Int)->[State] { Array(UnsafeBufferPointer(start:b.contents().bindMemory(to:State.self,capacity:n),count:n)) }
func encodeU64Offset(_ cb:MTLCommandBuffer,_ src:MTLBuffer,_ offset:Int,_ dst:MTLBuffer) { let e=cb.makeComputeCommandEncoder()!; e.setComputePipelineState(copyU64Pipe); e.setBuffer(src,offset:offset,index:0); e.setBuffer(dst,offset:0,index:1); e.dispatchThreads(MTLSize(width:1,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:1,height:1,depth:1)); e.endEncoding() }

// The same instrumented respropFull is run either as one production-style
// command buffer or one command buffer per rank. Observers are copy-only.
func gpuAudit(n:Int, tau:Int, P:Int, seed:UInt64, stratified:Bool, oneCommand:Bool, forceDeadAt:Int?=nil)->AuditRun {
    let N=n*n,G=groups(P),stateBytes=P*MemoryLayout<State>.stride
    var pc=UInt32(P),gg=UInt32(G),dead0:UInt32=0,smode:UInt32=stratified ? 1:0
    let pcBuf=device.makeBuffer(bytes:&pc,length:4,options:.storageModeShared)!,gBuf=device.makeBuffer(bytes:&gg,length:4,options:.storageModeShared)!,deadBuf=device.makeBuffer(bytes:&dead0,length:4,options:.storageModeShared)!,smodeBuf=device.makeBuffer(bytes:&smode,length:4,options:.storageModeShared)!
    let countBuf=device.makeBuffer(length:P*4,options:.storageModeShared)!,groupBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,cdfBuf=device.makeBuffer(length:G*8,options:.storageModeShared)!,totals=device.makeBuffer(length:N*8,options:.storageModeShared)!,moment=device.makeBuffer(length:G*8,options:.storageModeShared)!
    let x=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,y=device.makeBuffer(length:stateBytes,options:.storageModeShared)!; x.contents().bindMemory(to:State.self,capacity:P).initialize(repeating:State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0),count:P)
    var pars=[MTLBuffer](),ris=[MTLBuffer](); for rank in 1...N { var p=Par(n:UInt32(n),N:UInt32(N),tau:UInt32(tau),rank:UInt32(rank),seedlo:UInt32(truncatingIfNeeded:seed),seedhi:UInt32(truncatingIfNeeded:seed>>32)),ri=UInt32(rank-1); pars.append(device.makeBuffer(bytes:&p,length:MemoryLayout<Par>.stride,options:.storageModeShared)!); ris.append(device.makeBuffer(bytes:&ri,length:4,options:.storageModeShared)!) }
    let cObs=(0..<N).map{_ in device.makeBuffer(length:P*4,options:.storageModeShared)!},gObs=(0..<N).map{_ in device.makeBuffer(length:G*8,options:.storageModeShared)!},cdfObs=(0..<N).map{_ in device.makeBuffer(length:G*8,options:.storageModeShared)!},totalObs=(0..<N).map{_ in device.makeBuffer(length:8,options:.storageModeShared)!},deadObs=(0..<N).map{_ in device.makeBuffer(length:4,options:.storageModeShared)!},parentObs=(0..<N).map{_ in device.makeBuffer(length:P*4,options:.storageModeShared)!},stateObs=(0..<N).map{_ in device.makeBuffer(length:stateBytes,options:.storageModeShared)!}
    for b in parentObs { b.contents().bindMemory(to:UInt32.self,capacity:P).initialize(repeating:UInt32.max,count:P) }
    var src=x,dst=y; let single=oneCommand ? queue.makeCommandBuffer()! : nil
    for r in 0..<N {
        let cb=single ?? queue.makeCommandBuffer()!
        if forceDeadAt == r+1 { encode(cb,setDeadPipe,[(deadBuf,0)],1) }
        encode(cb,fullCountPipe,[(src,0),(countBuf,1),(pars[r],2),(deadBuf,3)],P); encode(cb,groupPipe,[(countBuf,0),(groupBuf,1),(pcBuf,2)],G); encode(cb,fullCDFPipe,[(groupBuf,0),(cdfBuf,1),(totals,2),(deadBuf,3),(gBuf,4),(ris[r],5)],1)
        encode(cb,copyU32Pipe,[(countBuf,0),(cObs[r],1)],P); encode(cb,copyU64Pipe,[(groupBuf,0),(gObs[r],1)],G); encode(cb,copyU64Pipe,[(cdfBuf,0),(cdfObs[r],1)],G); encodeU64Offset(cb,totals,r*8,totalObs[r]); encode(cb,copyU32Pipe,[(deadBuf,0),(deadObs[r],1)],1)
        if r+1 == N { encode(cb,copyStatePipe,[(src,0),(stateObs[r],1)],P); encode(cb,fullTermPipe,[(src,0),(countBuf,1),(moment,2),(pars[r],3),(pcBuf,4),(deadBuf,5)],G) } else { encode(cb,fullPropPipe,[(src,0),(countBuf,1),(cdfBuf,2),(totals,3),(dst,4),(pars[r],5),(pcBuf,6),(ris[r],7),(deadBuf,8),(smodeBuf,9),(parentObs[r],10)],P); encode(cb,copyStatePipe,[(dst,0),(stateObs[r],1)],P); swap(&src,&dst) }
        if !oneCommand { cb.commit(); cb.waitUntilCompleted(); precondition(cb.status == .completed,cb.error?.localizedDescription ?? "rank GPU failure") }
    }
    if let cb=single { cb.commit(); cb.waitUntilCompleted(); precondition(cb.status == .completed,cb.error?.localizedDescription ?? "single GPU failure") }
    let ts=read64(totals,N), isDead=ts.contains(0), ranks=(0..<N).map { r in RankAudit(counts:read32(cObs[r],P),groups:read64(gObs[r],G),cdf:read64(cdfObs[r],G),total:read64(totalObs[r],1)[0],parents:read32(parentObs[r],P),states:readStates(stateObs[r],P),dead:read32(deadObs[r],1)[0] != 0) }
    if isDead { return AuditRun(ranks:ranks,terminalTotal:ts.last!,z:0,u:0,dead:true) }
    var logZ=0.0; for (r,t) in ts.enumerated(){logZ += log(Double(t)/(Double(P)*Double(N-r)))}; let raw=read64(moment,G).reduce(0,&+),z=exp(logZ); return AuditRun(ranks:ranks,terminalTotal:ts.last!,z:z,u:z*Double(raw)/Double(ts.last!)/Double(N),dead:false)
}

func firstDiff<T:Equatable>(_ a:[T],_ b:[T])->String? { if a.count != b.count { return "length \(a.count) != \(b.count)" }; for i in a.indices where a[i] != b[i] { return "index \(i): \(a[i]) != \(b[i])" }; return nil }
func validate(_ expected:AuditRun,_ actual:AuditRun,_ label:String)->[String] { var out=[String](); if expected.dead != actual.dead { out.append("\(label) dead \(actual.dead) != \(expected.dead)") }; if expected.terminalTotal != actual.terminalTotal { out.append("\(label) terminal total \(actual.terminalTotal) != \(expected.terminalTotal)") }; if expected.z != actual.z || expected.u != actual.u { out.append("\(label) terminal Z/U \(actual.z),\(actual.u) != \(expected.z),\(expected.u)") }; for r in expected.ranks.indices { let e=expected.ranks[r],a=actual.ranks[r]; if let d=firstDiff(e.counts,a.counts){out.append("\(label) rank \(r+1) counts \(d)")}; if let d=firstDiff(e.groups,a.groups){out.append("\(label) rank \(r+1) groups \(d)")}; if let d=firstDiff(e.cdf,a.cdf){out.append("\(label) rank \(r+1) cdf \(d)")}; if e.total != a.total {out.append("\(label) rank \(r+1) total \(a.total) != \(e.total)")}; if let d=firstDiff(e.parents,a.parents){out.append("\(label) rank \(r+1) parents \(d)")}; for i in 0..<e.states.count where !stateEqual(e.states[i],a.states[i]) {out.append("\(label) rank \(r+1) state index \(i) differs");break}; if e.dead != a.dead {out.append("\(label) rank \(r+1) dead \(a.dead) != \(e.dead)")} }; return out }

struct AuditCase { let n:Int; let tau:Int; let p:Int; let forceDeadAt:Int? }
var cases=[AuditCase](); for t in 2...4 { cases.append(AuditCase(n:2,tau:t,p:513,forceDeadAt:nil)) }; for t in 2...9 { cases.append(AuditCase(n:3,tau:t,p:513,forceDeadAt:nil)) }; for t in [2,8,47,48] { cases.append(AuditCase(n:7,tau:t,p:17,forceDeadAt:nil)) }; cases.append(AuditCase(n:3,tau:5,p:513,forceDeadAt:3))
var failures=[String](), records=[String](), completed=0; let started=Date()
print("GPU full audit device: \(device.name); cases=\(cases.count), modes=2")
for c in cases { for stratified in [false,true] { let seed=0xA022000000000000 &+ UInt64(c.n*1_000_000+c.tau*10_000+(stratified ? 1 : 0)); let cpu=cpuReference(n:c.n,tau:c.tau,P:c.p,seed:seed,stratified:stratified,forceDeadAt:c.forceDeadAt),single=gpuAudit(n:c.n,tau:c.tau,P:c.p,seed:seed,stratified:stratified,oneCommand:true,forceDeadAt:c.forceDeadAt),step=gpuAudit(n:c.n,tau:c.tau,P:c.p,seed:seed,stratified:stratified,oneCommand:false,forceDeadAt:c.forceDeadAt); let tag="n\(c.n) tau\(c.tau) P\(c.p) \(stratified ? "strat" : "multi")\(c.forceDeadAt == nil ? "" : " forced-dead-r\(c.forceDeadAt!)")"; failures += validate(cpu,single,"single \(tag)"); failures += validate(cpu,step,"per-rank \(tag)"); failures += validate(single,step,"single-vs-per-rank \(tag)"); if c.forceDeadAt == nil { let production=runRepFull(n:c.n,tau:c.tau,P:c.p,seed:seed,stratified:stratified); if production.dead != single.dead || production.z != single.z || production.u != single.u { failures.append("production-single terminal mismatch \(tag)") } }; completed += 1; records.append("{\"case\":\"\(tag)\",\"seed_hex\":\"\(String(seed,radix:16))\",\"terminal_total\":\(single.terminalTotal),\"dead\":\(single.dead),\"Z\":\(single.z),\"U\":\(single.u)}"); print("PASS \(tag) terminalTotal=\(single.terminalTotal) dead=\(single.dead) Z=\(single.z) U=\(single.u)") } }
print("AUDIT_RESULT cases=\(completed) comparisons=\(completed*3) failures=\(failures.count) seconds=\(-started.timeIntervalSinceNow)")
let artifact="{\"device\":\"\(device.name)\",\"cases\":\(completed),\"cpu_single_stepwise_comparisons\":\(completed*3),\"production_terminal_comparisons\":\((cases.count-1)*2),\"failures\":\(failures.count),\"records\":[\(records.joined(separator:","))]}\n"
try! artifact.write(toFile:"/Volumes/data/CodexTemp/T2/a022-gpu-full-audit-final/audit-results.json",atomically:true,encoding:.utf8)
if !failures.isEmpty { for f in failures.prefix(100) { fputs("FAIL \(f)\n",stderr) }; exit(1) }
