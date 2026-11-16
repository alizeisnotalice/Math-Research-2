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
let oldPropPipe = try device.makeComputePipelineState(function: library.makeFunction(name: "resprop")!)
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

func sameBytes(_ a: MTLBuffer, _ b: MTLBuffer, _ n: Int) -> Bool { Data(bytes:a.contents(),count:n) == Data(bytes:b.contents(),count:n) }
func gpuAuditCase(n: Int, tau: Int, P: Int, seed: UInt64) -> [String:Any] {
    let N=n*n,G=groups(P),stateBytes=P*MemoryLayout<State>.stride,countBytes=P*MemoryLayout<UInt32>.stride
    var par=Par(n:UInt32(n),N:UInt32(N),tau:UInt32(tau),rank:1,seedlo:UInt32(truncatingIfNeeded:seed),seedhi:UInt32(truncatingIfNeeded:seed>>32))
    let parBuf=device.makeBuffer(bytes:&par,length:MemoryLayout<Par>.stride,options:.storageModeShared)!
    var pc=UInt32(P);let pcBuf=device.makeBuffer(bytes:&pc,length:4,options:.storageModeShared)!
    let oc=device.makeBuffer(length:countBytes,options:.storageModeShared)!,nc=device.makeBuffer(length:countBytes,options:.storageModeShared)!
    let gb=device.makeBuffer(length:G*MemoryLayout<UInt64>.stride,options:.storageModeShared)!,cdf=device.makeBuffer(length:G*MemoryLayout<UInt64>.stride,options:.storageModeShared)!
    var total:UInt64=0;let totalBuf=device.makeBuffer(bytes:&total,length:8,options:.storageModeShared)!
    let og=device.makeBuffer(length:P*MemoryLayout<UInt32>.stride,options:.storageModeShared)!,or=device.makeBuffer(length:P*MemoryLayout<UInt64>.stride,options:.storageModeShared)!
    let op=device.makeBuffer(length:P*MemoryLayout<UInt32>.stride,options:.storageModeShared)!,np=device.makeBuffer(length:P*MemoryLayout<UInt32>.stride,options:.storageModeShared)!
    let ox=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,oy=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,nx=device.makeBuffer(length:stateBytes,options:.storageModeShared)!,ny=device.makeBuffer(length:stateBytes,options:.storageModeShared)!
    for b in [ox,oy,nx,ny] { b.contents().bindMemory(to:State.self,capacity:P).initialize(repeating:State(m0:0,m1:0,m2:0,m3:0,m4:0,root1:0,root2:0),count:P) }
    var os=ox,od=oy,ns=nx,nd=ny;var checked=0
    for rank in 1...N {
        par.rank=UInt32(rank);parBuf.contents().bindMemory(to:Par.self,capacity:1).pointee=par
        _=dispatch(countPipe,[(os,0),(oc,1),(parBuf,2)],P);_=dispatch(countPipe,[(ns,0),(nc,1),(parBuf,2)],P)
        precondition(sameBytes(oc,nc,countBytes),"GPU count mismatch")
        let weights=Array(UnsafeBufferPointer(start:oc.contents().bindMemory(to:UInt32.self,capacity:P),count:P));total=weights.reduce(0){$0 &+ UInt64($1)}
        if total==0 { return ["n":n,"tau":tau,"particles":P,"ok":true,"dead_rank":rank,"checked_ranks":checked] }
        if rank==N {
            let om=device.makeBuffer(length:G*MemoryLayout<UInt64>.stride,options:.storageModeShared)!,nm=device.makeBuffer(length:G*MemoryLayout<UInt64>.stride,options:.storageModeShared)!
            _=dispatch(termPipe,[(os,0),(oc,1),(om,2),(parBuf,3),(pcBuf,4)],G);_=dispatch(termPipe,[(ns,0),(nc,1),(nm,2),(parBuf,3),(pcBuf,4)],G)
            precondition(sameBytes(om,nm,G*MemoryLayout<UInt64>.stride),"GPU terminal moment mismatch");checked += 1;break
        }
        var cum=[UInt64](),x:UInt64=0;for w in weights{x &+= UInt64(w);cum.append(x)}
        let ogp=og.contents().bindMemory(to:UInt32.self,capacity:P),orp=or.contents().bindMemory(to:UInt64.self,capacity:P)
        for child in 0..<P { var r=childKey(seed,rank,child,1);let u=below(&r,total);var lo=0,hi=G-1;while lo<hi{let m=(lo+hi)/2;if u<cum[m]{hi=m}else{lo=m+1}};let prev=lo==0 ? 0:cum[lo-1];ogp[child]=UInt32(lo);orp[child]=u-prev }
        cum.withUnsafeBufferPointer{cdf.contents().copyMemory(from:$0.baseAddress!,byteCount:G*MemoryLayout<UInt64>.stride)};totalBuf.contents().bindMemory(to:UInt64.self,capacity:1).pointee=total
        _=dispatch(oldPropPipe,[(os,0),(oc,1),(og,2),(or,3),(od,4),(op,5),(parBuf,6),(pcBuf,7)],P)
        _=dispatch(propPipe,[(ns,0),(nc,1),(cdf,2),(totalBuf,3),(nd,4),(np,5),(parBuf,6),(pcBuf,7)],P)
        precondition(sameBytes(op,np,P*MemoryLayout<UInt32>.stride),"GPU parent mismatch")
        precondition(sameBytes(od,nd,stateBytes),"GPU state mismatch")
        checked += 1;swap(&os,&od);swap(&ns,&nd)
    }
    return ["n":n,"tau":tau,"particles":P,"ok":true,"dead_rank":NSNull(),"checked_ranks":checked]
}

let args=CommandLine.arguments
if args.count==2 && args[1]=="audit" {
    var rows=[[String:Any]]()
    for tau in 2...4 { rows.append(gpuAuditCase(n:2,tau:tau,P:513,seed:0xA022000000000000 &+ UInt64(2_000_000+tau*10_000))) }
    for tau in 2...9 { rows.append(gpuAuditCase(n:3,tau:tau,P:513,seed:0xA022000000000000 &+ UInt64(3_000_000+tau*10_000))) }
    for tau in [2,8,47,48] { rows.append(gpuAuditCase(n:7,tau:tau,P:17,seed:0xA022000000000000 &+ UInt64(7_000_000+tau*10_000))) }
    let d=try! JSONSerialization.data(withJSONObject:["kind":"actual_gpu_state_parent_audit","cases":rows],options:[.prettyPrinted,.sortedKeys]);FileHandle.standardOutput.write(d);print("");exit(0)
}
guard args.count==7 else { fatalError("usage: tau_smc_integer n particles replicates tauLo tauHi output.jsonl") }
let n=Int(args[1])!, P=Int(args[2])!, R=Int(args[3])!, lo=Int(args[4])!, hi=Int(args[5])!
let out=URL(fileURLWithPath:args[6]); FileManager.default.createFile(atPath:out.path, contents:nil); let fh=try FileHandle(forWritingTo:out); let wallStart=Date()
func emit(_ s:String){try! fh.write(contentsOf:(s+"\n").data(using:.utf8)!)}
emit("{\"kind\":\"metadata\",\"device\":\"\(device.name)\",\"n\":\(n),\"particles\":\(P),\"replicates\":\(R),\"tau_lo\":\(lo),\"tau_hi\":\(hi),\"method\":\"integer group resampling; GPU parent CDF search; terminal no-resample\",\"timing\":\"gpu_hardware_seconds=sum Metal GPU timestamps; synchronous_dispatch_seconds=submit+wait on host\",\"distinct_parent_metric\":null}")
for tau in lo...hi { for rep in 0..<R { let seed=0xA022000000000000 &+ UInt64(n*1000000+tau*10000+rep); let q=runRep(n:n,tau:tau,P:P,seed:seed); let logZJSON = q.dead ? "null" : "\(q.logZ)"; emit("{\"kind\":\"replicate\",\"tau\":\(tau),\"replicate\":\(rep),\"logZ\":\(logZJSON),\"Z\":\(q.z),\"U\":\(q.u),\"min_distinct_parents\":null,\"dead\":\(q.dead),\"gpu_hardware_seconds\":\(q.gpuSeconds),\"synchronous_dispatch_seconds\":\(q.dispatchSeconds)}") } }
emit("{\"kind\":\"run_end\",\"wall_seconds\":\(-wallStart.timeIntervalSinceNow)}")
try fh.close()
