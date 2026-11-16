import Foundation
import Metal

// Fixed-tau original-law SIS.  A failed forced prefix emits alive=0 and is
// retained as a zero contribution in the denominator.
let msl = """
#include <metal_stdlib>
using namespace metal;
struct P { uint n,N,count,tau,seed; };
struct O { float logw; uint alive,e2,rho2,first,second; };
inline uint mix(uint x) { x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; x^=x>>16; return x|1u; }
inline uint rnd(thread uint& x) { x^=x<<13; x^=x>>17; x^=x<<5; return x; }
inline bool has(const thread uint* a,uint v) { return (a[v>>5]&(1u<<(v&31)))!=0; }
inline void add(thread uint* a,uint v) { a[v>>5]|=1u<<(v&31); }
inline bool touch(const thread uint* a,uint v,uint n) { uint x=v%n,y=v/n; return (x&&has(a,v-1))||(x+1<n&&has(a,v+1))||(y&&has(a,v-n))||(y+1<n&&has(a,v+n)); }
kernel void sample(constant P&p[[buffer(0)]],device O*out[[buffer(1)]],uint id[[thread_position_in_grid]]) {
 if(id>=p.count)return; thread uint a[5]={0,0,0,0,0};thread uint cands[144];uint q=mix(p.seed^(id*0x9e3779b9u));uint one=0xffffffffu,two=0xffffffffu;float lw=0;uint alive=1;
 for(uint zero=0;zero<p.N;zero++){uint rank=zero+1,c=0;for(uint v=0;v<p.N;v++)if(!has(a,v)){bool t=touch(a,v,p.n);bool use=(rank==1)||(rank<p.tau&&t)||(rank==p.tau&&!t)||(rank>p.tau&&t);if(use)cands[c++]=v;}if(!c){alive=0;break;}uint v=cands[rnd(q)%c];lw+=log(float(c));if(rank==1)one=v;else if(rank==p.tau)two=v;add(a,v);}
 uint e2=0,rho2=0;if(alive){int dx=int(one%p.n)-int(two%p.n),dy=int(one/p.n)-int(two/p.n);uint ax=uint(abs(dx)),ay=uint(abs(dy));e2=ax*ax+ay*ay;rho2=(ax+ay)*(ax+ay);}out[id]=(O){lw,alive,e2,rho2,one,two};
}
"""

struct P { var n,N,count,tau,seed: UInt32 }
struct O { var logw: Float; var alive,e2,rho2,first,second: UInt32 }
struct Target { let n: Int; let a2: String?; let q: Double? }
func mix(_ x: UInt32) -> UInt32 { var y=x;y^=y>>16;y&*=0x7feb352d;y^=y>>15;y&*=0x846ca68b;y^=y>>16;return y|1 }
func rnd(_ x: inout UInt32)->UInt32{x^=x&<<13;x^=x>>17;x^=x&<<5;return x}
func has(_ a:[UInt32],_ v:Int)->Bool{(a[v>>5]&(1<<UInt32(v&31))) != 0}
func touch(_ a:[UInt32],_ v:Int,_ n:Int)->Bool{let x=v%n,y=v/n;return(x>0&&has(a,v-1))||(x+1<n&&has(a,v+1))||(y>0&&has(a,v-n))||(y+1<n&&has(a,v+n))}
func logadd(_ a:Double,_ b:Double)->Double{if a == -.infinity{return b};if b == -.infinity{return a};let m=max(a,b);return m+log(exp(a-m)+exp(b-m))}
func jsonNumber(_ x:Double)->String{x.isFinite ? String(x):"null"}
func trace(n:Int,tau:Int,seed:UInt32)->O{let N=n*n;var a=[UInt32](repeating:0,count:5),q=mix(seed),lw=0.0,one=UInt32.max,two=UInt32.max,alive:UInt32=1;for z in 0..<N{let rank=z+1;var c=[Int]();for v in 0..<N where !has(a,v){let t=touch(a,v,n);if rank==1 || (rank<tau&&t) || (rank==tau && !t) || (rank>tau&&t){c.append(v)}};if c.isEmpty{alive=0;break};let v=c[Int(rnd(&q)%UInt32(c.count))];lw+=log(Double(c.count));if rank==1{one=UInt32(v)}else if rank==tau{two=UInt32(v)};a[v>>5]|=1<<UInt32(v&31)};var e2:UInt32=0,rho2:UInt32=0;if alive != 0{let dx=abs(Int(one%UInt32(n))-Int(two%UInt32(n))),dy=abs(Int(one/UInt32(n))-Int(two/UInt32(n)));e2=UInt32(dx*dx+dy*dy);rho2=UInt32((dx+dy)*(dx+dy))};return O(logw:Float(lw),alive:alive,e2:e2,rho2:rho2,first:one,second:two)}

let args=Array(CommandLine.arguments.dropFirst());guard args.count == 2 else { fatalError("usage: gpu_original_law_tau_adaptive OUTDIR ALLOCATION.json") };let dir=URL(fileURLWithPath:args[0]);let allocationURL=URL(fileURLWithPath:args[1]);try FileManager.default.createDirectory(at:dir,withIntermediateDirectories:true);let path=dir.appendingPathComponent("results.jsonl");guard !FileManager.default.fileExists(atPath:path.path) else{fatalError("choose empty output directory")};FileManager.default.createFile(atPath:path.path,contents:nil);let fh=try FileHandle(forWritingTo:path);func emit(_ s:String){try!fh.write(contentsOf:(s+"\n").data(using:.utf8)!)}
let dev=MTLCreateSystemDefaultDevice()!,lib=try dev.makeLibrary(source:msl,options:nil),pipe=try dev.makeComputePipelineState(function:lib.makeFunction(name:"sample")!),queue=dev.makeCommandQueue()!;let batch=65536
let allocationObject=try JSONSerialization.jsonObject(with:Data(contentsOf:allocationURL)) as! [String:Any];let n=allocationObject["n"] as! Int;guard n == 8 else { fatalError("this production executable is n=8-only") };let allocation=allocationObject["allocation"] as! [[String:Any]];let N=n*n;guard allocation.count == N-2 && Set(allocation.map{$0["tau"] as! Int}) == Set(2..<N) else { fatalError("allocation must cover every tau exactly once") };let pb=dev.makeBuffer(length:MemoryLayout<P>.stride,options:.storageModeShared)!,ob=dev.makeBuffer(length:batch*MemoryLayout<O>.stride,options:.storageModeShared)!;var totalA = -Double.infinity,totalE = -Double.infinity,totalR = -Double.infinity,gpuAll=0.0,wallAll=0.0
emit("{\"kind\":\"metadata\",\"device\":\"\(dev.name)\",\"proposal\":\"fixed descending second-birth rank tau; dead prefixes retained at zero\",\"batch\":\(batch),\"rng\":\"mix32 then xorshift32 deterministic streams\",\"seed_domain\":\"0xA022E*** independent of pilot\",\"allocation_file\":\"\(allocationURL.path)\",\"pilot_not_combined\":true}")
for item in allocation.sorted(by: { ($0["tau"] as! Int) < ($1["tau"] as! Int) }) {let tau=item["tau"] as! Int;let batches=item["batches"] as! Int;var p=P(n:UInt32(n),N:UInt32(N),count:1,tau:UInt32(tau),seed:0xa022e101);pb.contents().bindMemory(to:P.self,capacity:1).pointee=p;let cpu=trace(n:n,tau:tau,seed:p.seed);let c0=queue.makeCommandBuffer()!,e0=c0.makeComputeCommandEncoder()!;e0.setComputePipelineState(pipe);e0.setBuffer(pb,offset:0,index:0);e0.setBuffer(ob,offset:0,index:1);e0.dispatchThreads(MTLSize(width:1,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:1,height:1,depth:1));e0.endEncoding();c0.commit();c0.waitUntilCompleted();guard c0.status == .completed else{fatalError("trace failed")};let gpu=ob.contents().bindMemory(to:O.self,capacity:1).pointee;precondition(cpu.alive==gpu.alive&&cpu.first==gpu.first&&cpu.second==gpu.second&&cpu.e2==gpu.e2&&cpu.rho2==gpu.rho2&&((cpu.alive==0)||abs(Double(cpu.logw)-Double(gpu.logw))<1e-4),"trace mismatch");emit("{\"kind\":\"trajectory_check\",\"n\":\(n),\"tau\":\(tau),\"ok\":true,\"alive\":\(gpu.alive)}")
  p.count=UInt32(batch);pb.contents().bindMemory(to:P.self,capacity:1).pointee=p;var la = -Double.infinity,le = -Double.infinity,lr = -Double.infinity,lw2 = -Double.infinity,mx = -Double.infinity,alive=0;var gt=0.0,wt=0.0
  for b in 0..<batches {p.seed=0xa022e101&+UInt32(n*1009+tau*131+b*7919);pb.contents().bindMemory(to:P.self,capacity:1).pointee=p;let cb=queue.makeCommandBuffer()!,en=cb.makeComputeCommandEncoder()!;en.setComputePipelineState(pipe);en.setBuffer(pb,offset:0,index:0);en.setBuffer(ob,offset:0,index:1);en.dispatchThreads(MTLSize(width:batch,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1));en.endEncoding();let st=Date();cb.commit();cb.waitUntilCompleted();let wall=Date().timeIntervalSince(st);guard cb.status == .completed else{fatalError("batch failed")};let g=cb.gpuEndTime-cb.gpuStartTime;gt+=g;wt+=wall;let o=ob.contents().bindMemory(to:O.self,capacity:batch);for i in 0..<batch where o[i].alive != 0{alive += 1;let x=Double(o[i].logw);la=logadd(la,x);lw2=logadd(lw2,2*x);le=logadd(le,x+log(Double(o[i].e2)));lr=logadd(lr,x+log(Double(o[i].rho2)));mx=max(mx,x)};emit("{\"kind\":\"batch\",\"n\":\(n),\"tau\":\(tau),\"index\":\(b),\"seed\":\(p.seed),\"gpu_seconds\":\(g),\"wall_seconds\":\(wall)}")}
  let M=Double(batch*batches),logA=la-log(M),logE=le-log(M),logR=lr-log(M);totalA=logadd(totalA,logA);totalE=logadd(totalE,logE);totalR=logadd(totalR,logR);gpuAll+=gt;wallAll+=wt;let ess=exp(2*la-lw2),top=exp(mx-la);emit("{\"kind\":\"stratum\",\"n\":\(n),\"tau\":\(tau),\"samples\":\(batch*batches),\"alive\":\(alive),\"log_A2_tau_hat\":\(jsonNumber(logA)),\"Q_tau_euclidean\":\(jsonNumber(exp(le-la)/Double(N))),\"ESS\":\(jsonNumber(ess)),\"top_weight\":\(jsonNumber(top)),\"gpu_seconds\":\(gt),\"wall_seconds\":\(wt)}")
 }
let q=exp(totalE-totalA)/Double(N),qr=exp(totalR-totalA)/Double(N);emit("{\"kind\":\"summary\",\"n\":\(n),\"log_A2_hat\":\(jsonNumber(totalA)),\"A2_ratio_to_exact\":null,\"Q_euclidean_hat\":\(jsonNumber(q)),\"Q_manhattan_hat\":\(jsonNumber(qr)),\"Q_error_to_exact\":null,\"gpu_seconds\":\(gpuAll),\"submit_wait_seconds\":\(wallAll)}")
try fh.close()
