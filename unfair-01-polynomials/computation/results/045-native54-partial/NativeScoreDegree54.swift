import Foundation
import Metal

struct Input { let groups:Int; let nfreq:Int; let ncoef:Int; let lr:[Float]; let li:[Float]; let dc:[Float]; let ds:[Float] }
func flatten(_ x: Any) -> [Float] { (x as! [[Double]]).flatMap { $0.map(Float.init) } }
func load(_ p:String) throws -> Input {
  let d=try JSONSerialization.jsonObject(with:Data(contentsOf:URL(fileURLWithPath:p))) as! [String:Any]
  let lr=flatten(d["LR"]!), li=flatten(d["LI"]!), dc=flatten(d["DC"]!), ds=flatten(d["DS"]!)
  let g=d["groups"] as! Int; let nr=(d["DC"] as! [[Double]]).count; let nc=(d["DC"] as! [[Double]])[0].count
  guard g > 1 && g < 32 && nr <= 64 && nc <= 64 && lr.count == g*nr && li.count == g*nr && dc.count == nr*nc && ds.count == nr*nc && lr.allSatisfy({$0.isFinite}) && li.allSatisfy({$0.isFinite}) && dc.allSatisfy({$0.isFinite}) && ds.allSatisfy({$0.isFinite}) else { throw NSError(domain:"input",code:1) }
  return Input(groups:g,nfreq:nr,ncoef:nc,lr:lr,li:li,dc:dc,ds:ds)
}
let msl="""
#include <metal_stdlib>
using namespace metal;
kernel void scoreMasks(device const float* lr [[buffer(0)]], device const float* li [[buffer(1)]], device const float* dc [[buffer(2)]], device const float* ds [[buffer(3)]], device const uint* masks [[buffer(4)]], device float* scores [[buffer(5)]], device uint* finite [[buffer(6)]], constant uint& groups [[buffer(7)]], constant uint& nr [[buffer(8)]], constant uint& nc [[buffer(9)]], uint id [[thread_position_in_grid]]) {
 float ar[64],ai[64],br[64],bi[64]; uint mask=masks[id];
 float avr[64],avi[64],bvr[64],bvi[64]; uint ok=1;
 for(uint r=0;r<nr;r++){float x=0,y=0,u=0,v=0; for(uint g=0;g<groups;g++){uint bit=(g==0)?0u:((mask>>(groups-1-g))&1u);float R=lr[g*nr+r],I=li[g*nr+r]; x+=bit*R;y+=bit*I;u+=(1-bit)*R;v+=(1-bit)*I;} ar[r]=x;ai[r]=y;br[r]=u;bi[r]=v;avr[r]=exp(x)*cos(y);avi[r]=exp(x)*sin(y);bvr[r]=exp(u)*cos(v);bvi[r]=exp(u)*sin(v);if(!isfinite(x)||!isfinite(y)||!isfinite(u)||!isfinite(v)||!isfinite(avr[r])||!isfinite(avi[r])||!isfinite(bvr[r])||!isfinite(bvi[r])){ok=0;break;}}
 float score=0;
 float amin=INFINITY,bmin=INFINITY,amax=-INFINITY,bmax=-INFINITY,ima=0,imb=0;
 if(ok){for(uint j=0;j<nc;j++){float A=0,B=0,IA=0,IB=0;for(uint r=0;r<nr;r++){float C=dc[r*nc+j],S=ds[r*nc+j];A+=avr[r]*C+avi[r]*S;B+=bvr[r]*C+bvi[r]*S;IA+=avi[r]*C-avr[r]*S;IB+=bvi[r]*C-bvr[r]*S;}if(!isfinite(A)||!isfinite(B)||!isfinite(IA)||!isfinite(IB)){ok=0;break;}amin=min(amin,A);amax=max(amax,A);bmin=min(bmin,B);bmax=max(bmax,B);ima=max(ima,abs(IA));imb=max(imb,abs(IB));}} float q=ok?max(max(max(-amin,-bmin),max(amax-1,bmax-1)),max(ima,imb)):INFINITY;if(!isfinite(q)){ok=0;q=INFINITY;}score=q; if(mask==0) score=INFINITY; scores[id]=score;finite[id]=ok;
}
kernel void topK(device const float* scores [[buffer(0)]], device const uint* finite [[buffer(1)]], device uint* outIndex [[buffer(2)]], device float* outScore [[buffer(3)]], constant uint& n [[buffer(4)]], uint tid [[thread_position_in_grid]]) { if(tid) return; float best[4]={INFINITY,INFINITY,INFINITY,INFINITY};uint ix[4]={0,0,0,0};for(uint i=0;i<n;i++){float x=finite[i]?scores[i]:INFINITY;for(uint k=0;k<4;k++){if(x<best[k]){for(uint q=3;q>k;q--){best[q]=best[q-1];ix[q]=ix[q-1];}best[k]=x;ix[k]=i;break;}}}for(uint k=0;k<4;k++){outIndex[k]=ix[k];outScore[k]=best[k];}}
"""
func buf<T>(_ d:MTLDevice,_ a:[T],_ opt:MTLResourceOptions = .storageModeShared)->MTLBuffer { d.makeBuffer(bytes:a,length:MemoryLayout<T>.stride*a.count,options:opt)! }
let args=CommandLine.arguments; guard let ip=args.firstIndex(of:"--input"),let op=args.firstIndex(of:"--out") else { fatalError("args") }
let input=try load(args[ip+1]); guard let dev=MTLCreateSystemDefaultDevice(),let queue=dev.makeCommandQueue() else { fatalError("no Metal") }
let lib=try dev.makeLibrary(source:msl,options:nil);let score=try dev.makeComputePipelineState(function:lib.makeFunction(name:"scoreMasks")!);let top=try dev.makeComputePipelineState(function:lib.makeFunction(name:"topK")!)
let lr=buf(dev,input.lr),li=buf(dev,input.li),dc=buf(dev,input.dc),ds=buf(dev,input.ds); var g=UInt32(input.groups),nr=UInt32(input.nfreq),nc=UInt32(input.ncoef)
func oneBatch(_ masks:[UInt32]) -> [String:Any] {
  guard !masks.isEmpty && masks.allSatisfy({ $0 < (UInt32(1) << UInt32(input.groups-1)) }) else { fatalError("mask range") }
  let mb=buf(dev,masks),sb=dev.makeBuffer(length:4*masks.count,options:.storageModeShared)!,fb=dev.makeBuffer(length:4*masks.count,options:.storageModeShared)!,ib=dev.makeBuffer(length:16,options:.storageModeShared)!,tb=dev.makeBuffer(length:16,options:.storageModeShared)!;var n=UInt32(masks.count)
  let cb=queue.makeCommandBuffer()!;let e=cb.makeComputeCommandEncoder()!;e.setComputePipelineState(score);[lr,li,dc,ds,mb,sb,fb].enumerated().forEach{e.setBuffer($0.element,offset:0,index:$0.offset)};e.setBytes(&g,length:4,index:7);e.setBytes(&nr,length:4,index:8);e.setBytes(&nc,length:4,index:9);e.dispatchThreads(MTLSize(width:masks.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(256,masks.count),height:1,depth:1));e.endEncoding();let r=cb.makeComputeCommandEncoder()!;r.setComputePipelineState(top);[sb,fb,ib,tb].enumerated().forEach{r.setBuffer($0.element,offset:0,index:$0.offset)};r.setBytes(&n,length:4,index:4);r.dispatchThreads(MTLSize(width:1,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:1,height:1,depth:1));r.endEncoding();cb.commit();cb.waitUntilCompleted();guard cb.status == .completed else { fatalError("metal failure") }
  let scores=Array(UnsafeBufferPointer(start:sb.contents().assumingMemoryBound(to:Float.self),count:masks.count));let finite=Array(UnsafeBufferPointer(start:fb.contents().assumingMemoryBound(to:UInt32.self),count:masks.count));let ix=Array(UnsafeBufferPointer(start:ib.contents().assumingMemoryBound(to:UInt32.self),count:4));let ts=Array(UnsafeBufferPointer(start:tb.contents().assumingMemoryBound(to:Float.self),count:4));let jsonScores:[Any]=scores.count<=64 ? scores.map { x in if x.isFinite { return Double(x) as Any }; return NSNull() } : [];let jsonTop:[Any]=ts.map { x in if x.isFinite { return Double(x) as Any }; return NSNull() }
  return ["metal_device":dev.name,"groups":input.groups,"nfreq":input.nfreq,"ncoef":input.ncoef,"masks":masks.count<=64 ? masks.map(Int.init) : [],"range_count":masks.count,"scores":jsonScores,"finite":masks.count<=64 ? finite : [],"finite_count":finite.filter{$0==1}.count,"top_indices":ix.map(Int.init),"top_scores":jsonTop,"command_status":"completed"]
}
func arg(_ key:String) -> String? { guard let i=args.firstIndex(of:key), i+1<args.count else{return nil}; return args[i+1] }
var ranges:[(UInt32,UInt32)] = []
if let mp=arg("--masks") { let o=oneBatch(mp.split(separator:",").map{UInt32($0)!}); try JSONSerialization.data(withJSONObject:o,options:[.sortedKeys]).write(to:URL(fileURLWithPath:args[op+1])); exit(0) }
if let rp=arg("--range") { let z=rp.split(separator:",").map{UInt32($0)!}; guard z.count==2 else { fatalError("range") }; ranges=[(z[0],z[1])] }
else if let rp=arg("--ranges") { ranges=rp.split(separator:";").map { token in let z=token.split(separator:",").map{UInt32($0)!}; guard z.count==2 else { fatalError("ranges") }; return (z[0],z[1]) } }
else { fatalError("range") }
let journal=arg("--journal"), task=arg("--task"), inputSHA=arg("--input-sha"), sourceSHA=arg("--source-sha"), binarySHA=arg("--binary-sha")
if journal != nil && (task == nil || inputSHA == nil || sourceSHA == nil || binarySHA == nil) { fatalError("journal binding") }
let maxSeconds = arg("--max-seconds").flatMap(Double.init)
if maxSeconds != nil && (maxSeconds! <= 0 || !maxSeconds!.isFinite) { fatalError("max-seconds") }
if let j=journal { FileManager.default.createFile(atPath:j,contents:nil); }
let handle = journal == nil ? nil : try FileHandle(forWritingTo:URL(fileURLWithPath:journal!))
func emit(_ event:[String:Any]) throws { guard let h=handle else{return}; var row=event; row["event"]="batch"; let data=try JSONSerialization.data(withJSONObject:row,options:[.sortedKeys]); h.write(data); h.write(Data([10])); try h.synchronize() }
var batches:[[String:Any]]=[]
let started=Date(); var runStatus="completed"
for (start,count) in ranges { if let cap=maxSeconds, Date().timeIntervalSince(started) >= cap { runStatus="partial_time_cap"; break }; guard count>0 && start <= UInt32.max-count else { fatalError("range overflow") }; var b=oneBatch(Array(start..<(start+count))); b["range_start"]=Int(start); b["task"]=task; b["input_sha256"]=inputSHA; b["source_sha256"]=sourceSHA; b["binary_sha256"]=binarySHA; b["elapsed_seconds"]=Date().timeIntervalSince(started); try emit(b); batches.append(b) }
try handle?.close()
let out:[String:Any]=["command_status":runStatus,"batches":batches,"metal_device":dev.name,"groups":input.groups,"nfreq":input.nfreq,"ncoef":input.ncoef,"task":task as Any,"input_sha256":inputSHA as Any,"source_sha256":sourceSHA as Any,"binary_sha256":binarySHA as Any]
try JSONSerialization.data(withJSONObject:out,options:[.sortedKeys]).write(to:URL(fileURLWithPath:args[op+1]))
