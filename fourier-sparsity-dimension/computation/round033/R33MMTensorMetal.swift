import Foundation
import Metal
import CryptoKit

let kernelSource = """
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]], device int *output [[buffer(1)]],
                       constant uint &n [[buffer(2)]], constant uint &stride [[buffer(3)]],
                       uint gid [[thread_position_in_grid]]) {
    uint local=gid % n, base=gid-local, mate=base+(local^stride);
    int left=input[gid], right=input[mate]; output[gid]=(local&stride)?(right-left):(left+right);
}
"""

func sha256(_ u: URL)->String { SHA256.hash(data:try!Data(contentsOf:u)).map{String(format:"%02x",$0)}.joined() }
func diagnostic(_ s: String) { FileHandle.standardError.write(Data((s + "\n").utf8)) }
func parity(_ x:Int)->Int{x.nonzeroBitCount&1}
func parseHexWords(_ s:String)->[UInt64]{
    var hex=String(s.dropFirst(2));precondition(hex.count<=64);if hex.count<64{hex=String(repeating:"0",count:64-hex.count)+hex};var out=[UInt64]();
    for i in stride(from:hex.count,through:16,by:-16){out.append(UInt64(hex[hex.index(hex.startIndex,offsetBy:i-16)..<hex.index(hex.startIndex,offsetBy:i)],radix:16)!)};return out
}
func bit(_ words:[UInt64],_ i:Int)->Int{Int((words[i>>6]>>UInt64(i&63))&1)}
func cpuFWHT(_ values:[Int32])->[Int32]{var a=values;var h=1;while h<a.count{for s in stride(from:0,to:a.count,by:2*h){for i in s..<(s+h){let u=a[i],v=a[i+h];a[i]=u+v;a[i+h]=u-v}};h <<= 1};return a}
func gpuFWHT(_ device:MTLDevice,_ pipeline:MTLComputePipelineState,_ input:[Int32],_ phases:Int,_ n:Int)throws->[Int32]{
 let bytes=input.count*MemoryLayout<Int32>.stride,op:MTLResourceOptions = .storageModeShared
 let first=input.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:bytes,options:op)!};let second=device.makeBuffer(length:bytes,options:op)!;var a=first,b=second
 let c=device.makeCommandQueue()!.makeCommandBuffer()!;var h=1;while h<n{let e=c.makeComputeCommandEncoder()!;e.setComputePipelineState(pipeline);e.setBuffer(a,offset:0,index:0);e.setBuffer(b,offset:0,index:1);var n32=UInt32(n),h32=UInt32(h);e.setBytes(&n32,length:4,index:2);e.setBytes(&h32,length:4,index:3);let w=min(256,pipeline.maxTotalThreadsPerThreadgroup);e.dispatchThreads(MTLSize(width:phases*n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:w,height:1,depth:1));e.endEncoding();swap(&a,&b);h <<= 1};c.commit();c.waitUntilCompleted();guard c.status == .completed else{throw NSError(domain:"metal",code:Int(c.status.rawValue))};return Array(UnsafeBufferPointer(start:a.contents().bindMemory(to:Int32.self,capacity:input.count),count:input.count))
}

let args=CommandLine.arguments;guard args.count==4 else{fatalError("usage: R33MMTensorMetal result.json signedmask.json output.json")}
let resultURL=URL(fileURLWithPath:args[1]), signedURL=URL(fileURLWithPath:args[2]),outURL=URL(fileURLWithPath:args[3])
let result=try!JSONSerialization.jsonObject(with:Data(contentsOf:resultURL))as![String:Any];let sol=result["solution"]as![String:Any]
let f=parseHexWords(sol["f_hex"]as!String),fp=parseHexWords(sol["f_plus_h_hex"]as!String);let h=zip(f,fp).map { $0 ^ $1 }
let signed=try!JSONSerialization.jsonObject(with:Data(contentsOf:signedURL))as![String:Any];let points=signed["points"]as![Int]
let hSupport=(0..<256).filter{bit(h,$0)==1};let supportMatches=hSupport==points && hSupport.count==32 && hSupport.allSatisfy{$0<128}
guard let device=MTLCreateSystemDefaultDevice()else{fatalError("Metal unavailable")}
let library: MTLLibrary; do { library=try device.makeLibrary(source:kernelSource,options:nil) } catch { diagnostic("Metal library error: \(error)"); exit(2) }
let pipeline: MTLComputePipelineState; do { pipeline=try device.makeComputePipelineState(function:library.makeFunction(name:"fwht_stage")!) } catch { diagnostic("Metal pipeline error: \(error)"); exit(3) }
diagnostic("R33 inputs and Metal pipeline ready; supportMatches=\(supportMatches)");let begun=Date();var rows=[[String:Any]]()
for r in 0...4 {
 let dimension=8+2*r,N=1<<dimension,q=1<<(dimension/2),uvMask=(1<<r)-1;var vals=[Int32](repeating:0,count:2*N);var diffWeight=0
 for t in 0..<2 {for point in 0..<N {let x=point&255,u=(point>>8)&uvMask,v=point>>(8+r);let base=bit(f,x);let flip=bit(h,x)*(v==0 ? 1:0);let value=base ^ parity(u&v) ^ (t==0 ? 0:flip);vals[t*N+point]=value==0 ? 1:-1;if t==1 && flip==1{diffWeight += 1}}}
 let cpu0=cpuFWHT(Array(vals[0..<N])),cpu1=cpuFWHT(Array(vals[N..<(2*N)]));diagnostic("R33 CPU FWHT complete r=\(r)");let gpu:[Int32];do{gpu=try gpuFWHT(device,pipeline,vals,2,N)}catch{diagnostic("Metal execution error r=\(r): \(error)");exit(4)};var mismatch=0,badCPU=0,badGPU=0
 for i in 0..<(2*N){let c=i<N ? cpu0[i]:cpu1[i-N];if c != gpu[i]{mismatch += 1};if Int(abs(c)) != q{badCPU += 1};if Int(abs(gpu[i])) != q{badGPU += 1}}
 rows.append(["r":r,"dimension":dimension,"N":N,"q":q,"difference_weight":diffWeight,"expected_difference_weight":32*(1<<r),"gpu_cpu_mismatch_count":mismatch,"cpu_bad_absq":badCPU,"gpu_bad_absq":badGPU,"all_checks_passed":mismatch==0 && badCPU==0 && badGPU==0 && diffWeight==32*(1<<r)])
}
let all=rows.allSatisfy{$0["all_checks_passed"]as!Bool};let output:[String:Any]=["schema":"a28-r033-mm-tensor-metal-v1","result_input_sha256":sha256(resultURL),"signed_mask_input_sha256":sha256(signedURL),"difference_support_matches_S_times_y8zero":supportMatches,"difference_support":hSupport,"device":device.name,"registry_id":String(device.registryID),"cpu_backend":"exact Int32 iterative FWHT","gpu_backend":"Metal Int32 global staged ping-pong FWHT, elementwise checked","rows":rows,"all_checks_passed":all && supportMatches,"elapsed_seconds":Date().timeIntervalSince(begun)]
try!JSONSerialization.data(withJSONObject:output,options:[.prettyPrinted,.sortedKeys]).write(to:outURL);if !(output["all_checks_passed"]as!Bool){exit(10)};print("R33 pair and tensor lifts passed")
