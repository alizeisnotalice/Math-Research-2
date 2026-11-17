import Foundation
import Metal

struct Seed: Decodable { let schema:String; let old_domain_bits:String; let relative_bent_phase_masks:[UInt16]; let seed_count:Int; let product_family_count:Int; let product_domain_bits:String }
struct Chunk: Codable { let index:Int; let start:Int; let count:Int; let host_ms:Double; let bad_abs8:Int }
struct Calibration: Codable { let phases:Int; let host_ms:Double; let cpu_ms:Double; let equal:Bool }
struct Report: Codable { let schema:String; let input_sha256:String; let source_sha256:String; let device:String; let unified_memory:Bool; let thread_width:Int; let records_expected:Int; let records_completed:Int; let chunks_expected:Int; let chunks_completed:Int; let chunk_phases:Int; let output_bytes:Int; let total_host_ms:Double; let bad_abs8:Int; let checksum_fnv64:String; let all_command_buffers_completed:Bool; let command_buffer_error:String; let calibration:Calibration; let chunks:[Chunk]; let status:String }

let metal = #"""
#include <metal_stdlib>
using namespace metal;
kernel void fwht64(device const uint2 *masks [[buffer(0)]], device int *out [[buffer(1)]], uint3 tid3 [[thread_position_in_threadgroup]], uint3 group [[threadgroup_position_in_grid]]) {
 threadgroup int a[64]; uint p=group.x,x=tid3.x;
 uint f=((x&1u)*((x>>2u)&1u))^(((x>>1u)&1u)*((x>>3u)&1u))^(((x>>4u)&1u)*((x>>5u)&1u));
 uint2 m=masks[p];uint g=x<32u?((m.x>>x)&1u):((m.y>>(x-32u))&1u);a[x]=(f^g)==0u?1:-1;
 threadgroup_barrier(mem_flags::mem_threadgroup);
 for(uint s=1u;s<64u;s<<=1u){int l=a[x],r=a[x^s];threadgroup_barrier(mem_flags::mem_threadgroup);a[x]=(x&s)==0u?l+r:r-l;threadgroup_barrier(mem_flags::mem_threadgroup);} out[p*64u+x]=a[x];
}
"""#

func die(_ s:String)->Never { FileHandle.standardError.write(Data((s+"\n").utf8));exit(2) }
func f0(_ x:Int)->Int64 { Int64((((x>>0)&1)*((x>>2)&1))^(((x>>1)&1)*((x>>3)&1))^(((x>>4)&1)*((x>>5)&1))) }
func cpuFWHT(_ mask:UInt64)->[Int64] { var a=(0..<64).map{x->Int64 in ((mask>>UInt64(x))&1)==UInt64(f0(x)) ? 1:-1};var s=1;while s<64 {var i=0;while i<64 {for j in i..<(i+s){let l=a[j],r=a[j+s];a[j]=l+r;a[j+s]=l-r};i += 2*s};s *= 2};return a }
func makeMasks(_ h:[UInt16],_ start:Int,_ count:Int)->[SIMD2<UInt32>] { let q=h.count*h.count;return (start..<(start+count)).map {i in let a=i/q,rem=i%q,h0=UInt64(h[rem/h.count]),h1=UInt64(h[rem%h.count]);var g=h0|(h0<<16)|(h1<<32)|(h1<<48);if a==1{g ^= 0xffff0000ffff0000};return SIMD2(UInt32(truncatingIfNeeded:g),UInt32(truncatingIfNeeded:g>>32))} }
func joined(_ p:SIMD2<UInt32>)->UInt64 {UInt64(p.x)|(UInt64(p.y)<<32)}

guard CommandLine.arguments.count==5 else {die("usage INPUT OUTPUT INPUTSHA SOURCESHA")}
let inputURL=URL(fileURLWithPath:CommandLine.arguments[1]),outURL=URL(fileURLWithPath:CommandLine.arguments[2]),inputSHA=CommandLine.arguments[3],sourceSHA=CommandLine.arguments[4]
let seed=try! JSONDecoder().decode(Seed.self,from:Data(contentsOf:inputURL))
guard seed.schema=="r016-four-variable-affine-seed-v1" && seed.old_domain_bits=="x0,x1,y0,y1" && seed.product_domain_bits=="x0,x1,y0,y1,z,w" && seed.seed_count==896 && seed.relative_bent_phase_masks.count==896 && Set(seed.relative_bent_phase_masks).count==896 else {die("unexpected seed input")}
let expected=2*896*896;guard seed.product_family_count==expected else{die("unexpected family count")}
guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else {die("no Metal device")}
let library=try! device.makeLibrary(source:metal,options:nil),pipeline=try! device.makeComputePipelineState(function:library.makeFunction(name:"fwht64")!)
let chunk=2048,outBytes=chunk*64*MemoryLayout<Int32>.stride;guard outBytes<=512*1024 else{die("memory cap")}
func gpu(_ m:[SIMD2<UInt32>])->([Int32],Double){let ib=m.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!};let ob=device.makeBuffer(length:m.count*64*MemoryLayout<Int32>.stride,options:.storageModeShared)!;let cb=queue.makeCommandBuffer()!,ce=cb.makeComputeCommandEncoder()!;ce.setComputePipelineState(pipeline);ce.setBuffer(ib,offset:0,index:0);ce.setBuffer(ob,offset:0,index:1);ce.dispatchThreadgroups(MTLSize(width:m.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:64,height:1,depth:1));ce.endEncoding();let t=DispatchTime.now().uptimeNanoseconds;cb.commit();cb.waitUntilCompleted();guard cb.status == .completed else{die("command buffer failed \(cb.status.rawValue) \(cb.error?.localizedDescription ?? "unknown")")};let p=ob.contents().bindMemory(to:Int32.self,capacity:m.count*64);return(Array(UnsafeBufferPointer(start:p,count:m.count*64)),Double(DispatchTime.now().uptimeNanoseconds-t)/1e6)}
let warm=makeMasks(seed.relative_bent_phase_masks,0,chunk),wg=gpu(warm),ct=DispatchTime.now().uptimeNanoseconds;var calibrationOK=true
for i in 0..<warm.count {let ref=cpuFWHT(joined(warm[i]));for x in 0..<64 where Int64(wg.0[i*64+x]) != ref[x]{calibrationOK=false}}
let calibration=Calibration(phases:chunk,host_ms:wg.1,cpu_ms:Double(DispatchTime.now().uptimeNanoseconds-ct)/1e6,equal:calibrationOK);if !calibrationOK{die("first-chunk CPU/GPU mismatch")}
var start=0,chunkIndex=0,bad=0,totalMS=0.0,hash:UInt64=1469598103934665603,completed=0,rows=[Chunk]()
while start<expected {let n=min(chunk,expected-start),m=makeMasks(seed.relative_bent_phase_masks,start,n),r=gpu(m);var chunkBad=0;for v in r.0 {if abs(Int(v)) != 8{chunkBad += 1};hash=(hash ^ UInt64(UInt32(bitPattern:v))) &* 1099511628211};bad += chunkBad;totalMS += r.1;completed += n;rows.append(Chunk(index:chunkIndex,start:start,count:n,host_ms:r.1,bad_abs8:chunkBad));start += n;chunkIndex += 1}
let report=Report(schema:seed.schema,input_sha256:inputSHA,source_sha256:sourceSHA,device:device.name,unified_memory:device.hasUnifiedMemory,thread_width:pipeline.threadExecutionWidth,records_expected:expected,records_completed:completed,chunks_expected:(expected+chunk-1)/chunk,chunks_completed:rows.count,chunk_phases:chunk,output_bytes:outBytes,total_host_ms:totalMS,bad_abs8:bad,checksum_fnv64:String(hash,radix:16),all_command_buffers_completed:true,command_buffer_error:"none",calibration:calibration,chunks:rows,status:bad==0 ? "passed":"failed")
try! JSONEncoder().encode(report).write(to:outURL);print("status=\(report.status) records=\(completed) chunks=\(rows.count) host_ms=\(totalMS) bad=\(bad) checksum=\(report.checksum_fnv64)");exit(bad==0 ? 0:1)
