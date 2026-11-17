import Foundation
import Metal

struct Case: Decodable { let name:String; let bits:String; let f0_mask:UInt64; let basis_masks:[UInt64]; let dimension:Int; let phase_count:Int; let blocks:[[Int]] }
struct Input: Decodable { let schema:String; let cases:[Case] }
struct Chunk: Codable { let start:Int; let count:Int; let host_ms:Double; let bad_abs8:Int }
struct CaseReport: Codable { let name:String; let dimension:Int; let phases:Int; let rank_ok:Bool; let blocks_partition_ok:Bool; let blocks_disjoint_ok:Bool; let gpu_bad_abs8:Int; let cpu_entrywise_equal:Bool; let cpu_phases:Int; let chunks:[Chunk] }
struct Report: Codable { let schema:String; let input_sha256:String; let source_sha256:String; let device:String; let unified_memory:Bool; let thread_width:Int; let chunk_phases:Int; let output_bytes:Int; let total_phases:Int; let total_host_ms:Double; let all_command_buffers_completed:Bool; let command_buffer_error:String; let cases:[CaseReport]; let status:String }

let metal = #"""
#include <metal_stdlib>
using namespace metal;
kernel void fwht64(device const uint2 *masks [[buffer(0)]], device int *out [[buffer(1)]], uint3 tid3 [[thread_position_in_threadgroup]], uint3 group [[threadgroup_position_in_grid]]) {
 threadgroup int a[64];uint p=group.x,x=tid3.x;uint2 m=masks[p];uint b=x<32u?((m.x>>x)&1u):((m.y>>(x-32u))&1u);a[x]=b==0u?1:-1;
 threadgroup_barrier(mem_flags::mem_threadgroup);
 for(uint s=1u;s<64u;s<<=1u){int l=a[x],r=a[x^s];threadgroup_barrier(mem_flags::mem_threadgroup);a[x]=(x&s)==0u?l+r:r-l;threadgroup_barrier(mem_flags::mem_threadgroup);}out[p*64u+x]=a[x];
}
"""#
func die(_ s:String)->Never{FileHandle.standardError.write(Data((s+"\n").utf8));exit(2)}
func rank(_ a:[UInt64])->Int{var p=[Int:UInt64]();for var x in a{while x != 0{let q=63-x.leadingZeroBitCount;if let y=p[q]{x ^= y}else{p[q]=x;break}}};return p.count}
func cpu(_ m:UInt64)->[Int64]{var a=(0..<64).map{((m>>UInt64($0))&1)==0 ? Int64(1):Int64(-1)};var s=1;while s<64{var i=0;while i<64{for j in i..<(i+s){let l=a[j],r=a[j+s];a[j]=l+r;a[j+s]=l-r};i += 2*s};s *= 2};return a}
func masks(_ c:Case,_ start:Int,_ n:Int)->[SIMD2<UInt32>]{(start..<(start+n)).map{phase in var z=c.f0_mask;for j in 0..<c.dimension where ((phase>>j)&1)==1{z ^= c.basis_masks[j]};return SIMD2(UInt32(truncatingIfNeeded:z),UInt32(truncatingIfNeeded:z>>32))}}
func joined(_ p:SIMD2<UInt32>)->UInt64{UInt64(p.x)|(UInt64(p.y)<<32)}
func blocksOK(_ b:[[Int]])->Bool{let v=b.flatMap{$0};return v.count==64 && Set(v).count==64 && Set(v)==Set(0..<64)}
guard CommandLine.arguments.count==5 else{die("usage INPUT OUTPUT INPUTSHA SOURCESHA")}
let inputURL=URL(fileURLWithPath:CommandLine.arguments[1]),outURL=URL(fileURLWithPath:CommandLine.arguments[2]),inputSHA=CommandLine.arguments[3],sourceSHA=CommandLine.arguments[4]
let input=try! JSONDecoder().decode(Input.self,from:Data(contentsOf:inputURL))
guard input.schema=="a28-r017-nonparallel-exact-v1" && input.cases.count==2 else{die("unexpected input")}
guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else{die("no Metal device")}
let library=try! device.makeLibrary(source:metal,options:nil),pipeline=try! device.makeComputePipelineState(function:library.makeFunction(name:"fwht64")!)
let chunk=2048,outBytes=chunk*64*MemoryLayout<Int32>.stride;guard outBytes<=512*1024 else{die("memory cap")}
func gpu(_ m:[SIMD2<UInt32>])->([Int32],Double){let ib=m.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!};let ob=device.makeBuffer(length:m.count*64*MemoryLayout<Int32>.stride,options:.storageModeShared)!;let cb=queue.makeCommandBuffer()!,ce=cb.makeComputeCommandEncoder()!;ce.setComputePipelineState(pipeline);ce.setBuffer(ib,offset:0,index:0);ce.setBuffer(ob,offset:0,index:1);ce.dispatchThreadgroups(MTLSize(width:m.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:64,height:1,depth:1));ce.endEncoding();let t=DispatchTime.now().uptimeNanoseconds;cb.commit();cb.waitUntilCompleted();guard cb.status == .completed else{die("command buffer \(cb.status.rawValue) \(cb.error?.localizedDescription ?? "error")")};let p=ob.contents().bindMemory(to:Int32.self,capacity:m.count*64);return(Array(UnsafeBufferPointer(start:p,count:m.count*64)),Double(DispatchTime.now().uptimeNanoseconds-t)/1e6)}
var reports=[CaseReport](),total=0,totalMS=0.0,good=true
for c in input.cases{let rankOK=rank(c.basis_masks)==c.dimension && c.dimension==c.basis_masks.count && c.phase_count==(1<<c.dimension);let blocks=blocksOK(c.blocks);if !rankOK || !blocks{good=false};var start=0,bad=0,cpuOK=true,cpuN=0,chunks=[Chunk]();while start<c.phase_count{let n=min(chunk,c.phase_count-start),m=masks(c,start,n),r=gpu(m);var thisBad=0;for v in r.0 where abs(Int(v)) != 8{thisBad += 1};for j in 0..<n{let ref=cpu(joined(m[j]));cpuN += 1;for x in 0..<64 where Int64(r.0[j*64+x]) != ref[x]{cpuOK=false}};chunks.append(Chunk(start:start,count:n,host_ms:r.1,bad_abs8:thisBad));bad += thisBad;totalMS += r.1;start += n};if bad != 0 || !cpuOK{good=false};reports.append(CaseReport(name:c.name,dimension:c.dimension,phases:c.phase_count,rank_ok:rankOK,blocks_partition_ok:blocks,blocks_disjoint_ok:blocks,gpu_bad_abs8:bad,cpu_entrywise_equal:cpuOK,cpu_phases:cpuN,chunks:chunks));total += c.phase_count}
let report=Report(schema:input.schema,input_sha256:inputSHA,source_sha256:sourceSHA,device:device.name,unified_memory:device.hasUnifiedMemory,thread_width:pipeline.threadExecutionWidth,chunk_phases:chunk,output_bytes:outBytes,total_phases:total,total_host_ms:totalMS,all_command_buffers_completed:true,command_buffer_error:"none",cases:reports,status:good ? "passed":"failed")
try! JSONEncoder().encode(report).write(to:outURL);print("status=\(report.status) phases=\(total) host_ms=\(totalMS)");exit(good ? 0:1)
