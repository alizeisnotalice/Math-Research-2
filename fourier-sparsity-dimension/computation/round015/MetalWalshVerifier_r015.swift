import Foundation
import Metal

struct Code: Decodable { let name:String; let rank:Int; let basis_masks:[UInt64]; let E_rows:[[Int]] }
struct Input: Decodable { let schema:String; let domain_bits:String; let records:[Code] }
struct Calibration: Codable { let phases:Int; let gpu_ms:Double; let cpu_ms:Double; let equal:Bool }
struct Row: Codable { let index:Int;let name:String;let rank:Int;let phases:Int;let chunks:Int;let gpu_ms:Double;let bad_abs8:Int;let cpu_phases:Int;let cpu_equal:Bool;let stabilizer_count:Int }
struct Report: Codable { let schema:String;let domain_bits:String;let input_sha256:String;let source_sha256:String;let backend:String;let device:String;let unified_memory:Bool;let thread_width:Int;let chunk_phases:Int;let output_bytes:Int;let mps_fallback_env:String;let all_command_buffers_completed:Bool;let command_buffer_error:String?;let calibration:Calibration;let rows:[Row];let total_gpu_ms:Double;let total_bad_abs8:Int;let total_cpu_phases:Int;let status:String }

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
func fwht(_ mask:UInt64)->[Int64] { var a=(0..<64).map{x->Int64 in ((mask>>UInt64(x))&1)==UInt64(f0(x)) ? 1:-1};var s=1;while s<64 {var i=0;while i<64 {for j in i..<(i+s){let l=a[j],r=a[j+s];a[j]=l+r;a[j+s]=l-r};i += 2*s};s *= 2};return a }
func allMasks(_ b:[UInt64],_ start:Int,_ n:Int)->[SIMD2<UInt32>] { (start..<(start+n)).map {p in var z:UInt64=0;for j in 0..<b.count where ((p>>j)&1)==1{z ^= b[j]};return SIMD2(UInt32(truncatingIfNeeded:z),UInt32(truncatingIfNeeded:z>>32)) } }
func joined(_ p:SIMD2<UInt32>)->UInt64 {UInt64(p.x)|(UInt64(p.y)<<32)}
func picks(_ n:Int,_ seed:UInt64)->Set<Int>{if n<=8{return Set(0..<n)};var q=seed|1,a=Set<Int>();while a.count<8{q=q&*6364136223846793005&+1442695040888963407;a.insert(Int(q%UInt64(n)))};return a}

guard CommandLine.arguments.count==5 else {die("usage INPUT OUTPUT INPUTSHA SOURCESHA")}
let inputURL=URL(fileURLWithPath:CommandLine.arguments[1]),outURL=URL(fileURLWithPath:CommandLine.arguments[2]),inputSHA=CommandLine.arguments[3],sourceSHA=CommandLine.arguments[4]
let input=try! JSONDecoder().decode(Input.self,from:Data(contentsOf:inputURL))
guard input.schema=="r015-q4-selector-affine-v1" && input.domain_bits=="x,z0,y,w0,z1,w1" && input.records.count==3 else {die("unexpected schema")}
guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else {die("no Metal device")}
let library=try! device.makeLibrary(source:metal,options:nil)
let pipeline=try! device.makeComputePipelineState(function:library.makeFunction(name:"fwht64")!)
let chunk=2048,outBytes=chunk*64*MemoryLayout<Int32>.stride
guard outBytes<=512*1024 else {die("output budget")}
func gpu(_ m:[SIMD2<UInt32>])->([Int32],Double){let ib=m.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!};let ob=device.makeBuffer(length:m.count*64*MemoryLayout<Int32>.stride,options:.storageModeShared)!;let cb=queue.makeCommandBuffer()!,ce=cb.makeComputeCommandEncoder()!;ce.setComputePipelineState(pipeline);ce.setBuffer(ib,offset:0,index:0);ce.setBuffer(ob,offset:0,index:1);ce.dispatchThreadgroups(MTLSize(width:m.count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:64,height:1,depth:1));ce.endEncoding();let t=DispatchTime.now().uptimeNanoseconds;cb.commit();cb.waitUntilCompleted();guard cb.status == .completed else {die("command buffer failed \(cb.status.rawValue) \(cb.error?.localizedDescription ?? "unknown")")};let p=ob.contents().bindMemory(to:Int32.self,capacity:m.count*64);return(Array(UnsafeBufferPointer(start:p,count:m.count*64)),Double(DispatchTime.now().uptimeNanoseconds-t)/1e6)}
let warm=allMasks(input.records[0].basis_masks,0,min(8,1<<input.records[0].rank)),wg=gpu(warm),ct=DispatchTime.now().uptimeNanoseconds;var calibrationOK=true
for i in 0..<warm.count {let c=fwht(joined(warm[i]));for x in 0..<64 where Int64(wg.0[i*64+x]) != c[x]{calibrationOK=false}}
let calibration=Calibration(phases:warm.count,gpu_ms:wg.1,cpu_ms:Double(DispatchTime.now().uptimeNanoseconds-ct)/1e6,equal:calibrationOK);if !calibrationOK {die("calibration mismatch")}
var rows=[Row](),totalGPU=0.0,totalBad=0,totalCPU=0,good=true
for (i,c) in input.records.enumerated(){guard c.rank==c.basis_masks.count && c.rank<=15 else{die("bad rank")};let phases=1<<c.rank,need=Set(0..<phases),all=allMasks(c.basis_masks,0,phases),set=Set(all.map(joined));var stab=0;for a in set{if c.basis_masks.allSatisfy({set.contains(a & $0)}){stab += 1}};var start=0,chunks=0,bad=0,cpuN=0,cpuOK=true,ms=0.0;while start<phases{let n=min(chunk,phases-start),mm=Array(all[start..<(start+n)]),r=gpu(mm);ms += r.1;chunks += 1;for v in r.0 where abs(Int(v)) != 8{bad += 1};for j in 0..<n where need.contains(start+j){let ref=fwht(joined(mm[j]));cpuN += 1;for x in 0..<64 where Int64(r.0[j*64+x]) != ref[x]{cpuOK=false}};start += n};if bad != 0 || !cpuOK || stab != 2{good=false};rows.append(Row(index:i,name:c.name,rank:c.rank,phases:phases,chunks:chunks,gpu_ms:ms,bad_abs8:bad,cpu_phases:cpuN,cpu_equal:cpuOK,stabilizer_count:stab));totalGPU += ms;totalBad += bad;totalCPU += cpuN}
let report=Report(schema:input.schema,domain_bits:input.domain_bits,input_sha256:inputSHA,source_sha256:sourceSHA,backend:"Swift Metal int32 threadgroup FWHT",device:device.name,unified_memory:device.hasUnifiedMemory,thread_width:pipeline.threadExecutionWidth,chunk_phases:chunk,output_bytes:outBytes,mps_fallback_env:ProcessInfo.processInfo.environment["PYTORCH_ENABLE_MPS_FALLBACK"] ?? "unset",all_command_buffers_completed:true,command_buffer_error:nil,calibration:calibration,rows:rows,total_gpu_ms:totalGPU,total_bad_abs8:totalBad,total_cpu_phases:totalCPU,status:good ? "passed":"failed")
try! JSONEncoder().encode(report).write(to:outURL);print("status=\(report.status) rows=\(rows.count) gpu_ms=\(totalGPU) bad=\(totalBad) cpu=\(totalCPU)");exit(good ? 0:1)

