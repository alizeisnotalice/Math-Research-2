import Foundation
import Metal
import CryptoKit

let n = 16, q = 4, truthCount = 1 << n
let kernelSource = """
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]], device int *output [[buffer(1)]],
                       constant uint &stride [[buffer(2)]], uint gid [[thread_position_in_grid]]) {
    uint local = gid & 15u, base = gid - local, mate = base + (local ^ stride);
    int left = input[gid], right = input[mate];
    output[gid] = (local & stride) ? (right-left) : (left+right);
}
"""
func parity(_ x: Int) -> Int { x.nonzeroBitCount & 1 }
func sha256(_ url: URL) -> String { SHA256.hash(data: try! Data(contentsOf: url)).map { String(format:"%02x",$0) }.joined() }
func directWalsh(_ masks: [UInt16]) -> [Int32] {
    var out=[Int32](repeating:0,count:masks.count*n)
    for p in 0..<masks.count { for xi in 0..<n {
        var s=0; for x in 0..<n { s += (Int((masks[p] >> UInt16(x)) & 1) ^ parity(x&xi)) == 0 ? 1 : -1 }
        out[p*n+xi]=Int32(s)
    } }
    return out
}
func metalFWHT(_ device: MTLDevice, _ pipeline: MTLComputePipelineState, _ masks: [UInt16]) throws -> [Int32] {
    var values=[Int32](repeating:0,count:masks.count*n)
    for p in 0..<masks.count { for x in 0..<n { values[p*n+x]=((masks[p] >> UInt16(x))&1)==0 ? 1:-1 } }
    let length=values.count*MemoryLayout<Int32>.stride, options: MTLResourceOptions = .storageModeShared
    let first=values.withUnsafeBytes { device.makeBuffer(bytes:$0.baseAddress!,length:length,options:options)! }
    let second=device.makeBuffer(length:length,options:options)!; var input=first, output=second
    let command=device.makeCommandQueue()!.makeCommandBuffer()!
    for h:UInt32 in [1,2,4,8] {
        let e=command.makeComputeCommandEncoder()!; e.setComputePipelineState(pipeline)
        e.setBuffer(input,offset:0,index:0);e.setBuffer(output,offset:0,index:1);var stride=h;e.setBytes(&stride,length:4,index:2)
        let width=min(256,pipeline.maxTotalThreadsPerThreadgroup)
        e.dispatchThreads(MTLSize(width:masks.count*n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:width,height:1,depth:1));e.endEncoding();swap(&input,&output)
    }
    command.commit();command.waitUntilCompleted();guard command.status == .completed else { throw NSError(domain:"metal",code:Int(command.status.rawValue)) }
    return Array(UnsafeBufferPointer(start:input.contents().bindMemory(to:Int32.self,capacity:values.count),count:values.count))
}
func liftW(_ h: Int) -> UInt16 {
    var out:UInt16=0;for x in 0..<n where ((h>>(x>>2))&1)==1 {out |= UInt16(1)<<UInt16(x)};return out
}
func isMM(_ m: UInt16) -> Bool {
    var slopes=[Int]()
    for w in 0..<q {
        var support=[Int]()
        for a in 0..<q { var s=0;for z in 0..<q {let bit=Int((m>>UInt16(z+(w<<2)))&1)^parity(z&a);s += bit==0 ? 1:-1};if s != 0 {if abs(s) != q{return false};support.append(a)} }
        if support.count != 1{return false};slopes.append(support[0])
    }
    return Set(slopes).count == q
}
func permutations(_ a:[Int]) -> [[Int]] { if a.count<2{return [a]};var out=[[Int]]();for i in 0..<a.count {var rest=a;let x=rest.remove(at:i);for p in permutations(rest){out.append([x]+p)}};return out }
func formulaCounterexample() -> UInt16 { var g:UInt16=0;for w in 0..<q {for z in 0..<q {let bit=parity(z&w)^((z&1)*((z>>1)&1));if bit==1{g|=UInt16(1)<<UInt16(z+(w<<2))}}};return g }

let args=CommandLine.arguments;guard args.count==2 else{fatalError("usage: R25SidonMetal output.json")}
guard let device=MTLCreateSystemDefaultDevice() else{fatalError("Metal unavailable")}
let library=try!device.makeLibrary(source:kernelSource,options:nil);let pipeline=try!device.makeComputePipelineState(function:library.makeFunction(name:"fwht_stage")!)
let begun=Date();let masks=(0..<truthCount).map{UInt16($0)}
let cpu=directWalsh(masks),gpu=try!metalFWHT(device,pipeline,masks)
var mismatch=0,bent=[Bool](repeating:false,count:truthCount),bentCount=0
for p in 0..<truthCount {var good=true;for xi in 0..<n {let i=p*n+xi;if cpu[i] != gpu[i]{mismatch += 1};if Int(abs(cpu[i])) != q{good=false}};bent[p]=good;if good{bentCount += 1}}
let hMasks=(0..<16).map(liftW)
var compatible=[UInt16]()
for m in 0..<truthCount where bent[m] {
    let word = UInt16(m)
    if hMasks.allSatisfy({ h in bent[Int(word ^ h)] }) { compatible.append(word) }
}
let mmRecognized: [UInt16] = (0..<truthCount).filter { isMM(UInt16($0)) }.map { UInt16($0) }
var constructed=Set<UInt16>();for permutation in permutations([0,1,2,3]) {for offset in 0..<16 {var word:UInt16=0;for w in 0..<q {for z in 0..<q {if (parity(z&permutation[w])^((offset>>w)&1))==1{word|=UInt16(1)<<UInt16(z+(w<<2))}}};constructed.insert(word)}}
let g=formulaCounterexample(), expectedG:UInt16=0xe428
let affineW=[0,15,10,5,12,3,6,9]
let counterAffine=affineW.allSatisfy{bent[Int(g^liftW($0))]}
let counterDetected = !bent[Int(g^liftW(0x0e))]
let complete=Set(compatible),recognized=Set(mmRecognized)
let output:[String:Any]=[
 "schema":"a28-r025-sidon-metal-v1","all_truth_tables":truthCount,"signed_walsh_coefficients_checked":truthCount*n,
 "cpu_oracle":"direct 16x16 character sums","gpu_backend":"Metal int32 global staged ping-pong FWHT, four stages","device":device.name,"registry_id":String(device.registryID),
 "gpu_cpu_mismatch_count":mismatch,"bent_count":bentCount,"H":"all 16 functions of w alone","compatible_masks_hex":compatible.map{String(format:"0x%04x",$0)},"compatible_count":compatible.count,
 "fixed_split_MM_count_by_recognizer":mmRecognized.count,"fixed_split_MM_count_constructed":constructed.count,"compatible_equals_MM":complete==recognized && complete==constructed,
 "counterexample_mask_hex":String(format:"0x%04x",g),"counterexample_formula_matches_e428":g==expectedG,"counterexample_affine_w_all_bent":counterAffine,"counterexample_not_MM":!isMM(g),"counterexample_w_mask_0e_detects_illegal":counterDetected,
 "all_checks_passed":mismatch==0 && bentCount==896 && compatible.count==384 && complete==recognized && complete==constructed && counterAffine && !isMM(g) && counterDetected,
 "elapsed_seconds":Date().timeIntervalSince(begun)]
try!JSONSerialization.data(withJSONObject:output,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:args[1]));if !(output["all_checks_passed"] as!Bool){exit(10)};print("R25 passed 65536 truth tables")
