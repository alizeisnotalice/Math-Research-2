import Foundation
import Metal

let source = """
#include <metal_stdlib>
using namespace metal;

uint normal(uint id) {
    return ((id&1u)<<3)|((id&14u)<<4)|((id&2032u)<<5)|((id&67106816u)<<6);
}
int fw(int s, uint lane) {
    for(ushort h=1;h<32;h*=2) {
        int partner=simd_shuffle_xor(s,h);
        s=(lane&h)?partner-s:partner+s;
    }
    return s;
}

struct Params { uint first, count, phases; };

kernel void norms(device const uint *affine [[buffer(0)]],
                  device ushort *norms [[buffer(1)]],
                  device int *samples [[buffer(2)]],
                  constant Params &p [[buffer(3)]],
                  uint global [[thread_position_in_grid]],
                  uint lane [[thread_index_in_simdgroup]]) {
    uint group=global>>5;
    if(group>=p.count)return;
    uint id=p.first+group;
    uint mask=normal(id);
    int value=fw(((mask>>lane)&1u)?-1:1,lane);
    uint total=simd_sum(uint(abs(value)));
    if(group==0||group+1u==p.count) samples[(group==0?0:32)+lane]=value;
    if(lane==0) norms[group]=ushort(total);
}
kernel void orbit_norms(device const uint *affine [[buffer(0)]],
                        device ushort *norms [[buffer(1)]],
                        device int *samples [[buffer(2)]],
                        constant Params &p [[buffer(3)]],
                        uint3 grid [[threadgroup_position_in_grid]],
                        uint phase_group [[simdgroup_index_in_threadgroup]],
                        uint lane [[thread_index_in_simdgroup]],
                        uint tid [[thread_index_in_threadgroup]]) {
    threadgroup uint low[8],high[8];
    uint id=p.first+grid.x;
    uint mask=normal(id);
    uint lower=380,upper=0;
    for(uint phase=phase_group;phase<64;phase+=8) {
        uint word=mask^affine[phase];
        int value=fw(((word>>lane)&1u)?-1:1,lane);
        uint total=simd_sum(uint(abs(value)));
        lower=min(lower,total);upper=max(upper,total);
        if(phase==0&&(grid.x==0||grid.x+1u==p.count))
            samples[(grid.x==0?0:32)+lane]=value;
    }
    if(lane==0) { low[phase_group]=lower; high[phase_group]=upper; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(tid==0) {
        for(uint j=1;j<8;j++){lower=min(lower,low[j]);upper=max(upper,high[j]);}
        norms[grid.x]=(lower==upper)?ushort(lower):ushort(999);
    }
}
"""

struct Params { var first:UInt32; var count:UInt32; var phases:UInt32 }
struct Failure: Error { let message:String }
func check(_ x:Bool,_ why:String) throws { guard x else {throw Failure(message:why)} }

func normal(_ id:UInt32)->UInt32 {
    ((id&1)<<3) | ((id&14)<<4) | ((id&2032)<<5) | ((id&67106816)<<6)
}
func checksum(_ id:UInt64,_ norm:UInt64)->UInt64 {
    var x=(id<<8)|norm
    x=(x^(x>>30)) &* 0xbf58476d1ce4e5b9
    x=(x^(x>>27)) &* 0x94d049bb133111eb
    return x^(x>>31)
}
func direct(_ id:Int)->[Int32] {
    let mask=normal(UInt32(id))
    return (0..<32).map { w in
        Int32((0..<32).reduce(0){ sum,x in
            let bit = Int((mask>>x)&1)^((x&w).nonzeroBitCount&1)
            return sum+(bit==0 ? 1:-1)
        })
    }
}

final class Slot {
    let norms:MTLBuffer
    let samples:MTLBuffer
    var command:MTLCommandBuffer?
    var first=0
    var count=0
    init(_ d:MTLDevice,_ n:Int) {
        norms=d.makeBuffer(length:n*2,options:.storageModeShared)!
        samples=d.makeBuffer(length:64*4,options:.storageModeShared)!
    }
}

func run() throws {
    let arg=CommandLine.arguments
    try check(arg.count==6,"usage: flat_metal output.json first count batch phases")
    guard let first=Int(arg[2]),let count=Int(arg[3]),let batch=Int(arg[4]),let phases=Int(arg[5]),
          first>=0,count>0,first+count<=1<<26,batch>0,batch<=262144,phases==1||phases==64 else {
        throw Failure(message:"Invalid corpus range")
    }
    guard let device=MTLCreateSystemDefaultDevice() else {throw Failure(message:"No Metal device")}
    let lib=try device.makeLibrary(source:source,options:nil)
    let function=lib.makeFunction(name:phases==1 ? "norms":"orbit_norms")!
    let p=try device.makeComputePipelineState(function:function)
    try check(p.threadExecutionWidth==32 && p.maxTotalThreadsPerThreadgroup>=256,"Unexpected device SIMD width")
    var columns=[UInt32](repeating:0,count:64)
    for b in 0..<64 {
        for x in 0..<32 where ((b&x).nonzeroBitCount&1)^(b>>5) != 0 {
            columns[b] |= UInt32(1)<<x
        }
    }
    let affine=columns.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
    let slots=(0..<4).map { _ in Slot(device,batch) }
    let queues=(0..<2).map { _ in device.makeCommandQueue()! }
    var histogram=[UInt64](repeating:0,count:193)
    var below:[[Int]]=[]
    var keyed:UInt64=0
    var sampled=0
    var generated=0
    var gpuSeconds:Double=0
    var intervals:[[Double]]=[]
    let uptime=ProcessInfo.processInfo.systemUptime
    let shift=Date().timeIntervalSince1970-uptime
    func collect(_ slot:Slot) throws {
        guard let command=slot.command else{return}
        command.waitUntilCompleted()
        if let error=command.error {throw error}
        try check(command.status == .completed,"GPU queue failed")
        gpuSeconds += command.gpuEndTime-command.gpuStartTime
        intervals.append([command.gpuStartTime+shift,command.gpuEndTime+shift])
        let bins=slot.norms.contents().bindMemory(to:UInt16.self,capacity:slot.count)
        for i in 0..<slot.count {
            let norm=Int(bins[i])
            try check(norm<193,"Out-of-range Walsh norm or modulation invariance violation")
            histogram[norm]+=1
            keyed ^= checksum(UInt64(slot.first+i),UInt64(norm))
            if norm<96{below.append([slot.first+i,norm])}
        }
        let data=slot.samples.contents().bindMemory(to:Int32.self,capacity:64)
        for (offset,index) in [slot.first,slot.first+slot.count-1].enumerated() where offset==0 || slot.count>1 {
            let W=direct(index)
            try check(W.indices.allSatisfy {W[$0]==data[32*offset+$0]},"Direct CPU/GPU coefficients differ")
            sampled+=32
        }
        generated+=slot.count
        slot.command=nil
    }
    var start=first, submitted=0
    let threads=MTLSize(width:256,height:1,depth:1)
    while start<first+count {
        let slot=slots[submitted%slots.count]
        try collect(slot)
        let length=min(batch,first+count-start)
        let cmd=queues[submitted%queues.count].makeCommandBuffer()!
        let e=cmd.makeComputeCommandEncoder()!
        e.setComputePipelineState(p)
        e.setBuffer(affine,offset:0,index:0)
        e.setBuffer(slot.norms,offset:0,index:1)
        e.setBuffer(slot.samples,offset:0,index:2)
        var options=Params(first:UInt32(start),count:UInt32(length),phases:UInt32(phases))
        e.setBytes(&options,length:MemoryLayout<Params>.stride,index:3)
        if phases==1 {
            e.dispatchThreads(MTLSize(width:32*length,height:1,depth:1),threadsPerThreadgroup:threads)
        } else {
            e.dispatchThreadgroups(MTLSize(width:length,height:1,depth:1),threadsPerThreadgroup:threads)
        }
        e.endEncoding()
        slot.command=cmd;slot.first=start;slot.count=length
        cmd.commit()
        start+=length;submitted+=1
    }
    for slot in slots{try collect(slot)}
    try check(generated==count && histogram.reduce(0,+)==UInt64(count),"GPU enumeration lost candidates")
    let h=Dictionary(uniqueKeysWithValues:histogram.enumerated().filter{$0.element != 0}.map{(String($0.offset),$0.element)})
    let r:[String:Any]=[
        "schema":"a28-r056-five-variable-exhaustive-metal-v1",
        "first_id":first,"count":count,"device":device.name,"pid":ProcessInfo.processInfo.processIdentifier,
        "simd_width":p.threadExecutionWidth,"affine_variants_transformed_for_each_representative":phases,
        "distinct_raw_truth_tables_transformed":count*phases,"walsh_coefficients_evaluated":count*phases*32,
        "norm_position_checksum64":String(format:"%016llx",keyed),
        "norm_histogram_units32":h,
        "below_three_id_and_norm":below.sorted{$0[0]<$1[0]},
        "direct_CPU_walsh_coefficients_crosschecked":sampled,
        "gpu_metal_command_seconds_summed":gpuSeconds,"gpu_command_intervals_epoch":intervals,
        "gpu_metal_command_buffers":submitted,"gpu_wall_seconds":ProcessInfo.processInfo.systemUptime-uptime,
        "all_checks_passed":true
    ]
    let out=try JSONSerialization.data(withJSONObject:r,options:[.prettyPrinted,.sortedKeys])
    try (out+Data([10])).write(to:URL(fileURLWithPath:arg[1]))
    print("Metal scan on \(device.name) passed \(count*phases) actual functions.")
}
do { try run() }
catch {
    FileHandle.standardError.write(Data(("R56 Metal error: \(error)\n").utf8))
    exit(1)
}
