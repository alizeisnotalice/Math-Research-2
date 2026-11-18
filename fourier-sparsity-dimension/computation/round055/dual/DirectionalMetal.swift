import Foundation
import Metal

let kernel = """
#include <metal_stdlib>
using namespace metal;
uint mix(uint x) { x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16); }
uint translate(uint set, uint word) {
    uint out=0;
    while(set) { uint bit=ctz(set); set &= set-1; out|=1u<<(bit^word); }
    return out;
}
void update(thread uint *record,uint rank,uint mass) {
    record[0]++;
    int raw=int(mass)-(64-(64>>rank));
    if(raw<0)record[2]++;
    if(raw==0)record[1]++;
    if(raw>0) {
        if(raw<(rank>=4?16:8))record[3]++;
        record[4]=min(record[4],uint(raw));
    }
}
void publish(threadgroup atomic_uint *work, uint rank, thread uint *record) {
    uint bin=(rank-2)*5;
    for(uint i=0;i<4;i++)
        if(record[i])atomic_fetch_add_explicit(work+bin+i,record[i],memory_order_relaxed);
    if(record[4]<999)atomic_fetch_min_explicit(work+bin+4,record[4],memory_order_relaxed);
}
struct Param { uint first, count, nconf; };
kernel void directional(device const uint *configs [[buffer(0)]],
                        device uint *output [[buffer(1)]],
                        device int *samples [[buffer(2)]],
                        constant Param &param [[buffer(3)]],
                        uint3 g [[threadgroup_position_in_grid]],
                        uint lane [[thread_index_in_threadgroup]]) {
    threadgroup int W[32];
    threadgroup uint absW[32];
    threadgroup atomic_uint metrics[25];
    uint table=mix(param.first+g.x);
    if(lane<32)W[lane]=(table>>lane)&1u ? -1 : 1;
    if(lane<25)atomic_store_explicit(metrics+lane,lane%5==4?999u:0u,memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for(uint h=1;h<32;h*=2) {
        if(lane<16) {
            uint left=((lane & ~(h-1u))<<1)|(lane&(h-1u));
            int a=W[left],b=W[left+h];
            W[left]=a+b; W[left+h]=a-b;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if(lane<32) {
        absW[lane]=uint(abs(W[lane]));
        if(g.x==0||g.x+1==param.count) samples[(g.x==0?0:32)+lane]=W[lane];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint total=0;
    for(uint w=0;w<32;w++)total+=absW[w];
    for(uint slot=lane;slot<param.nconf;slot+=128) {
        const device uint *cfg=configs+33*slot;
        uint r=cfg[0],m=1u<<r;
        uint data[32]={0}, live=0;
        for(uint w=0;w<32;w++)if(absW[w]) {
            uint z=cfg[1+w];
            data[z]+=absW[w];live|=1u<<z;
        }
        uint stored[5]={0,0,0,0,999};
        bool full=popcount(live)>m/2;
        uint af=0;
        uint rankSize=0;
        if(!full) {
            uint pivot[5]={0,0,0,0,0};
            uint o=ctz(live);
            for(uint rest=live;rest;rest&=rest-1) {
                uint x=ctz(rest)^o;
                while(x) {
                    uint pos=31u-clz(x);
                    if(pivot[pos]==0u) {pivot[pos]=x;break;}
                    x^=pivot[pos];
                }
            }
            uint span=1;
            for(uint j=0;j<5;j++)if(pivot[j])span|=translate(span,pivot[j]);
            rankSize=popcount(span);
            full=rankSize==m;
            af=translate(span,o);
        }
        if(full) {
            uint phase[5]={0,0,0,0,999};
            update(phase,r+1,total); publish(metrics,r+1,phase);
        }
        if(r>=2 && (full||rankSize==m/2)) {
            for(uint off=0;off<m;off++)
                if(full||((af>>off)&1u)==0u)update(stored,r,total-data[off]);
            publish(metrics,r,stored);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(lane<25)output[25*g.x+lane]=atomic_load_explicit(metrics+lane,memory_order_relaxed);
}
"""

struct Config: Decodable { let rank: UInt32; let labels: [UInt32] }
struct Params { var first: UInt32; var count: UInt32; var nconf: UInt32 }
struct Failure: Error { let message:String }
func check(_ value:Bool,_ explanation:String) throws {
    guard value else { throw Failure(message:explanation) }
}
func mix(_ input:UInt32)->UInt32{
    var x=input
    x ^= x >> 16; x = x &* 0x7feb352d
    x ^= x >> 15; x = x &* 0x846ca68b
    return x ^ (x >> 16)
}
func direct(_ input: UInt32)->[Int32]{
    (0..<32).map { frequency in
        Int32((0..<32).reduce(0) { total,point in
            let bit=Int((input>>point)&1) ^ ((frequency&point).nonzeroBitCount&1)
            return total+(bit==0 ? 1 : -1)
        })
    }
}

final class Slot {
    let data:MTLBuffer
    let samples:MTLBuffer
    var command: MTLCommandBuffer?
    var first=0
    var count=0
    init(_ d:MTLDevice,_ batch:Int){
        data=d.makeBuffer(length:batch*25*4,options:.storageModeShared)!
        samples=d.makeBuffer(length:64*4,options:.storageModeShared)!
    }
}
func run() throws {
    let arg=CommandLine.arguments
    try check(arg.count==6,"usage: directional_metal configurations.json output.json first count batch")
    guard let first=Int(arg[3]),let count=Int(arg[4]),let batch=Int(arg[5]),first>=0,count>0,
          batch>0,batch<=32768,first+count<=Int(UInt32.max) else {throw Failure(message:"Invalid range")}
    let configuration=try JSONDecoder().decode([Config].self,from:Data(contentsOf:URL(fileURLWithPath:arg[1])))
    try check(configuration.count==373 && configuration.allSatisfy {$0.labels.count==32},"Invalid configuration data")
    guard let d=MTLCreateSystemDefaultDevice() else {throw Failure(message:"Metal unavailable")}
    let library=try d.makeLibrary(source:kernel,options:nil)
    let pipe=try d.makeComputePipelineState(function:library.makeFunction(name:"directional")!)
    let packed=configuration.flatMap{[$0.rank]+$0.labels}
    let cfg=packed.withUnsafeBytes{d.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
    let slots=(0..<4).map{_ in Slot(d,batch)}
    let queues=(0..<2).map{_ in d.makeCommandQueue()!}
    var counts=[UInt64](repeating:0,count:25)
    for i in stride(from:4,to:25,by:5) {counts[i]=999}
    var intervals:[[Double]]=[]
    var actual=0, sampled=0
    let initTime=ProcessInfo.processInfo.systemUptime
    let shift=Date().timeIntervalSince1970-initTime
    func collect(_ slot:Slot) throws {
        guard let command=slot.command else{return}
        command.waitUntilCompleted()
        if let e=command.error{throw e}
        try check(command.status == .completed,"GPU command did not complete")
        intervals.append([command.gpuStartTime+shift,command.gpuEndTime+shift])
        let ptr=slot.data.contents().bindMemory(to:UInt32.self,capacity:slot.count*25)
        for row in 0..<slot.count {
            for i in 0..<25 {
                let value=UInt64(ptr[row*25+i])
                if i%5==4{counts[i]=min(counts[i],value)}else{counts[i]+=value}
            }
        }
        let out=slot.samples.contents().bindMemory(to:Int32.self,capacity:64)
        for (lane,id) in [slot.first,slot.first+slot.count-1].enumerated() where lane==0 || slot.count>1 {
            let coeff=direct(mix(UInt32(id)))
            try check(coeff.indices.allSatisfy{ coeff[$0]==out[lane*32+$0] },"Direct CPU/GPU Walsh mismatch")
            sampled+=32
        }
        actual+=slot.count
        slot.command=nil
    }
    var a=first,submitted=0
    while a<first+count {
        let slot=slots[submitted%slots.count]
        try collect(slot)
        let n=min(batch,first+count-a)
        let cmd=queues[submitted%queues.count].makeCommandBuffer()!
        let e=cmd.makeComputeCommandEncoder()!
        e.setComputePipelineState(pipe)
        e.setBuffer(cfg,offset:0,index:0)
        e.setBuffer(slot.data,offset:0,index:1)
        e.setBuffer(slot.samples,offset:0,index:2)
        var param=Params(first:UInt32(a),count:UInt32(n),nconf:UInt32(configuration.count))
        e.setBytes(&param,length:MemoryLayout<Params>.stride,index:3)
        e.dispatchThreadgroups(MTLSize(width:n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:128,height:1,depth:1))
        e.endEncoding()
        slot.command=cmd; slot.first=a; slot.count=n
        cmd.commit()
        a+=n;submitted+=1
    }
    for s in slots{try collect(s)}
    try check(actual==count,"GPU lost functions")
    let blocks=(2...6).map { rank -> (String,[String:UInt64]) in
        let j=(rank-2)*5
        return (String(rank),["full":counts[j],"endpoint":counts[j+1],"below_bound":counts[j+2],
                              "positive_below_requested_gap":counts[j+3],
                              "minimum_positive_units32":counts[j+4]==999 ? 0:counts[j+4]])
    }
    let result:[String:Any]=[
        "schema":"a28-r055-five-variable-full-affine-evaluation-metal-v1",
        "first_id":first,"count":count,"pid":ProcessInfo.processInfo.processIdentifier,
        "device":d.name,"gpu_command_intervals_epoch":intervals,"gpu_command_buffers":submitted,
        "gpu_command_seconds_sum":intervals.reduce(0){$0+$1[1]-$1[0]},
        "gpu_wall_seconds":ProcessInfo.processInfo.systemUptime-initTime,
        "direct_independent_cpu_coefficients":sampled,
        "by_rank":Dictionary(uniqueKeysWithValues:blocks),
        "all_checks_passed":true
    ]
    let out=try JSONSerialization.data(withJSONObject:result,options:[.prettyPrinted,.sortedKeys])
    try (out+Data([10])).write(to:URL(fileURLWithPath:arg[2]))
    print("Verified \(count) five-variable truth tables on \(d.name).")
}
do{try run()}
catch{
    FileHandle.standardError.write(Data(("R55 Metal error \(error)\n").utf8));exit(1)
}
