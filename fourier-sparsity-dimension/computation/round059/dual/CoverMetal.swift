import Foundation
import Metal

let source = """
#include <metal_stdlib>
using namespace metal;
struct Cover { uint masks[48],weights[48],m,N,n,family,target; };
uint mix(uint x) {x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);}
uint next_word(thread uint &s) {s=s*1664525u+1013904223u;return s>>12;}
void add(thread Cover &c,uint mask,uint weight){
    for(uint i=0;i<c.m;i++) if(c.masks[i]==mask){c.weights[i]+=weight;return;}
    c.masks[c.m]=mask;c.weights[c.m]=weight;c.m+=1;
}
void partition(thread Cover &c,uint id,uint w,device const uint *bank){
    uint base=17u*((c.n-3u)*64u+id);
    for(uint i=1;i<=bank[base];i++)add(c,bank[base+i],w);
}
void cycles(thread Cover &c,thread uint &seed){
    uint permutation[16];
    uint sizes[4]={3,5,3,5};
    for(uint i=0;i<c.N;i++)permutation[i]=i;
    for(uint i=c.N-1;i>0;i--){
        uint j=next_word(seed)%(i+1),t=permutation[i];permutation[i]=permutation[j];permutation[j]=t;
    }
    uint cursor=0;
    for(uint i=0;i<c.N/4;i++){
        uint count=sizes[i];
        for(uint j=0;j<count;j++)
            add(c,(1u<<permutation[cursor+j])|(1u<<permutation[cursor+(j+1)%count]),1u);
        cursor+=count;
    }
}
Cover generate(uint id,device const uint *bank){
    Cover c;
    c.m=0;c.n=3+((id>>2)&1u);c.family=id&3u;c.target=0;
    if(id<4)c.n=3+(id&1u);
    c.N=1u<<c.n;
    uint seed=mix(id^0xa2805901u);
    if(id<2){
        for(uint j=0;j<c.n;j++){
            uint a=0;
            for(uint x=0;x<c.N;x++)if((x>>j)&1u)a|=1u<<x;
            add(c,a,j+1);add(c,((1u<<c.N)-1)^a,j+1);c.target+=j+1;
        }
    } else if(id<4||c.family==2){
        cycles(c,seed);c.target=2;
    } else if(c.family==3){
        cycles(c,seed);uint i=next_word(seed)%64;partition(c,i,2,bank);c.target=4;
    } else {
        for(uint j=0;j<2+c.family;j++){
            uint i=next_word(seed)%64;
            uint w=next_word(seed)%5+1;
            partition(c,i,w,bank);c.target+=w;
        }
    }
    return c;
}
uint binary_add(uint x,thread uint *pivot){
    while(x){
        uint i=31-clz(x);
        if(!pivot[i]){pivot[i]=x;return 1;}
        x^=pivot[i];
    }
    return 0;
}
uint reduce(ulong z) {
    const ulong p=2147483647ul;
    ulong r=(z&p)+(z>>31);r=(r&p)+(r>>31);
    return uint(r>=p?r-p:r);
}
uint real_rank(thread const Cover &c){
    const uint p=2147483647u;
    uint e[225];
    for(uint j=0;j<225;j++)e[j]=0;
    uint rank=0, cols=c.N-1;
    for(uint row=0;row<c.m;row++){
        uint u[15];
        int base=int(c.masks[row]&1u);
        for(uint j=0;j<cols;j++){
            int v=int((c.masks[row]>>(j+1))&1u)-base;
            u[j]=v<0?p-1:uint(v);
        }
        for(uint j=0;j<cols;j++)if(u[j]){
            uint a=e[j*15+j],b=u[j];
            if(!a){for(uint k=j;k<cols;k++)e[j*15+k]=u[k];rank+=1;break;}
            for(uint k=j;k<cols;k++)u[k]=reduce(ulong(u[k])*a+ulong(p-b)*e[j*15+k]);
        }
    }
    return rank;
}
struct Params { uint first,length; };
kernel void check_covers(device const uint *bank [[buffer(0)]],
                         device uint *results [[buffer(1)]],
                         constant Params &p [[buffer(2)]],
                         uint threadID [[thread_position_in_grid]]){
    if(threadID>=p.length)return;
    Cover c=generate(p.first+threadID,bank);
    int spectrum[16]={0};
    uint directions[4]={0};
    uint t=0,maximum=0;
    for(uint i=0;i<c.m;i++){
        for(uint v=0;v<c.N;v++){
            int coefficient=0;
            for(uint z=c.masks[i];z;z&=z-1)
                coefficient+=popcount(ctz(z)&v)&1u ? -1:1;
            spectrum[v]+=coefficient*int(c.weights[i]);
            if(v!=0&&coefficient)t+=binary_add(v,directions);
        }
    }
    for(uint x=0;x<c.N;x++){
        uint active=0;
        for(uint i=0;i<c.m;i++)active+=(c.masks[i]>>x)&1u;
        maximum=max(maximum,active);
    }
    uint r=real_rank(c);
    uint code=c.m|(maximum<<6)|(t<<12)|(r<<15)|((c.n-3)<<20)|(c.family<<21);
    bool valid=spectrum[0]==int(c.N*c.target)&&t<=r&&t<=c.m-maximum;
    for(uint v=1;v<c.N;v++)valid &= spectrum[v]==0;
    results[threadID]=valid?code:0xffffffffu;
}
"""

struct Params { var first:UInt32; var length:UInt32 }
struct Failure:Error { let message:String }
func check(_ condition:Bool,_ description:String)throws {
    guard condition else {throw Failure(message:description)}
}
func checksum(_ id:Int,_ code:UInt32)->UInt64 {
    var x=(UInt64(id)<<24)|UInt64(code)
    x=(x^(x>>30)) &* 0xbf58476d1ce4e5b9
    x=(x^(x>>27)) &* 0x94d049bb133111eb
    return x^(x>>31)
}
final class Slot {
    let buffer:MTLBuffer
    var command:MTLCommandBuffer?
    var head=0,length=0
    init(_ device:MTLDevice,_ batch:Int){buffer=device.makeBuffer(length:4*batch,options:.storageModeShared)!}
}
func run()throws {
    let args=CommandLine.arguments
    try check(args.count==6,"usage: cover_metal bank.bin result.json first count batch")
    guard let first=Int(args[3]),let count=Int(args[4]),let batch=Int(args[5]),
          first>=0,count>0,batch>0,batch<=32768 else{throw Failure(message:"Invalid range")}
    let bank=try Data(contentsOf:URL(fileURLWithPath:args[1]))
    try check(bank.count==2*64*17*4,"Bad candidate bank")
    guard let device=MTLCreateSystemDefaultDevice()else{throw Failure(message:"No Metal")}
    let pipe=try device.makeComputePipelineState(function:device.makeLibrary(source:source,options:nil).makeFunction(name:"check_covers")!)
    let queues=(0..<2).map{_ in device.makeCommandQueue()!}
    let ring=(0..<4).map{_ in Slot(device,batch)}
    let sourceData=bank.withUnsafeBytes{device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
    var hist=Array(repeating:Array(repeating:0,count:4),count:2)
    var overlap=0,tightActivity=0,tightRank=0,cases=0
    var sum:UInt64=0
    var time:Double=0
    var header:[String:UInt32]=[:]
    var intervals:[[Double]]=[]
    let clock=ProcessInfo.processInfo.systemUptime
    let offset=Date().timeIntervalSince1970-clock
    func read(_ slot:Slot)throws{
        guard let c=slot.command else{return}
        c.waitUntilCompleted()
        if let error=c.error{throw error}
        try check(c.status == .completed,"GPU failed")
        time+=c.gpuEndTime-c.gpuStartTime
        intervals.append([c.gpuStartTime+offset,c.gpuEndTime+offset])
        let ptr=slot.buffer.contents().bindMemory(to:UInt32.self,capacity:slot.length)
        for i in 0..<slot.length {
            let code=ptr[i]
            try check(code != UInt32.max,"Metal invalid weighted cover or inequality")
            let m=Int(code&63),maximum=Int((code>>6)&63),t=Int((code>>12)&7),r=Int((code>>15)&15)
            let dimension=Int((code>>20)&1), family=Int((code>>21)&3)
            let id=slot.head+i
            if id<4{ header[String(id)]=code }
            if maximum>1{overlap+=1}
            if t==m-maximum{tightActivity+=1}
            if t==r{tightRank+=1}
            hist[dimension][family]+=1
            sum^=checksum(id,code)
        }
        cases+=slot.length
        slot.command=nil
    }
    var next=first,submitted=0
    while next<first+count {
        let slot=ring[submitted%ring.count]
        try read(slot)
        let length=min(batch,first+count-next)
        let command=queues[submitted%queues.count].makeCommandBuffer()!
        let enc=command.makeComputeCommandEncoder()!
        enc.setComputePipelineState(pipe)
        enc.setBuffer(sourceData,offset:0,index:0)
        enc.setBuffer(slot.buffer,offset:0,index:1)
        var params=Params(first:UInt32(next),length:UInt32(length))
        enc.setBytes(&params,length:MemoryLayout<Params>.stride,index:2)
        enc.dispatchThreads(MTLSize(width:length,height:1,depth:1),
                            threadsPerThreadgroup:MTLSize(width:64,height:1,depth:1))
        enc.endEncoding()
        slot.command=command
        slot.head=next;slot.length=length
        command.commit()
        next+=length
        submitted+=1
    }
    for slot in ring {try read(slot)}
    try check(cases==count,"Lost case")
    let result:[String:Any]=[
        "schema":"a28-r059-affine-exact-cover-metal-v1","first_id":first,"count":count,"device":device.name,
        "pid":ProcessInfo.processInfo.processIdentifier,"by_dimension_and_family":hist,"cases":cases,
        "overlapping_cases":overlap,"t_equals_m_minus_maximum":tightActivity,"t_equals_centered_real_rank":tightRank,
        "positional_checksum64":String(format:"%016llx",sum),"witness_codes":header,
        "gpu_command_intervals_epoch":intervals,"gpu_command_buffers":submitted,"gpu_command_seconds":time,
        "metal_wall_seconds":ProcessInfo.processInfo.systemUptime-clock,"all_checks_passed":true
    ]
    let out=try JSONSerialization.data(withJSONObject:result,options:[.prettyPrinted,.sortedKeys])
    try (out+Data([10])).write(to:URL(fileURLWithPath:args[2]))
    print("R59 Metal checked \(count) exact affine covers on \(device.name).")
}
do{try run()}catch{FileHandle.standardError.write(Data(("R59 Metal error: \(error)\n").utf8));exit(1)}
