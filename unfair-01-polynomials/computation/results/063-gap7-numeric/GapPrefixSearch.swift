import Foundation
import Metal

// Independent projected Adam trajectories for the fixed low convolution prefix.
let source = """
#include <metal_stdlib>
using namespace metal;

struct Config {
    uint m, n, count, start, steps;
    float learning;
};

inline bool aLocked(uint i, uint m, uint n) {
    return i == m-1 || i == m-2 || i == m-5 ||
           i == n-1 || i == n-2 || i == n-5;
}

inline bool bLocked(uint i, uint m, uint n) {
    return i == m-1 || i == m-3 || i == m-4 ||
           i == n-1 || i == n-3 || i == n-4;
}

kernel void search(device float* coefficients [[buffer(0)]],
                   device float* first [[buffer(1)]],
                   device float* second [[buffer(2)]],
                   device float* residuals [[buffer(3)]],
                   constant Config& config [[buffer(4)]],
                   uint id [[thread_position_in_grid]]) {
    if (id >= config.count) return;
    const uint m = config.m;
    const uint n = config.n;
    const uint stride = m+n+2;
    const uint offset = id*stride;
    float a[64], b[64], r[127], gradient[128];
    for (uint i=0; i<=m; ++i) a[i] = coefficients[offset+i];
    for (uint j=0; j<=n; ++j) b[j] = coefficients[offset+m+1+j];
    for (uint iter=0; iter<config.steps; ++iter) {
        for (uint k=8; k<=m+n; ++k) {
            float c = 0;
            const uint i0 = k>n ? k-n : 0;
            const uint i1 = min(m,k);
            for (uint i=i0; i<=i1; ++i) c += a[i]*b[k-i];
            const bool forced = (k <= 11 || k == m || k == n);
            const float chosen = forced || c >= 0.5f ? 1.0f : 0.0f;
            r[k] = (forced ? 3.0f : 1.0f)*(c-chosen);
        }
        for (uint i=8; i<m; ++i) {
            float g = 0;
            for (uint j=0; j<=n; ++j) g += b[j]*r[i+j];
            gradient[i] = g;
        }
        for (uint j=8; j<n; ++j) {
            float g = 0;
            for (uint i=0; i<=m; ++i) g += a[i]*r[i+j];
            gradient[m+1+j] = g;
        }
        const float iteration = float(config.start+iter+1);
        const float rate = config.learning*sqrt(1.0f-pow(0.999f,iteration))/(1.0f-pow(0.9f,iteration));
        for (uint i=8; i<m; ++i) {
            if (aLocked(i,m,n)) continue;
            const uint pos = offset+i;
            const float g = gradient[i];
            const float f = 0.9f*first[pos]+0.1f*g;
            const float s = 0.999f*second[pos]+0.001f*g*g;
            first[pos]=f; second[pos]=s;
            a[i]=clamp(a[i]-rate*f/(sqrt(s)+1.0e-8f),0.0f,1.0f);
        }
        for (uint j=8; j<n; ++j) {
            if (bLocked(j,m,n)) continue;
            const uint pos = offset+m+1+j;
            const float g = gradient[m+1+j];
            const float f = 0.9f*first[pos]+0.1f*g;
            const float s = 0.999f*second[pos]+0.001f*g*g;
            first[pos]=f; second[pos]=s;
            b[j]=clamp(b[j]-rate*f/(sqrt(s)+1.0e-8f),0.0f,1.0f);
        }
    }
    for (uint i=8; i<m; ++i) coefficients[offset+i]=a[i];
    for (uint j=8; j<n; ++j) coefficients[offset+m+1+j]=b[j];
    float maxError=0;
    for (uint k=8; k<=m+n; ++k) {
        float c=0;
        const uint i0=k>n ? k-n : 0;
        const uint i1=min(m,k);
        for (uint i=i0; i<=i1; ++i) c += a[i]*b[k-i];
        const bool forced = (k <= 11 || k == m || k == n);
        const float chosen = forced || c >= 0.5f ? 1.0f : 0.0f;
        maxError=max(maxError,abs(c-chosen));
    }
    residuals[id]=maxError;
}
"""

struct Config {
    var m, n, count, start, steps: UInt32
    var learning: Float
}
struct SplitMix {
    var state: UInt64
    mutating func sample() -> Float {
        state &+= 0x9e3779b97f4a7c15
        var z = state
        z = (z ^ (z >> 30)) &* 0xbf58476d1ce4e5b9
        z = (z ^ (z >> 27)) &* 0x94d049bb133111eb
        z ^= z >> 31
        return Float(z >> 40) / 16777216
    }
}
struct Case: Codable {
    let A, B: [Float]
    let residual: Float
}
struct Record: Codable {
    let step: Int
    let minimum: Float
    let seconds: Double
}
struct Result: Codable {
    let status, device: String
    let m, n, trajectories, steps: Int
    let seed: UInt64
    let completed: Int
    let history: [Record]
    let cases: [Case]
}

let args=CommandLine.arguments
guard args.count==7, let m=Int(args[1]), let n=Int(args[2]), let count=Int(args[3]),
      let total=Int(args[4]), let seed=UInt64(args[5]),
      (15..<64).contains(m), (15..<64).contains(n), m != n, (16...65536).contains(count), total>0
else { fatalError("usage: m n trajectories steps seed new-output-json") }
let out=URL(fileURLWithPath:args[6])
guard !FileManager.default.fileExists(atPath:out.path) else { fatalError("existing output") }
let stride=m+n+2
let bytes=count*stride*MemoryLayout<Float>.stride
let device=MTLCreateSystemDefaultDevice()!
let queue=device.makeCommandQueue()!
let library=try device.makeLibrary(source:source,options:nil)
let pipeline=try device.makeComputePipelineState(function:library.makeFunction(name:"search")!)
let coeff=device.makeBuffer(length:bytes,options:.storageModeShared)!
let first=device.makeBuffer(length:bytes,options:.storageModeShared)!
let second=device.makeBuffer(length:bytes,options:.storageModeShared)!
let scores=device.makeBuffer(length:count*MemoryLayout<Float>.stride,options:.storageModeShared)!
memset(first.contents(),0,bytes)
memset(second.contents(),0,bytes)
let state=coeff.contents().bindMemory(to:Float.self,capacity:count*stride)
let score=scores.contents().bindMemory(to:Float.self,capacity:count)
var random=SplitMix(state:seed)
let t=0.719429655179540
let b2=1-t+t*t
let a3=1-t*b2
let a4=1-(1-t)*a3
let a0=[1,t,0,a3,a4,0,0,0].map(Float.init)
let b0=[1,1-t,b2,0,0,t*t*b2,0,0].map(Float.init)
let aZero=Set([m-1,m-2,m-5,n-1,n-2,n-5].filter { $0 <= m })
let bZero=Set([m-1,m-3,m-4,n-1,n-3,n-4].filter { $0 <= n })
guard aZero.allSatisfy({ $0 >= 8 && $0 < m }),
      bZero.allSatisfy({ $0 >= 8 && $0 < n }) else { fatalError("forced zero conflicts with fixed input") }
for id in 0..<count {
    let base=id*stride
    for i in 0...m {
        state[base+i]=i<8 ? a0[i] : i==m ? 1 : aZero.contains(i) ? 0 : 0.22*random.sample()
    }
    for j in 0...n {
        state[base+m+1+j]=j<8 ? b0[j] : j==n ? 1 : bZero.contains(j) ? 0 : 0.22*random.sample()
    }
}
let encoder=JSONEncoder()
encoder.outputFormatting=[.prettyPrinted,.sortedKeys]
let begin=DispatchTime.now().uptimeNanoseconds
func elapsed() -> Double { Double(DispatchTime.now().uptimeNanoseconds-begin)/1_000_000_000 }
var history=[Record]()
var best=[Case]()
var step=0
let block=16
while step<total && elapsed()<14400 {
    let amount=min(total-step,block)
    var config=Config(m:UInt32(m),n:UInt32(n),count:UInt32(count),start:UInt32(step),steps:UInt32(amount),learning:0.026)
    let command=queue.makeCommandBuffer()!
    let encode=command.makeComputeCommandEncoder()!
    encode.setComputePipelineState(pipeline)
    for (index,b) in [coeff,first,second,scores].enumerated() { encode.setBuffer(b,offset:0,index:index) }
    encode.setBytes(&config,length:MemoryLayout<Config>.stride,index:4)
    encode.dispatchThreads(MTLSize(width:count,height:1,depth:1),
                           threadsPerThreadgroup:MTLSize(width:min(128,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
    encode.endEncoding()
    command.commit()
    command.waitUntilCompleted()
    guard command.status == .completed else { fatalError("Metal command failed") }
    step+=amount
    if step==block || step%512==0 || step==total {
        let indices=(0..<count).sorted{ score[$0]<score[$1] }.prefix(16)
        best=indices.map { id in
            let address=state.advanced(by:id*stride)
            return Case(A:Array(UnsafeBufferPointer(start:address,count:m+1)),
                        B:Array(UnsafeBufferPointer(start:address.advanced(by:m+1),count:n+1)),
                        residual:score[id])
        }
        guard best.allSatisfy({ $0.residual.isFinite }) else { fatalError("nonfinite result") }
        let event=Record(step:step,minimum:best.first!.residual,seconds:elapsed())
        history.append(event)
        let checkpoint=Result(status:"checkpoint",device:device.name,m:m,n:n,trajectories:count,
                              steps:total,seed:seed,completed:step,history:history,cases:best)
        try encoder.encode(checkpoint).write(to:out.appendingPathExtension("checkpoint"),options:.atomic)
        print(String(decoding:try encoder.encode(event),as:UTF8.self))
        fflush(stdout)
    }
}
let result=Result(status:step==total ? "finished_numeric_search" : "partial_time",device:device.name,
                  m:m,n:n,trajectories:count,steps:total,seed:seed,completed:step,history:history,cases:best)
try encoder.encode(result).write(to:out,options:.atomic)
