import CryptoKit
import Dispatch
import Foundation
import Metal

let kernelSource = """
#include <metal_stdlib>
using namespace metal;

kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &stride [[buffer(2)]],
                       uint index [[thread_position_in_grid]]) {
    int own = input[index];
    int other = input[index ^ stride];
    output[index] = (index & stride) == 0 ? own + other : other - own;
}
"""

struct CheckFailure: Error {
    let reason: String
}

func check(_ condition: Bool, _ reason: String) throws {
    guard condition else { throw CheckFailure(reason: reason) }
}

func hash(_ data: Data) -> String {
    SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func parity(_ word: Int) -> Int {
    word.nonzeroBitCount & 1
}

func sign(_ bit: Int) -> Int32 {
    bit == 0 ? 1 : -1
}

func directWalsh(_ table: [Int32]) -> [Int32] {
    (0..<table.count).map { frequency in
        table.indices.reduce(into: Int32(0)) { sum, point in
            sum += table[point] * sign(parity(frequency & point))
        }
    }
}

final class MetalWalsh {
    let device: MTLDevice
    let queue: MTLCommandQueue
    let pipeline: MTLComputePipelineState
    var dispatches = 0
    var gpuSeconds: Double = 0
    var commandSpans: [(Double, Double)] = []

    init() throws {
        guard let device = MTLCreateSystemDefaultDevice(),
              let queue = device.makeCommandQueue()
        else { throw CheckFailure(reason: "Metal device or command queue unavailable") }
        let library = try device.makeLibrary(source: kernelSource, options: nil)
        guard let function = library.makeFunction(name: "fwht_stage") else {
            throw CheckFailure(reason: "Metal FWHT shader unavailable")
        }
        self.device = device
        self.queue = queue
        self.pipeline = try device.makeComputePipelineState(function: function)
    }

    func transform(_ input: [Int32], rowSize: Int) throws -> [Int32] {
        let bytes = input.count * MemoryLayout<Int32>.stride
        let first = input.withUnsafeBytes {
            device.makeBuffer(bytes: $0.baseAddress!, length: bytes, options: .storageModeShared)
        }
        guard var read = first,
              var write = device.makeBuffer(length: bytes, options: .storageModeShared),
              let command = queue.makeCommandBuffer()
        else { throw CheckFailure(reason: "Metal buffer allocation failed") }
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        var stride = 1
        while stride < rowSize {
            guard let encoder = command.makeComputeCommandEncoder() else {
                throw CheckFailure(reason: "Metal command encoder unavailable")
            }
            var step = UInt32(stride)
            encoder.setComputePipelineState(pipeline)
            encoder.setBuffer(read, offset: 0, index: 0)
            encoder.setBuffer(write, offset: 0, index: 1)
            encoder.setBytes(&step, length: MemoryLayout<UInt32>.stride, index: 2)
            encoder.dispatchThreads(
                MTLSize(width: input.count, height: 1, depth: 1),
                threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
            )
            encoder.endEncoding()
            swap(&read, &write)
            stride <<= 1
        }
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error { throw error }
        try check(command.status == .completed, "Metal command did not complete")
        dispatches += 1
        commandSpans.append((command.gpuStartTime, command.gpuEndTime))
        gpuSeconds += command.gpuEndTime - command.gpuStartTime
        let address = read.contents().bindMemory(to: Int32.self, capacity: input.count)
        return Array(UnsafeBufferPointer(start: address, count: input.count))
    }
}

struct RNG {
    var state: UInt64
    mutating func next() -> UInt64 {
        state &+= 0x9e3779b97f4a7c15
        var output = state
        output = (output ^ (output >> 30)) &* 0xbf58476d1ce4e5b9
        output = (output ^ (output >> 27)) &* 0x94d049bb133111eb
        return output ^ (output >> 31)
    }
}

func cpuFWHT(_ input: [Int32]) -> [Int32] {
    var result = input
    var width = 1
    while width < result.count {
        for start in stride(from: 0, to: result.count, by: 2 * width) {
            for offset in start..<(start + width) {
                let left = result[offset]
                let right = result[offset + width]
                result[offset] = left + right
                result[offset + width] = left - right
            }
        }
        width <<= 1
    }
    return result
}

let orbitShader = """
#include <metal_stdlib>
using namespace metal;

struct OrbitTask {
    uint d0, d1, d2;
    uint rank, pivot, offset;
};

kernel void orbit_histograms(device const int *truth [[buffer(0)]],
                             device const OrbitTask *tasks [[buffer(1)]],
                             device uint *histograms [[buffer(2)]],
                             uint3 workgroup [[threadgroup_position_in_grid]],
                             uint lane [[thread_index_in_threadgroup]]) {
    threadgroup atomic_uint hist[256];
    atomic_store_explicit(hist+lane,0,memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    OrbitTask job=tasks[workgroup.x];
    uint dirs[8]={0,job.d0,job.d1,job.d0^job.d1,job.d2,job.d0^job.d2,job.d1^job.d2,job.d0^job.d1^job.d2};
    uint vertices = 1u << job.rank;
    uint count = (1u<<14u)/vertices;
    for (uint i=lane; i<count; i+=256u) {
        // Choose the unique orbit representative having all pivot bits zero.
        uint x=0, original=0;
        for (uint j=0;j<14u;++j)
            if ((job.pivot & (1u<<j)) == 0u) {
                x |= ((i>>original)&1u)<<j;
                original += 1u;
            }
        uint word=0;
        for (uint j=0;j<vertices;++j) {
            uint a=dirs[j]>>1u;
            uint c=dirs[j]&1u;
            uint phase=(truth[job.offset+(x^a)]<0?1u:0u)^c;
            word|=phase<<j;
        }
        atomic_fetch_add_explicit(hist+word,1,memory_order_relaxed);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    histograms[workgroup.x*256u+lane]=atomic_load_explicit(hist+lane,memory_order_relaxed);
}
"""

struct Input: Decodable {
    struct Row: Decodable {
        let index: Int
        let support: [Int]
        let signs: [Int32]
        let period_dimension: Int
    }
    let source_sha256: String
    let rows: [Row]
}

struct OrbitTask {
    let d0: UInt32
    let d1: UInt32
    let d2: UInt32
    let rank: UInt32
    let pivot: UInt32
    let offset: UInt32
}

struct Code {
    let task: OrbitTask
    let basis: [Int]
    let elements: [Int]
    let words: [UInt64]
    let supportUnion: Int
    let row: Int
}

func makePivots(_ vectors: [Int], _ dimension: Int) -> [Int] {
    var rows = [Int](repeating:0,count:dimension)
    for x in vectors {
        var v = x
        while v != 0 {
            let top = Int.bitWidth - v.leadingZeroBitCount - 1
            if rows[top] == 0 { rows[top] = v; break }
            v ^= rows[top]
        }
    }
    return rows
}

func affineCertificate(_ points: [Int], _ dimension: Int) -> (rank: Int, anchor: Int, basis: [Int]) {
    let a = points[0]
    let c = makePivots(points.map { $0^a },dimension).filter { $0 != 0 }
    return (c.count,a,c)
}

func xorSpan(_ basis: [Int]) -> [Int] {
    var result = [0]
    for v in basis { result.append(contentsOf: result.map { $0^v }) }
    return result
}

func minimumWords(_ support: [Int]) throws -> ([Int],[UInt64],[Int:Int]) {
    var words: [UInt64] = [0]
    for j in 0..<15 {
        var column: UInt64 = 0
        for i in support.indices where (support[i]>>j)&1 != 0 { column |= UInt64(1)<<i }
        words.append(contentsOf: words.map { $0^column })
    }
    try check(words.count==32768 && Set(words).count==32768, "Affine phase-direction code is not injective")
    let minWeight = words.dropFirst().map { $0.nonzeroBitCount }.min()!
    try check(minWeight==8, "Unexpected minimum")
    let directions = words.indices.filter { words[$0].nonzeroBitCount==minWeight }
    var weights: [Int:Int] = [:]
    for word in words { weights[word.nonzeroBitCount,default:0] += 1 }
    return (directions,words,weights)
}

func codeCatalog(_ min: [Int], _ words: [UInt64], _ index: Int, _ offset: Int) throws -> [Code] {
    let n=min.count
    var rank2=Set<[Int]>()
    var rank3=Set<[Int]>()
    for i in 0..<n {
        for j in (i+1)..<n {
            let a=min[i], b=min[j]
            rank2.insert([a,b,a^b].sorted())
            for k in (j+1)..<n {
                let c=min[k]
                if c != a^b { rank3.insert([a,b,c,a^b,a^c,b^c,a^b^c].sorted()) }
            }
        }
    }
    let codes = min.map { [$0] } + rank2.sorted { $0.lexicographicallyPrecedes($1) }
                                  + rank3.sorted { $0.lexicographicallyPrecedes($1) }
    try check(min.count==24 && rank2.count==244 && rank3.count==1372,"Wrong exhaustive code counts")
    var result:[Code]=[]
    for members in codes {
        var basis:[Int]=[]
        var space=[0]
        for d in members where words[d].nonzeroBitCount==8 {
            if !space.contains(d) {
                basis.append(d)
                space.append(contentsOf:space.map { $0^d })
            }
        }
        try check(Array(space.dropFirst()).sorted()==members,"Code lacks a minimum-word basis")
        let physical=makePivots(basis.map { $0>>1 },14)
        try check(physical.filter { $0 != 0 }.count==basis.count,"Phase-only nonzero direction in code")
        let p=physical.indices.reduce(0) { $0 | (physical[$1] != 0 ? 1<<$1 : 0) }
        let unionWeight=members.reduce(UInt64(0)) { $0|words[$1] }.nonzeroBitCount
        let lower=[0,8,12,16][basis.count]
        try check(unionWeight>=lower,"Spectral support union below dim 1/2/3 bound")
        let padded=basis+[0,0,0]
        let task=OrbitTask(d0:UInt32(padded[0]),d1:UInt32(padded[1]),d2:UInt32(padded[2]),
                           rank:UInt32(basis.count),pivot:UInt32(p),offset:UInt32(offset))
        result.append(Code(task:task,basis:basis,elements:space,words:space.map{words[$0]},
                           supportUnion:unionWeight,row:index))
    }
    return result
}

func directHistogram(_ code: Code, _ truth: [Int32]) -> [Int] {
    let p=Int(code.task.pivot)
    let offset=Int(code.task.offset)
    var hist=[Int](repeating:0,count:256)
    for x in 0..<16384 where p&x==0 {
        let mask=code.elements.enumerated().reduce(0) { w, item in
            let d=item.element
            let bit=(truth[offset+(x^(d>>1))]<0 ? 1:0)^(d&1)
            return w|(bit<<item.offset)
        }
        hist[mask]+=1
    }
    return hist
}

struct Pattern {
    let nonconstant: Int
    let energyNumerator: Int
    let derivativeCounts: [Int]
    let allBasisCoincident: Int
    let allDirectionsCoincident: Int
}

func patterns(_ rank: Int) -> [Pattern] {
    let n=1<<rank
    return (0..<(1<<n)).map { mask in
        let w=mask.nonzeroBitCount
        func crossed(_ x: Int,_ d: Int) -> Bool {
            ((mask>>x)^(mask>>(x^d)))&1 != 0
        }
        let derivativeCounts=(1..<n).map { d in (0..<n).filter { x in crossed(x,d) }.count }
        let independent=(0..<rank).map { 1<<$0 }
        let common=(0..<n).filter { x in independent.allSatisfy { crossed(x,$0) } }.count
        let entire=(0..<n).filter { x in (1..<n).allSatisfy { crossed(x,$0) } }.count
        return Pattern(nonconstant:w>0 && w<n ? 1:0,
                       energyNumerator:4*w*(n-w),derivativeCounts:derivativeCounts,
                       allBasisCoincident:common,allDirectionsCoincident:entire)
    }
}

struct Distribution {
    var codeCount=0
    var smallestActivity=16385
    var largestActivity=0
    var minCommon=16385
    var maxCommon=0
    var minAll=16385
    var maxAll=0
    var activityEqualities=0
    var zeroCommon=0
    var zeroAll=0

    mutating func add(_ activity: Int,_ common: Int,_ all: Int,_ equality: Bool) {
        codeCount+=1
        smallestActivity=min(smallestActivity,activity); largestActivity=max(largestActivity,activity)
        minCommon=min(minCommon,common); maxCommon=max(maxCommon,common)
        minAll=min(minAll,all); maxAll=max(maxAll,all)
        if equality { activityEqualities+=1 }
        if common==0 { zeroCommon+=1 }
        if all==0 { zeroAll+=1 }
    }

    var json: [String:Any] {
        [
            "subcodes":codeCount,"active_point_range":[smallestActivity,largestActivity],
            "common_independent_minimal_derivative_point_range":[minCommon,maxCommon],
            "common_all_nonzero_derivative_point_range":[minAll,maxAll],
            "active_density_equals_spectral_union_mass":activityEqualities,
            "zero_common_independent_minimal_derivative_cases":zeroCommon,
            "zero_common_all_nonzero_derivative_cases":zeroAll,
        ]
    }
}

func run() throws {
    let args=CommandLine.arguments
    try check(args.count==5,"usage: R53FlatOrbitMetal r11.json output.json firstRow endRow")
    let inputBytes=try Data(contentsOf:URL(fileURLWithPath:args[1]))
    let input=try JSONDecoder().decode(Input.self,from:inputBytes)
    guard let first=Int(args[3]),let end=Int(args[4]),0<=first,first<end,end<=input.rows.count else {
        throw CheckFailure(reason:"Invalid shard range")
    }
    let rows=Array(input.rows[first..<end])
    let time=ProcessInfo.processInfo.systemUptime
    let metal=try MetalWalsh()
    var sparse=[Int32](repeating:0,count:16384*rows.count)
    for (i,row) in rows.enumerated() {
        try check(row.support.count==64 && row.signs.count==64,"Invalid flat input")
        try check(Set(row.support).count==64 && row.support.allSatisfy { $0&1==1 },"Missing phase coordinate")
        try check(row.signs.allSatisfy { abs($0)==1 },"Nonflat seed input")
        for (j,frequency) in row.support.enumerated() { sparse[i*16384+(frequency>>1)]=row.signs[j] }
    }
    let truth=try metal.transform(sparse,rowSize:16384)
    for i in rows.indices {
        let r=(16384*i)..<(16384*(i+1))
        try check(cpuFWHT(Array(sparse[r]))==Array(truth[r]) && truth[r].allSatisfy { abs($0)==8 },
                  "CPU/Metal inverse FWHT mismatch or non-Boolean truth")
    }
    let orbitLib=try metal.device.makeLibrary(source:orbitShader,options:nil)
    let pipeline=try metal.device.makeComputePipelineState(function:orbitLib.makeFunction(name:"orbit_histograms")!)
    let truthBuffer=truth.withUnsafeBytes {
        metal.device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
    }

    var allCodes:[Code]=[]
    var rowSummary:[[String:Any]]=[]
    var phaseHistogram=[0,0]
    var lifts:[String:[String:Any]]=[:]
    for (i,row) in rows.enumerated() {
        let offset=i*16384
        let (min,words,weights)=try minimumWords(row.support)
        for d in min {
            let wordPoints=words[d].nonzeroBitCount
            let activeIndices=row.support.indices.filter { (words[d]>>$0)&1==1 }
            let frequency=activeIndices.map { row.support[$0]>>1 }
            let aff=affineCertificate(frequency,14)
            let intrinsic=affineCertificate(activeIndices,6)
            try check(wordPoints==8 && aff.rank==3 && intrinsic.rank==3,"Minimum spectral word does not give a flat")
            let support=(0..<16384).filter { x in
                truth[offset+x] != sign(d&1)*truth[offset+(x^(d>>1))]
            }
            try check(support.count==2048,"Minimum derivative physical density is wrong")
            let physical=affineCertificate(support,14)
            try check(physical.rank==11,"Minimum derivative active domain is not an affine flat")
            let supportForOtherPhase=(0..<16384).filter { x in
                truth[offset+x] == sign(d&1)*truth[offset+(x^(d>>1))]
            }
            try check(supportForOtherPhase.count==14336 && words[d^1].nonzeroBitCount==56,
                      "Phase-complement check failed")
            phaseHistogram[d&1]+=1
            let key=d&1==0 ? "phase_zero" : "phase_one"
            if lifts[key]==nil {
                lifts[key]=[
                    "row":row.index,"phase_direction":d,"physical_shift":d>>1,"phase":d&1,
                    "spectral_flat_anchor":aff.anchor,"spectral_flat_basis":aff.basis,
                    "physical_active_flat_anchor":physical.anchor,"physical_active_flat_basis":physical.basis,
                    "spectral_word_points":wordPoints,"physical_active_size":support.count
                ]
            }
        }
        let codes=try codeCatalog(min,words,row.index,offset)
        allCodes.append(contentsOf:codes)
        rowSummary.append([
            "row_index":row.index,"period_dimension":row.period_dimension,
            "minimum_phase_directions":min,
            "full_direction_weight_histogram":Dictionary(uniqueKeysWithValues:weights.map { (String($0.key),$0.value) }),
            "row_codes_dim1_dim2_dim3":(1...3).map { rank in codes.filter { $0.basis.count==rank }.count },
        ])
    }
    let inputGPU=allCodes.map { $0.task }
    try check(MemoryLayout<OrbitTask>.stride==24,"Unexpected GPU task layout")
    let jobs=inputGPU.withUnsafeBytes {
        metal.device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
    }
    let hist=metal.device.makeBuffer(length:allCodes.count*256*4,options:.storageModeShared)!
    let command=metal.queue.makeCommandBuffer()!
    let encoder=command.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(truthBuffer,offset:0,index:0)
    encoder.setBuffer(jobs,offset:0,index:1)
    encoder.setBuffer(hist,offset:0,index:2)
    encoder.dispatchThreadgroups(MTLSize(width:allCodes.count,height:1,depth:1),
                                  threadsPerThreadgroup:MTLSize(width:256,height:1,depth:1))
    encoder.endEncoding()
    command.commit()
    command.waitUntilCompleted()
    if let error=command.error { throw error }
    try check(command.status == .completed,"Metal orbit evaluation failed")
    let bins=hist.contents().bindMemory(to:UInt32.self,capacity:allCodes.count*256)

    let shape=[[Pattern]()] + (1...3).map(patterns)
    var stats:[String:Distribution]=[:]
    var witnesses:[String:[String:Any]]=[:]
    var subsetCounts=[0,0,0,0]
    var directCounts=[0,0,0,0]
    var orbitCounts=[0,0,0,0]
    var truthAndPhaseCounts=[0,0,0,0]
    var derivativeIdentities=0
    for (i,code) in allCodes.enumerated() {
        let r=code.basis.count, m=1<<r, n=16384/m
        let start=256*i
        let patterns=shape[r]
        let key="dim\(r)_union\(code.supportUnion)"
        var activity=0, common=0, commonAll=0, energy=0, sum=0
        var derivativeCounts=[Int](repeating:0,count:m-1)
        for (pattern,p) in patterns.enumerated() {
            let count=Int(bins[start+pattern])
            if count==0 { continue }
            sum+=count
            activity+=count*p.nonconstant*m
            common+=count*p.allBasisCoincident
            commonAll+=count*p.allDirectionsCoincident
            energy+=count*p.energyNumerator
            for j in 0..<(m-1) { derivativeCounts[j]+=count*p.derivativeCounts[j] }
        }
        try check(sum==n,"Incomplete orbit histogram")
        try check(64*energy==code.supportUnion*n*m*m,"Orbit mean identity failed")
        try check(activity*64>=16384*code.supportUnion,"Physical activity below spectral union density")
        for j in 0..<(m-1) {
            try check(derivativeCounts[j]*64==code.words[j+1].nonzeroBitCount*16384,
                      "Actual derivative activity and spectral fraction differ")
            derivativeIdentities+=1
        }
        let equality=activity*64==16384*code.supportUnion
        stats[key,default:Distribution()].add(activity,common,commonAll,equality)
        subsetCounts[r]+=1
        orbitCounts[r]+=n
        truthAndPhaseCounts[r]+=n*m
        if r==1 || i%89==0 || witnesses[key]==nil {
            let cpu=directHistogram(code,truth)
            try check(cpu.indices.allSatisfy {cpu[$0]==Int(bins[start+$0])}, "CPU/Metal orbit histograms differ")
            directCounts[r]+=1
        }
        if witnesses[key]==nil {
            witnesses[key]=[
                "row":code.row,"minimum_phase_direction_basis":code.basis,
                "nonzero_directions":Array(code.elements.dropFirst()),
                "weights":Array(code.words.dropFirst()).map { $0.nonzeroBitCount },
                "spectral_support_union":code.supportUnion,
                "physical_orbit_count":n,"twisted_truth_pattern_histogram":Dictionary(
                    uniqueKeysWithValues:patterns.indices.filter { bins[start+$0]>0 }
                    .map { (String($0),Int(bins[start+$0])) }),
                "active_physical_points":activity,
                "points_where_all_chosen_minimum_derivatives_act":common,
                "points_where_all_nonzero_derivatives_act":commonAll,
            ]
        }
    }
    let output:[String:Any]=[
        "schema":"a28-r053-flat-spectra-minimum-orbit-audit-v1",
        "input_sha256":sha256Data(inputBytes),"archive_input_sha256":input.source_sha256,
        "row_interval_half_open":[first,end],"device":metal.device.name,
        "record_count":rows.count,"whole_boolean_truth_cpu_metal_values":truth.count,
        "enumerated_phase_translation_words_including_zero":rows.count*32768,
        "minimum_directions_by_phase_bit":phaseHistogram,
        "minimum_spectral_flats_and_physical_flats_checked":phaseHistogram.reduce(0,+),
        "minimum_phase_witnesses":lifts,
        "complete_subcode_counts_dim1_dim2_dim3":Array(subsetCounts.dropFirst()),
        "metal_orbits_dim1_dim2_dim3":Array(orbitCounts.dropFirst()),
        "metal_twisted_truth_values_dim1_dim2_dim3":Array(truthAndPhaseCounts.dropFirst()),
        "direct_cpu_orbit_histograms_crosschecked_dim1_dim2_dim3":Array(directCounts.dropFirst()),
        "derivative_density_equalities_all_nonzero_code_directions":derivativeIdentities,
        "orbit_statistics_by_dimension_and_spectral_union":stats.mapValues { $0.json },
        "orbit_exact_witnesses":witnesses,"rows":rowSummary,
        "gpu_fwht_command_seconds":metal.gpuSeconds,
        "gpu_orbit_command_seconds":command.gpuEndTime-command.gpuStartTime,
        "verification_wall_seconds":ProcessInfo.processInfo.systemUptime-time,
        "all_checks_passed":true,
    ]
    let rendered=try JSONSerialization.data(withJSONObject:output,options:[.prettyPrinted,.sortedKeys])
    try (rendered+Data([10])).write(to:URL(fileURLWithPath:args[2]))
    print("R53 verified \(rows.count) seed rows and \(allCodes.count) distinct subcodes on \(metal.device.name).")
}

func sha256Data(_ bytes: Data) -> String {
    SHA256.hash(data: bytes).map { String(format:"%02x",$0) }.joined()
}

do { try run() }
catch {
    FileHandle.standardError.write(Data(("R53 verification failed: \(error)\n").utf8))
    exit(1)
}
