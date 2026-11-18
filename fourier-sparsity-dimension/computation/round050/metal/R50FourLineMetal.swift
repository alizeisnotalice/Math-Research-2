import CryptoKit
import Foundation
import Metal

let shader = """
#include <metal_stdlib>
using namespace metal;
kernel void joined_walsh(
    device const int *bank [[buffer(0)]],
    device const uint *address_to_line [[buffer(1)]],
    device int *spectra [[buffer(2)]],
    constant uint &first [[buffer(3)]],
    uint3 group [[threadgroup_position_in_grid]],
    uint lane [[thread_index_in_threadgroup]]) {
    threadgroup int row[128];
    uint id = first + group.x;
    for (uint i=lane; i<128u; i+=64u) {
        uint loc = address_to_line[i & 7u];
        uint line = loc >> 1;
        uint t = loc & 1u;
        uint choice = (id >> (5u*line)) & 31u;
        row[i] = bank[(line*32u+choice)*32u + t*16u + (i>>3)];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint h=1; h<128u; h*=2u) {
        uint left = ((lane & ~(h-1u)) << 1u) | (lane & (h-1u));
        uint right = left + h;
        int a = row[left], b = row[right];
        row[left] = a+b;
        row[right] = a-b;
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    spectra[128u*group.x+lane] = row[lane];
    spectra[128u*group.x+lane+64u] = row[lane+64u];
}
"""

struct CheckFailure: Error { let explanation: String }
func check(_ condition: Bool, _ description: String) throws {
    guard condition else { throw CheckFailure(explanation: description) }
}

func parity(_ x: Int) -> Int { x.nonzeroBitCount & 1 }
func sign(_ b: Int) -> Int32 { b == 0 ? 1 : -1 }

struct Phase {
    let label: Int
    let p: Int
    let q: Int
    let constant: Int
    let inputShift: Int

    func evaluate(_ y: Int) -> Int32 {
        let t = y ^ inputShift
        let term = parity(p&t) & parity(q&t)
        return sign(term ^ parity(label&y) ^ constant)
    }
}

struct Choice {
    let truth: [Int32]
    let lifted: [Int32]
    let dataAddressCounts: [Int]
    let intrinsicPreimages: [Int]
}

let linePoints = [[0,1], [2,4], [3,7], [5,6]]
let ambientSize = 128

func walsh(_ truth: [Int32]) -> [Int32] {
    (0..<truth.count).map { alpha in
        truth.indices.reduce(into: Int32(0)) { sum, x in
            sum += sign(parity(alpha&x)) * truth[x]
        }
    }
}

func makeBank() throws -> [[Choice]] {
    var output: [[Choice]] = []
    for line in 0..<4 {
        let p = [1,1,4,2][line]
        let q = [2,4,8,8][line]
        let outside = [4,2,1,1][line]
        let endpoints = linePoints[line]
        let direction = endpoints[0] ^ endpoints[1]
        var choices: [Choice] = []
        for id in 0..<32 {
            let mode = id / 4
            let variant = id % 4
            let offset = line ^ (variant << 2)
            func character(_ label: Int, _ c: Int = 0) -> Phase {
                Phase(label: label,p: 0,q: 0,constant: c,inputShift: 0)
            }
            func bent(_ label: Int, _ quadP: Int, _ quadQ: Int,
                      _ c: Int = 0, _ shift: Int = 0) -> Phase {
                Phase(label: label,p: quadP,q: quadQ,constant: c,inputShift: shift)
            }
            let pair: [Phase]
            switch mode {
            case 0: pair = [character(offset), character(offset)]
            case 1: pair = [character(offset), character(offset ^ outside)]
            case 2: pair = [bent(offset,p,q), bent(offset ^ p,p,q,1)]
            case 3: pair = [bent(offset,p,q), bent(offset ^ outside,p,q)]
            case 4: pair = [bent(offset,p,q), bent(offset,p ^ outside,q,0,1)]
            case 5: pair = [character(offset), bent(offset,p,q)]
            case 6: pair = [bent(offset,p,q), character(offset ^ outside)]
            default: pair = [bent(offset ^ q,p ^ outside,q),
                             bent(offset ^ p,p,q,0,1)]
            }
            let byEndpoint = pair.map { phase in (0..<16).map(phase.evaluate(_:)) }
            let spectra = byEndpoint.map(walsh)
            for i in 0..<2 {
                let expected = pair[i].p == 0 ? 1 : 4
                try check(spectra[i].filter { $0 != 0 }.count == expected,
                          "Choice endpoint is not the advertised character or four-point bent spectrum")
            }
            var lifted: [Int32] = []
            for z in 0..<ambientSize {
                let alpha = z & 7
                let beta = z >> 3
                lifted.append(sign(parity(alpha & endpoints[0])) *
                              (spectra[0][beta] + sign(parity(alpha & direction)) * spectra[1][beta]))
            }
            let preimageCounts = (0..<16).map { beta -> Int in
                let a = spectra[0][beta], b = spectra[1][beta]
                return (a+b != 0 ? 1 : 0) + (a-b != 0 ? 1 : 0)
            }
            let addresses = (0..<16).map { beta in spectra.filter { $0[beta] != 0 }.count }
            var zeroExtended = [Int32](repeating: 0, count: ambientSize)
            for y in 0..<16 {
                for t in 0..<2 { zeroExtended[endpoints[t] + 8*y] = byEndpoint[t][y] }
            }
            try check(walsh(zeroExtended) == lifted, "Lifted bank spectrum differs from direct transform")
            for beta in 0..<16 {
                let count = lifted[(8*beta)..<(8*beta+8)].filter { $0 != 0 }.count
                try check(count == 4*preimageCounts[beta], "Incorrect local coset multiplicity")
            }
            choices.append(Choice(truth: byEndpoint.flatMap { $0 }, lifted: lifted,
                                  dataAddressCounts: addresses,intrinsicPreimages: preimageCounts))
        }
        try check(Set(choices.map { $0.truth }).count == 32, "Duplicate local choice")
        output.append(choices)
    }
    return output
}

func sha(_ bytes: Data) -> String {
    SHA256.hash(data: bytes).map { String(format: "%02x", $0) }.joined()
}

final class Slot {
    let buffer: MTLBuffer
    var command: MTLCommandBuffer?
    var first = 0
    var count = 0

    init(_ device: MTLDevice, _ batch: Int) {
        buffer = device.makeBuffer(length: batch*ambientSize*MemoryLayout<Int32>.stride,options:.storageModeShared)!
    }
}

struct Totals {
    var cases: Int = 0
    var bands: Int = 0
    var activeBands: Int = 0
    var privateBands: Int = 0
    var privateTwoPreimages: Int = 0
    var collidingBands: Int = 0
    var cancelingBands: Int = 0
    var cancellations: Int = 0
    var bandMassEqualities: Int = 0
    var actualAddressMassEqualities: Int = 0
    var globalMassEqualities: Int = 0
    var fullDirectCases: Int = 0
    var spectrumCoefficients: Int = 0
    var gpuSeconds: Double = 0
    var gpuWindows: [[Double]] = []
    var histogramBySparsity = [Int](repeating: 0,count: 129)
    var selectedWitnesses: [String:[String:Any]] = [:]
    var directWitnesses: [[String:Any]] = []
}

func run() throws {
    let args = CommandLine.arguments
    try check(args.count == 5, "usage: R50FourLineMetal output.json first count batch")
    guard let first = Int(args[2]), let count = Int(args[3]), let batch = Int(args[4]),
          first >= 0, count > 0, first + count <= 1 << 20, batch > 0 && batch <= 32768 else {
        throw CheckFailure(explanation: "Invalid range or batch")
    }
    let bank = try makeBank()
    let entries = bank.flatMap { $0 }
    let words = entries.flatMap { $0.truth }
    var addressToLine = [UInt32](repeating: 99, count: 8)
    for i in 0..<4 {
        for t in 0..<2 { addressToLine[linePoints[i][t]] = UInt32(2*i+t) }
    }
    try check(Set(addressToLine) == Set((0..<8).map(UInt32.init)), "Addresses are not a partition")
    guard let device = MTLCreateSystemDefaultDevice() else { throw CheckFailure(explanation: "Metal unavailable") }
    let library = try device.makeLibrary(source: shader,options:nil)
    let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name:"joined_walsh")!)
    let queues = (0..<2).map { _ in device.makeCommandQueue()! }
    let slots = (0..<3).map { _ in Slot(device,batch) }
    let bankGPU = words.withUnsafeBytes {
        device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
    }
    let geometryGPU = addressToLine.withUnsafeBytes {
        device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!
    }
    let uptimeStart = ProcessInfo.processInfo.systemUptime
    let epochOffset = Date().timeIntervalSince1970 - uptimeStart
    var totals = Totals()

    func collect(_ slot: Slot) throws {
        guard let command = slot.command else { return }
        command.waitUntilCompleted()
        if let error = command.error { throw error }
        try check(command.status == .completed, "Metal buffer failed")
        totals.gpuSeconds += command.gpuEndTime-command.gpuStartTime
        totals.gpuWindows.append([command.gpuStartTime+epochOffset, command.gpuEndTime+epochOffset])
        let gpu = slot.buffer.contents().bindMemory(to:Int32.self,capacity:slot.count*ambientSize)
        for row in 0..<slot.count {
            let id = slot.first+row
            let c0 = entries[id&31]
            let c1 = entries[32+((id>>5)&31)]
            let c2 = entries[64+((id>>10)&31)]
            let c3 = entries[96+((id>>15)&31)]
            let rowStart = row*ambientSize
            var globalS = 0
            var densityLowerNumerator = 0
            var mismatch = false
            var energy: Int = 0

            for beta in 0..<16 {
                let r0 = c0.intrinsicPreimages[beta], r1 = c1.intrinsicPreimages[beta]
                let r2 = c2.intrinsicPreimages[beta], r3 = c3.intrinsicPreimages[beta]
                let lineCount = (r0>0 ? 1:0)+(r1>0 ? 1:0)+(r2>0 ? 1:0)+(r3>0 ? 1:0)
                let spatial = c0.dataAddressCounts[beta]+c1.dataAddressCounts[beta]
                            + c2.dataAddressCounts[beta]+c3.dataAddressCounts[beta]
                var bandS = 0
                var union = 0
                for alpha in 0..<8 {
                    let z = 8*beta+alpha
                    let a = c0.lifted[z], b = c1.lifted[z], c = c2.lifted[z], d = c3.lifted[z]
                    let value = a+b+c+d
                    if gpu[rowStart+z] != value { mismatch = true }
                    if value != 0 { bandS += 1 }
                    if a != 0 || b != 0 || c != 0 || d != 0 { union += 1 }
                    energy += Int(value)*Int(value)
                }
                totals.bands += 1
                globalS += bandS
                if lineCount == 0 {
                    try check(bandS == 0 && spatial == 0,"Empty incidence has nonempty true band")
                    continue
                }
                totals.activeBands += 1
                try check(bandS*lineCount >= 4 && bandS*spatial >= 8,
                          "Density lower bound fails")
                densityLowerNumerator += 48/lineCount
                if bandS*lineCount == 4 { totals.bandMassEqualities += 1 }
                if bandS*spatial == 8 { totals.actualAddressMassEqualities += 1 }
                let countBefore = 4*(r0+r1+r2+r3)
                if countBefore > union { totals.collidingBands += 1 }
                if union > bandS {
                    totals.cancelingBands += 1
                    totals.cancellations += union-bandS
                }
                if lineCount == 1 {
                    totals.privateBands += 1
                    try check(bandS == countBefore,"Private band lift is not exact")
                    if r0+r1+r2+r3 == 2 { totals.privateTwoPreimages += 1 }
                }
                let kind: String?
                if lineCount == 1 && bandS == 8 { kind = "private_with_two_local_preimages" }
                else if lineCount > 1 && union > bandS { kind = "overlapping_bands_with_cancellation" }
                else { kind = nil }
                if let key = kind, totals.selectedWitnesses[key] == nil {
                    totals.selectedWitnesses[key] = [
                        "case_id": id,"data_frequency": beta,"active_lines": lineCount,
                        "local_preimages": [r0,r1,r2,r3],"supported_addresses": spatial,
                        "precollision_occurrences": countBefore,"precollision_union": union,
                        "true_band_sparsity": bandS,
                        "address_walsh": (0..<8).map {gpu[rowStart+8*beta+$0]},
                    ]
                }
            }
            try check(!mismatch, "All-coefficient CPU/GPU mismatch")
            try check(energy == ambientSize*ambientSize, "Global Boolean Parseval fails")
            try check(12*globalS >= densityLowerNumerator, "Total weighted band lower bound fails")
            if 12*globalS == densityLowerNumerator { totals.globalMassEqualities += 1 }
            totals.histogramBySparsity[globalS] += 1
            if row == 0 || row == slot.count-1 {
                let choices = [c0,c1,c2,c3]
                let truth = (0..<ambientSize).map { point -> Int32 in
                    let y = point>>3
                    let loc = Int(addressToLine[point&7])
                    return choices[loc>>1].truth[16*(loc&1)+y]
                }
                let direct = walsh(truth)
                try check(direct.indices.allSatisfy { direct[$0] == gpu[rowStart+$0] },
                          "Direct whole-function CPU audit fails")
                totals.fullDirectCases += 1
                if totals.directWitnesses.count < 8 {
                    let le = direct.map { $0.littleEndian }
                    totals.directWitnesses.append(["case_id":id, "coefficients":direct.count,
                                                  "spectrum_sha256_le_int32":le.withUnsafeBytes { sha(Data($0)) }])
                }
            }
        }
        totals.cases += slot.count
        totals.spectrumCoefficients += slot.count*ambientSize
        slot.command = nil
    }

    var next = first
    var submitted = 0
    while next < first+count {
        let slot = slots[submitted % slots.count]
        try collect(slot)
        let span = min(batch,first+count-next)
        let command = queues[submitted % queues.count].makeCommandBuffer()!
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(bankGPU,offset:0,index:0)
        encoder.setBuffer(geometryGPU,offset:0,index:1)
        encoder.setBuffer(slot.buffer,offset:0,index:2)
        var start = UInt32(next)
        encoder.setBytes(&start,length:MemoryLayout<UInt32>.stride,index:3)
        encoder.dispatchThreadgroups(MTLSize(width:span,height:1,depth:1),
                                     threadsPerThreadgroup:MTLSize(width:64,height:1,depth:1))
        encoder.endEncoding()
        slot.command = command
        slot.first = next
        slot.count = span
        command.commit()
        next += span
        submitted += 1
    }
    for slot in slots { try collect(slot) }
    try check(totals.cases == count && totals.spectrumCoefficients == count*ambientSize,"Lost functions")
    let elapsed = ProcessInfo.processInfo.systemUptime-uptimeStart
    let raw: [String:Any] = [
        "schema": "a28-r050-four-line-twisted-bank-metal-v1",
        "device":device.name,"pid":ProcessInfo.processInfo.processIdentifier,
        "address_lines":linePoints,"data_bits":4,"per_line_choices":32,
        "case_interval_half_open":[first,first+count],
        "verified_cases":totals.cases,"cpu_gpu_integer_coefficients_verified":totals.spectrumCoefficients,
        "tested_bands_including_empty":totals.bands,"tested_nonempty_bands":totals.activeBands,
        "private_bands":totals.privateBands,"private_bands_with_two_local_preimages":totals.privateTwoPreimages,
        "colliding_bands":totals.collidingBands, "bands_with_complete_coefficient_cancellation":totals.cancelingBands,
        "completely_canceled_coefficient_occurrences":totals.cancellations,
        "band_line_mass_equalities":totals.bandMassEqualities,"band_actual_address_mass_equalities":totals.actualAddressMassEqualities,
        "global_line_mass_equalities":totals.globalMassEqualities,
        "direct_cpu_full_walsh_cases":totals.fullDirectCases,"direct_cpu_full_walsh_sample_evidence":totals.directWitnesses,
        "walsh_sparsity_histogram":totals.histogramBySparsity,"witnesses":totals.selectedWitnesses,
        "bank_truth_sha256_le_int32":words.map { $0.littleEndian }.withUnsafeBytes { sha(Data($0)) },
        "gpu_command_buffers":submitted,"gpu_command_seconds_sum":totals.gpuSeconds,
        "gpu_command_utc_unix_intervals":totals.gpuWindows,"seconds_excluding_compile":elapsed,
        "all_checks_passed":true,
    ]
    let out = try JSONSerialization.data(withJSONObject:raw,options:[.prettyPrinted,.sortedKeys])
    try (out+Data([10])).write(to:URL(fileURLWithPath:args[1]))
    print("R50 verified \(count) functions and \(count*ambientSize) coefficients on \(device.name)")
}

do { try run() }
catch {
    FileHandle.standardError.write(Data(("R50 verification failed: \(error)\n").utf8))
    exit(1)
}
