import Foundation
import Metal
import CryptoKit

let kernelSource = #"""
#include <metal_stdlib>
using namespace metal;
kernel void fwht_stage(device const int *input [[buffer(0)]],
                       device int *output [[buffer(1)]],
                       constant uint &n [[buffer(2)]],
                       constant uint &stride [[buffer(3)]],
                       uint gid [[thread_position_in_grid]]) {
    uint local = gid % n;
    uint base = gid - local;
    uint mate = base + (local ^ stride);
    int left = input[gid];
    int right = input[mate];
    output[gid] = (local & stride) ? (right - left) : (left + right);
}
"""#

struct CaseRow {
    let n: Int
    let k: Int
    let f: [UInt64]
    let g: [UInt64]
}

func sha256(_ url: URL) -> String {
    let data = try! Data(contentsOf: url)
    return SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
}

func parseWord(_ item: Any) -> UInt64 {
    let s = item as! String
    return UInt64(s.dropFirst(2), radix: 16)!
}

func parseCases(_ url: URL) -> [CaseRow] {
    let root = try! JSONSerialization.jsonObject(with: Data(contentsOf: url)) as! [String: Any]
    let rows = root["cases"] as! [[String: Any]]
    return rows.map { row in
        let n = row["N"] as! Int
        let k = row["k"] as! Int
        let f = (row["f_truth_words_le"] as! [Any]).map(parseWord)
        let g = (row["g_truth_words_le"] as! [Any]).map(parseWord)
        precondition(f.count == (n + 63) / 64 && g.count == (n + 63) / 64)
        return CaseRow(n: n, k: k, f: f, g: g)
    }
}

func bit(_ words: [UInt64], _ index: Int) -> Bool {
    return ((words[index >> 6] >> UInt64(index & 63)) & 1) != 0
}

func cpuFWHT(_ words: [UInt64], _ n: Int) -> [Int64] {
    var a = (0..<n).map { bit(words, $0) ? Int64(-1) : Int64(1) }
    var h = 1
    while h < n {
        var start = 0
        while start < n {
            for i in start..<(start + h) {
                let u = a[i], v = a[i + h]
                a[i] = u + v
                a[i + h] = u - v
            }
            start += 2 * h
        }
        h *= 2
    }
    return a
}

func runBatch(device: MTLDevice, pipeline: MTLComputePipelineState, wordsRows: [[UInt64]], n: Int) throws -> [[Int32]] {
    let phases = wordsRows.count
    var input = [Int32](repeating: 0, count: phases * n)
    for p in 0..<phases {
        for x in 0..<n { input[p * n + x] = bit(wordsRows[p], x) ? -1 : 1 }
    }
    let bytes = input.count * MemoryLayout<Int32>.stride
    let options: MTLResourceOptions = .storageModeShared
    let inBuffer = input.withUnsafeBytes { raw in
        device.makeBuffer(bytes: raw.baseAddress!, length: bytes, options: options)!
    }
    let outBuffer = device.makeBuffer(length: bytes, options: options)!
    var currentIn = inBuffer, currentOut = outBuffer
    let queue = device.makeCommandQueue()!
    let command = queue.makeCommandBuffer()!
    var stride = 1
    while stride < n {
        let encoder = command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(currentIn, offset: 0, index: 0)
        encoder.setBuffer(currentOut, offset: 0, index: 1)
        var n32 = UInt32(n), s32 = UInt32(stride)
        encoder.setBytes(&n32, length: MemoryLayout<UInt32>.stride, index: 2)
        encoder.setBytes(&s32, length: MemoryLayout<UInt32>.stride, index: 3)
        let total = phases * n
        let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: total, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding() // tracked hazards order this stage before the next encoder.
        swap(&currentIn, &currentOut)
        stride <<= 1
    }
    command.commit()
    command.waitUntilCompleted()
    guard command.status == .completed else { throw NSError(domain: "metal", code: Int(command.status.rawValue)) }
    let values = currentIn.contents().bindMemory(to: Int32.self, capacity: phases * n)
    return (0..<phases).map { p in Array(UnsafeBufferPointer(start: values + p * n, count: n)) }
}

let args = CommandLine.arguments
guard args.count == 4 else { fatalError("usage: R20Metal smoke|full cases.json output.json") }
let mode = args[1]
let inputURL = URL(fileURLWithPath: args[2])
let outputURL = URL(fileURLWithPath: args[3])
let allCases = parseCases(inputURL)
let cases = mode == "smoke" ? allCases.filter { $0.n == 16 } : allCases
precondition(mode == "smoke" || mode == "full")
if mode == "smoke" { precondition(cases.count == 24) } else { precondition(cases.count == 120) }
guard let device = MTLCreateSystemDefaultDevice() else { fatalError("no Metal device") }
let library = try! device.makeLibrary(source: kernelSource, options: nil)
let pipeline = try! device.makeComputePipelineState(function: library.makeFunction(name: "fwht_stage")!)
let begun = Date()
var resultGroups: [[String: Any]] = []
var phasesChecked = 0, coefficientsChecked = 0, mismatch = 0, badMagnitude = 0
for n in [16, 64, 256, 1024, 4096] {
    let group = cases.filter { $0.n == n }
    if group.isEmpty { continue }
    for field in ["f", "g"] {
        let masks = group.map { field == "f" ? $0.f : $0.g }
        let cpu = masks.map { cpuFWHT($0, n) }
        let gpu = try! runBatch(device: device, pipeline: pipeline, wordsRows: masks, n: n)
        for p in 0..<gpu.count {
            for i in 0..<n {
                if Int64(gpu[p][i]) != cpu[p][i] { mismatch += 1 }
                if abs(Int64(gpu[p][i])) != Int64(1 << group[p].k) { badMagnitude += 1 }
            }
        }
        phasesChecked += group.count
        coefficientsChecked += group.count * n
    }
    resultGroups.append(["N": n, "cases": group.count, "phases": group.count * 2, "coefficients": group.count * 2 * n])
}
let result: [String: Any] = [
    "schema": "a28-r020-metal-fwht-v1", "mode": mode,
    "device": device.name, "registry_id": String(device.registryID),
    "input_binding": sha256(inputURL), "case_count": cases.count,
    "groups": resultGroups, "phases_checked": phasesChecked,
    "coefficients_checked": coefficientsChecked, "gpu_cpu_mismatch_count": mismatch,
    "bad_abs_walsh_count": badMagnitude,
    "all_checks_passed": mismatch == 0 && badMagnitude == 0,
    "elapsed_seconds": Date().timeIntervalSince(begun),
    "backend": "Metal int32 global staged ping-pong FWHT; Swift int64 CPU elementwise oracle"
]
try! JSONSerialization.data(withJSONObject: result, options: [.prettyPrinted, .sortedKeys]).write(to: outputURL)
if mismatch != 0 || badMagnitude != 0 { exit(10) }
