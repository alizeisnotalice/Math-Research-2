#!/usr/bin/env python3
"""Cross-check Apple Metal integer Walsh coefficients against a CPU butterfly."""

from array import array
import ctypes as C
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import platform
import random
import time


SEED = 770829
BATCHES = {7: 128, 8: 64, 9: 64, 10: 32, 11: 16, 12: 8}
P = C.c_void_p
U = C.c_ulong
D = C.c_double
MSL = '''
#include <metal_stdlib>
using namespace metal;

kernel void walsh_half(const device uint *truth_bits [[buffer(0)]],
                       device int *out [[buffer(1)]],
                       constant uint &size [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint frequency = task % size;
    const uint base = (task / size) * (size / 32u);
    int total = 0;
    for (uint x = 0u; x < size; ++x) {
        const uint f = (truth_bits[base + x / 32u] >> (x % 32u)) & 1u;
        total += 1 - 2 * int(f ^ (popcount(x & frequency) & 1u));
    }
    out[task] = total / 2;
}
'''


class MetalSize(C.Structure):
    _fields_ = [("width", U), ("height", U), ("depth", U)]


class Metal:
    def __init__(self):
        self.objc = C.CDLL("/usr/lib/libobjc.A.dylib")
        self.foundation = C.CDLL("/System/Library/Frameworks/Foundation.framework/Foundation")
        self.framework = C.CDLL("/System/Library/Frameworks/Metal.framework/Metal")
        self.objc.sel_registerName.argtypes = [C.c_char_p]
        self.objc.sel_registerName.restype = P
        self.objc.objc_getClass.argtypes = [C.c_char_p]
        self.objc.objc_getClass.restype = P
        self.msg_send = C.cast(self.objc.objc_msgSend, P).value
        self.nsstring = self.objc.objc_getClass(b"NSString")
        pool_class = self.objc.objc_getClass(b"NSAutoreleasePool")
        self.pool = self.call(self.call(pool_class, "alloc"), "init")
        self.framework.MTLCreateSystemDefaultDevice.argtypes = []
        self.framework.MTLCreateSystemDefaultDevice.restype = P
        self.device = self.framework.MTLCreateSystemDefaultDevice()
        if not self.device:
            raise RuntimeError("Metal default device is unavailable")
        self.device_name = self.string(self.call(self.device, "name"))
        error = P()
        self.library = self.call(self.device, "newLibraryWithSource:options:error:", P,
                                 (P, P, C.POINTER(P)), self.ns_str(MSL), None, C.byref(error))
        if not self.library:
            raise RuntimeError(self.error(error))
        self.function = self.call(self.library, "newFunctionWithName:", P,
                                  (P,), self.ns_str("walsh_half"))
        self.pipeline = self.call(self.device, "newComputePipelineStateWithFunction:error:", P,
                                  (P, C.POINTER(P)), self.function, C.byref(error))
        if not self.pipeline:
            raise RuntimeError(self.error(error))
        self.queue = self.call(self.device, "newCommandQueue")

    @lru_cache(maxsize=None)
    def method(self, selector, return_type, argument_types):
        function = C.CFUNCTYPE(return_type, P, P, *argument_types)(self.msg_send)
        return function, self.objc.sel_registerName(selector.encode("ascii"))

    def call(self, target, selector, return_type=P, argument_types=(), *args):
        function, registered_selector = self.method(selector, return_type, argument_types)
        return function(target, registered_selector, *args)

    def ns_str(self, value):
        return self.call(self.nsstring, "stringWithUTF8String:", P, (C.c_char_p,), value.encode())

    def string(self, value):
        return self.call(value, "UTF8String", C.c_char_p).decode()

    def error(self, value):
        return self.string(self.call(value, "localizedDescription")) if value else "Unknown Metal error"

    def compute(self, truth_bits, size, count):
        count_frequencies = size * count
        bytes_output = C.sizeof(C.c_int32) * count_frequencies
        bytes_input = len(truth_bits) * truth_bits.itemsize
        data = (C.c_uint32 * len(truth_bits)).from_buffer(truth_bits)
        input_buffer = self.call(self.device, "newBufferWithBytes:length:options:", P,
                                 (P, U, U), C.cast(data, P), bytes_input, 0)
        output_buffer = self.call(self.device, "newBufferWithLength:options:", P, (U, U), bytes_output, 0)
        if not (input_buffer and output_buffer):
            raise RuntimeError("Failed to allocate Metal shared buffers")
        try:
            command = self.call(self.queue, "commandBuffer")
            encoder = self.call(command, "computeCommandEncoder")
            self.call(encoder, "setComputePipelineState:", None, (P,), self.pipeline)
            self.call(encoder, "setBuffer:offset:atIndex:", None, (P, U, U), input_buffer, 0, 0)
            self.call(encoder, "setBuffer:offset:atIndex:", None, (P, U, U), output_buffer, 0, 1)
            size_arg = C.c_uint32(size)
            self.call(encoder, "setBytes:length:atIndex:", None, (P, U, U), C.byref(size_arg), 4, 2)
            group_width = self.call(self.pipeline, "threadExecutionWidth", U)
            self.call(encoder, "dispatchThreads:threadsPerThreadgroup:", None, (MetalSize, MetalSize),
                      MetalSize(count_frequencies, 1, 1), MetalSize(group_width * 4, 1, 1))
            self.call(encoder, "endEncoding", None)
            host_start = time.perf_counter()
            self.call(command, "commit", None)
            self.call(command, "waitUntilCompleted", None)
            host_elapsed = time.perf_counter() - host_start
            status = self.call(command, "status", U)
            if status != 4:
                raise RuntimeError((status, self.error(self.call(command, "error"))))
            gpu_start = self.call(command, "GPUStartTime", D)
            gpu_end = self.call(command, "GPUEndTime", D)
            if not (0 < gpu_start < gpu_end):
                raise RuntimeError((gpu_start, gpu_end, "invalid Metal GPU timestamps"))
            result = C.string_at(self.call(output_buffer, "contents"), bytes_output)
        finally:
            self.call(input_buffer, "release", None)
            self.call(output_buffer, "release", None)
        return result, gpu_end - gpu_start, host_elapsed

    def close(self):
        for owned in (self.queue, self.pipeline, self.function, self.library):
            self.call(owned, "release", None)
        self.call(self.pool, "drain", None)
        self.method.cache_clear()


def cpu_reference(size, truth):
    values = [1 - 2 * ((truth >> x) & 1) for x in range(size)]
    stride = 1
    while stride < size:
        for start in range(0, size, 2 * stride):
            for left in range(start, start + stride):
                right = left + stride
                a, b = values[left], values[right]
                values[left], values[right] = a + b, a - b
        stride *= 2
    return [value // 2 for value in values]


def verify_residue_polynomials(n, table_index, spectrum):
    size = 1 << n
    axes = [1 << i for i in range(n)]
    signs_mod4 = [int(value % 4 == 3) for value in spectrum]
    global_sign_bit = signs_mod4[0]
    affine_mask = sum(axis for axis in axes if signs_mod4[axis] != global_sign_bit)
    magnitude_bits = [int(abs(value) % 8 in (3, 5)) for value in spectrum]
    constant = magnitude_bits[0]
    linear_mask = sum(axis for axis in axes if magnitude_bits[axis] != constant)
    quadratic_masks = [i | j for i, j in combinations(axes, 2)
                       if magnitude_bits[i | j] ^ magnitude_bits[i] ^ magnitude_bits[j] ^ constant]
    for u in range(size):
        if spectrum[u] % 2 != 1:
            raise AssertionError((n, table_index, u, "unexpected even Walsh/2"))
        expected_sign_bit = global_sign_bit ^ (bin(u & affine_mask).count("1") % 2)
        if signs_mod4[u] != expected_sign_bit:
            raise AssertionError((n, table_index, u, "non-affine sign"))
        expected_q = constant ^ (bin(u & linear_mask).count("1") % 2)
        for pair in quadratic_masks:
            expected_q ^= int(u & pair == pair)
        if magnitude_bits[u] != expected_q:
            raise AssertionError((n, table_index, u, "non-quadratic mod-8 magnitude bit"))


def check_dimension(metal, n, count, rng):
    size = 1 << n
    truths = []
    truth_bits = array("I")
    for index in range(count):
        tail = rng.getrandbits(size - 1)
        truth = tail | ((1 ^ (bin(tail).count("1") % 2)) << (size - 1))
        truths.append(truth)
        truth_bits.extend((truth >> offset) & 0xffffffff for offset in range(0, size, 32))
    if truth_bits.itemsize != 4:
        raise RuntimeError("Unsigned packed truth data must use 32-bit words")
    result, gpu_seconds, host_seconds = metal.compute(truth_bits, size, count)
    actual = array("i")
    actual.frombytes(result)
    expected = [half for truth in truths for half in cpu_reference(size, truth)]
    if len(actual) != len(expected):
        raise AssertionError((n, len(expected), len(actual), "coefficient count"))
    for index, (got, want) in enumerate(zip(actual, expected)):
        if got != want:
            raise AssertionError((n, index // size, index % size, want, got))
    for table_index in range(count):
        offset = table_index * size
        verify_residue_polynomials(n, table_index, actual[offset:offset + size])
    return {
        "dimension": n,
        "sample_tables": count,
        "checked_coefficients": len(expected),
        "checked_affine_and_quadratic_residues": len(expected),
        "gpu_command_seconds": gpu_seconds,
        "commit_and_wait_seconds": host_seconds,
        "command_fraction_of_commit_and_wait": gpu_seconds / host_seconds,
        "int32_walsh_half_sha256": hashlib.sha256(result).hexdigest(),
    }


def main():
    rng = random.Random(SEED)
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    metal = Metal()
    try:
        results = [check_dimension(metal, n, count, rng) for n, count in BATCHES.items()]
    finally:
        metal.close()
    wall_seconds = time.perf_counter() - wall_start
    cpu_seconds = time.process_time() - cpu_start
    gpu_seconds = sum(record["gpu_command_seconds"] for record in results)
    wait_seconds = sum(record["commit_and_wait_seconds"] for record in results)
    print(json.dumps({
        "status": "pass",
        "seed": SEED,
        "arithmetic": "Metal signed int32 versus Python exact integers",
        "backend": "Metal compute, System Default Device",
        "device": metal.device_name,
        "python": platform.python_version(),
        "host": platform.node(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "wall_seconds": wall_seconds,
        "process_cpu_seconds": cpu_seconds,
        "total_gpu_command_seconds": gpu_seconds,
        "total_commit_and_wait_seconds": wait_seconds,
        "gpu_command_fraction_of_wall": gpu_seconds / wall_seconds,
        "gpu_command_fraction_of_commit_and_wait": gpu_seconds / wait_seconds,
        "records": results,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
