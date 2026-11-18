#!/usr/bin/env python3
"""Integer Metal audit computing Walsh support and affine rank independently."""
import argparse
from array import array
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time


SOURCE = '''
#include <metal_stdlib>
using namespace metal;

kernel void walsh_half(const device uint *param [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint n = param[0];
    const uint total_flats = param[1];
    const uint size = 1u << n;
    const uint truth = param[2u + 2u * total_flats + task];
    uint support = 0;
    for (uint u = 0; u < size; ++u) {
        int coefficient = 0;
        for (uint x = 0; x < size; ++x) {
            const uint fx = (truth >> x) & 1u;
            coefficient += 1 - 2 * int(fx ^ (popcount(x & u) & 1u));
        }
        support |= uint(coefficient != 0) << u;
    }
    uint signature = support;
    for (uint index = 0; index < total_flats; ++index) {
        uint sample = support & param[2u + 2u * index];
        const uint codim = param[3u + 2u * index];
        const uint q = popcount(sample);
        uint basis[6] = {0u, 0u, 0u, 0u, 0u, 0u};
        uint rank = 0;
        while (sample != 0u) {
            uint column = size | ctz(sample);
            sample &= sample - 1;
            while (column != 0u) {
                const uint pivot = 31u - clz(column);
                if (basis[pivot] == 0u) {
                    basis[pivot] = column;
                    ++rank;
                    break;
                }
                column ^= basis[pivot];
            }
        }
        if (q > 2u * (q - rank) + (1u << (codim + 1u))) {
            out[task] = 0xFFFFFFFFu;
            return;
        }
        signature = signature * 1664525u + q + 257u * rank + 65537u * codim;
    }
    out[task] = signature;
}
'''


def import_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError('module unavailable: ' + name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def print_bundle(audit_path, backend_path):
    for name,path in [('audit',audit_path),('metal_base',backend_path)]:
        source = Path(path).read_text()
        print('import types, sys')
        print(f'module = types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}] = module')
        print(f'exec(compile({source!r}, {name!r}, "exec"), module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__ ==',1)[0])
    print('run()')


def run():
    import audit
    import metal_base
    metal_base.MSL = SOURCE
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    metal = metal_base.Metal()
    records = []
    try:
        for n in (4,5):
            exhaustive = n == 4
            flats = audit.affine_flats(n)
            truths = list(audit.dataset(n,exhaustive))
            params = array('I',[n,len(flats)])
            for t,flat in flats:
                params.extend([flat,t])
            params.extend(truths)
            raw,gpu_seconds,wait_seconds = metal.compute(params,1,len(truths))
            got = array('I')
            got.frombytes(raw)
            cache = {}
            for truth,signature in zip(truths,got):
                if signature == 0xFFFFFFFF:
                    raise AssertionError(('GPU reported a violated inequality',n,truth))
                support = audit.signing_from_table(n,truth)
                if support not in cache:
                    cache[support] = audit.evaluate_support(n,support,flats)[0]
                if signature != cache[support]:
                    raise AssertionError(('CPU and Metal signatures differ',n,truth,signature,cache[support]))
            if len(got) != len(truths):
                raise AssertionError('Metal returned wrong dimensions')
            records.append(dict(n=n,tables=len(truths),flats=len(flats),incidences=len(truths)*len(flats),
                                gpu_command_seconds=gpu_seconds,commit_and_wait_seconds=wait_seconds,
                                table_signature_sha256=hashlib.sha256(raw).hexdigest()))
    finally:
        metal.close()
    wall_seconds = time.perf_counter()-wall_start
    gpu_seconds = sum(r['gpu_command_seconds'] for r in records)
    waits = sum(r['commit_and_wait_seconds'] for r in records)
    print(json.dumps(dict(status='pass',records=records,backend='Metal int32 direct Walsh and GF2 rank',
                          host=platform.node(),python=platform.python_version(),device=metal.device_name,
                          wall_seconds=wall_seconds,process_cpu_seconds=time.process_time()-cpu_start,
                          gpu_command_seconds=gpu_seconds,gpu_command_fraction_of_wait=gpu_seconds/waits,
                          gpu_command_fraction_of_wall=gpu_seconds/wall_seconds),indent=2,sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('audit_path')
    parser.add_argument('metal_backend_path')
    args = parser.parse_args()
    if args.bundle:
        print_bundle(args.audit_path,args.metal_backend_path)
    else:
        import_module('audit',args.audit_path)
        import_module('metal_base',args.metal_backend_path)
        run()
