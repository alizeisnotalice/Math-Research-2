#!/usr/bin/env python3
"""Direct integer Walsh and one-cut affine matroid components on Apple Metal."""
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

uint root(thread uint *links, uint i) {
    while (links[i] != i) i = links[i];
    return i;
}

kernel void walsh_half(const device uint *params [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint n = params[0];
    const uint size = 1u << n;
    const uint truth = params[1u + task];
    uint support = 0;
    for (uint u = 0; u < size; ++u) {
        int coefficient = 0;
        for (uint x = 0; x < size; ++x) {
            const uint bit = ((truth >> x) & 1u) ^ (popcount(x & u) & 1u);
            coefficient += 1 - 2 * int(bit);
        }
        support |= uint(coefficient != 0) << u;
    }
    uint signature = support;
    for (uint direction = 1; direction < size; ++direction) {
        uint rows[6] = {0u, 0u, 0u, 0u, 0u, 0u};
        uint expressions[6] = {0u, 0u, 0u, 0u, 0u, 0u};
        uint links[32];
        uint active = support;
        uint components = 0;
        while (active != 0u) {
            const uint u = ctz(active);
            active &= active - 1u;
            if ((popcount(u & direction) & 1u) == 0u) continue;
            ++components;
            links[u] = u;
            uint column = u | size;
            uint expression = 1u << u;
            while (column != 0u) {
                const uint p = 31u - clz(column);
                if (rows[p] == 0u) {
                    rows[p] = column;
                    expressions[p] = expression;
                    break;
                }
                column ^= rows[p];
                expression ^= expressions[p];
            }
            if (column != 0u) continue;
            while (expression != 0u) {
                const uint v = ctz(expression);
                expression &= expression - 1u;
                const uint r = root(links, v);
                if (r == u) continue;
                links[r] = u;
                --components;
            }
        }
        if (components > 2u) {
            out[task] = 0xFFFFFFFFu;
            return;
        }
        signature = signature * 1664525u + components + 31u * direction;
    }
    out[task] = signature;
}
'''


def import_module(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot import '+name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)


def print_bundle(audit_path,backend_path):
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
    wall,cpu = time.perf_counter(),time.process_time()
    gpu = metal_base.Metal()
    records = []
    try:
        for n in (4,5):
            truths = list(audit.datasets(n,n==4))
            params = array('I',[n,*truths])
            payload,gpu_seconds,wait_seconds = gpu.compute(params,1,len(truths))
            signatures = array('I')
            signatures.frombytes(payload)
            if len(signatures) != len(truths):
                raise AssertionError('wrong output length')
            cache = {}
            for truth,got in zip(truths,signatures):
                if got == 0xffffffff:
                    raise AssertionError(('Metal rejects one-cut inequality',n,truth))
                support = audit.walsh_mask(n,truth)
                if support not in cache:
                    expected = support & audit.MASK32
                    for _,direction,cut in audit.direction_cut(n):
                        count = len(audit.small_components(n,support&cut))
                        expected = (expected*1664525+count+31*direction) & audit.MASK32
                    cache[support] = expected
                if got != cache[support]:
                    raise AssertionError(('different CPU and Metal component signatures',n,truth,got,cache[support]))
            records.append(dict(n=n,tables=len(truths),directions=(1<<n)-1,
                                incidences=len(truths)*((1<<n)-1),
                                direction_signature_sha256=hashlib.sha256(payload).hexdigest(),
                                gpu_command_seconds=gpu_seconds,commit_and_wait_seconds=wait_seconds))
    finally:
        gpu.close()
    stop = time.perf_counter()
    gpu_seconds = sum(r['gpu_command_seconds'] for r in records)
    waits = sum(r['commit_and_wait_seconds'] for r in records)
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),
                          backend='Metal direct int32 Walsh and fundamental-circuit connectivity',
                          device=gpu.device_name,records=records,wall_seconds=stop-wall,
                          process_cpu_seconds=time.process_time()-cpu,gpu_command_seconds=gpu_seconds,
                          gpu_command_fraction_of_wait=gpu_seconds/waits,
                          gpu_command_fraction_of_wall=gpu_seconds/(stop-wall)),indent=2,sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('audit_path')
    parser.add_argument('backend_path')
    args = parser.parse_args()
    if args.bundle:
        print_bundle(args.audit_path,args.backend_path)
    else:
        import_module('audit',args.audit_path)
        import_module('metal_base',args.backend_path)
        run()
