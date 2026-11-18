#!/usr/bin/env python3
"""True Walsh and connectivity sweep on both Apple Metal devices."""
import argparse
from array import array
from collections import Counter
import hashlib
import importlib.util
from itertools import combinations
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

uint component_count(uint extent, uint size) {
    uint rows[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    uint expressions[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    uint links[32];
    uint count = 0;
    while (extent != 0u) {
        const uint u = ctz(extent);
        extent &= extent - 1;
        ++count;
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
            expression &= expression - 1;
            const uint r = root(links,v);
            if (r == u) continue;
            links[r] = u;
            --count;
        }
    }
    return count;
}

uint direct_walsh_support(uint truth, uint size) {
    uint support = 0;
    for (uint u = 0; u < size; ++u) {
        int coefficient = 0;
        for (uint x = 0; x < size; ++x) {
            uint bit = ((truth >> x) & 1u) ^ (popcount(x & u) & 1u);
            coefficient += 1 - 2*int(bit);
        }
        support |= uint(coefficient != 0) << u;
    }
    return support;
}

kernel void walsh_half(const device uint *params [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint mode = params[0];
    const uint size = 1u << params[1];
    if (mode == 0u) {
        const uint kernels = params[2];
        const uint truth = params[3u + 2u*kernels + task];
        const uint support = direct_walsh_support(truth, size);
        const uint base = support & ((1u << (size/4u)) - 1u);
        uint signature = base;
        for (uint j=0; j<kernels; ++j) {
            const uint H = params[3u+2u*j];
            const uint codim = params[4u+2u*j];
            const uint count = component_count(support & ~H,size);
            if ((base & (base-1u)) && count>1u) {
                out[task]=0xFFFFFFFFu;
                return;
            }
            signature = signature*1664525u + count + 31u*codim;
        }
        out[task]=signature;
    } else {
        const uint truth = params[2];
        const uint H = params[3u+2u*task];
        const uint D = params[4u+2u*task];
        const uint support = direct_walsh_support(truth,size);
        out[task] = component_count(support & ~(H | D),size);
    }
}
'''


def import_module(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot import '+name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)


def print_bundle(paths):
    for name,path in zip(('check_components','bent_tensor','metal_base'),paths):
        source = Path(path).read_text()
        print('import types, sys')
        print(f'module = types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}] = module')
        print(f'exec(compile({source!r}, {name!r}, "exec"), module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print('run()')


def audit_primary(gpu,comp,tensor,n):
    truths=list(tensor.datasets(n,True))
    kernels=comp.spaces(n+2)
    params=array('I',[0,n+2,len(kernels)])
    for codim,H in kernels:
        params.extend([H,codim])
    params.extend(tensor.truth_product(n,t) for t in truths)
    payload,gpu_seconds,wait_seconds=gpu.compute(params,1,len(truths))
    result=array('I')
    result.frombytes(payload)
    if len(result)!=len(truths):
        raise AssertionError('incomplete Metal output')
    direct={}
    for t,digest in zip(truths,result):
        S,T=tensor.tensor_support(n,t)
        if S not in direct:
            direct[S]=tensor.all_kernel_connectivity(n,S,T)[2]
        if digest!=direct[S]:
            raise AssertionError((n,t,digest,direct[S],'Metal and CPU disagree'))
    return dict(operation='all subspaces',base_dimension=n,ambient_dimension=n+2,tables=len(truths),
                kernel_cases=len(kernels)*len(truths),signature_sha256=hashlib.sha256(payload).hexdigest(),
                gpu_command_seconds=gpu_seconds,commit_and_wait_seconds=wait_seconds)


def audit_deletions(gpu,comp,tensor):
    n=4
    truth=tensor.truth_product(2,1)
    S=comp.walsh_mask(n,truth)
    cases=[]
    expected=[]
    params=array('I',[1,n,truth])
    for codim,H in comp.spaces(n):
        for d in (0,1,2,3):
            for deleted in combinations(comp.points(S&~H),d):
                D=sum(1<<u for u in deleted)
                cases.append((n-codim,d))
                params.extend([H,D])
                expected.append(len(comp.small_components(n,S&~(H|D))))
    raw,gpu_seconds,wait_seconds=gpu.compute(params,1,len(cases))
    result=array('I')
    result.frombytes(raw)
    if list(result)!=expected:
        raise AssertionError('Metal and CPU deletion sweeps differ')
    failures=Counter((H,d) for (H,d),k in zip(cases,result) if k>1)
    first=min(failures,default=None)
    return dict(operation='all subspaces and up to three deletions',tables=1,ambient_dimension=n,
                kernel_cases=len(cases),first_failure=first,
                failures=[dict(dimH=H,D=d,count=c) for (H,d),c in sorted(failures.items())],
                signature_sha256=hashlib.sha256(raw).hexdigest(),gpu_command_seconds=gpu_seconds,
                commit_and_wait_seconds=wait_seconds)


def run():
    import bent_tensor as tensor
    import check_components as comp
    import metal_base
    metal_base.MSL=SOURCE
    wall,cpu=time.perf_counter(),time.process_time()
    gpu=metal_base.Metal()
    try:
        records=[audit_primary(gpu,comp,tensor,n) for n in (1,2,3)]
        records.append(audit_deletions(gpu,comp,tensor))
    finally:
        gpu.close()
    gpu_time=sum(r['gpu_command_seconds'] for r in records)
    wait=sum(r['commit_and_wait_seconds'] for r in records)
    print(json.dumps(dict(status='pass',device=gpu.device_name,backend='Metal direct int32 Walsh and GF(2) circuits',
                          host=platform.node(),python=platform.python_version(),records=records,
                          wall_seconds=time.perf_counter()-wall,process_cpu_seconds=time.process_time()-cpu,
                          total_gpu_command_seconds=gpu_time,gpu_fraction_of_commit_and_wait=gpu_time/wait),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('components_path')
    parser.add_argument('tensor_path')
    parser.add_argument('backend_path')
    args=parser.parse_args()
    paths=[args.components_path,args.tensor_path,args.backend_path]
    if args.bundle:
        print_bundle(paths)
    else:
        for name,path in zip(('check_components','bent_tensor','metal_base'),paths):
            import_module(name,path)
        run()
