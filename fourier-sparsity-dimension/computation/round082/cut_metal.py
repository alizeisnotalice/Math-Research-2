#!/usr/bin/env python3
"""Independent integer Metal Walsh and exact-cut audit on two Apple GPUs."""
import argparse
from array import array
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import sys
import time

SOURCE = '''
#include <metal_stdlib>
using namespace metal;

uint insert_column(thread uint *basis, uint v) {
    while (v) {
        const uint pivot = 31u - clz(v);
        const uint row = basis[pivot];
        if (row == 0u) {
            basis[pivot] = v;
            return 1u;
        }
        v ^= row;
    }
    return 0u;
}

kernel void walsh_half(const device uint *parameters [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint n = parameters[0];
    const uint size = 1u << n;
    const uint offset = 1u + task * (size/32u);
    int W[256];
    uint full_basis[9] = {0u};
    uint R = 0u, s = 0u;
    for (uint u=0; u<size; ++u) {
        int total = 0;
        for (uint x=0; x<size; ++x) {
            const uint bit = (parameters[offset + x/32u] >> (x % 32u)) & 1u;
            total += 1 - (2 * int(bit ^ (popcount(u&x)&1u)));
        }
        W[u] = total;
        if (total != 0) {
            R += insert_column(full_basis,size | u);
            ++s;
        }
    }
    uint signature = 0u;
    for (uint query=1; query<size; ++query) {
        const uint leading = 1u << (31u-clz(query));
        uint basis0[9] = {0u}, basis1[9] = {0u}, union_basis[9] = {0u};
        uint R0=0u,R1=0u,Runion=0u,s0=0u,s1=0u,k=0u,e=0u;
        for (uint q=0; q<size; ++q) {
            if (q & leading) continue;
            const int a = W[q];
            const int b = W[q^query];
            const int left = a+b, right=a-b;
            if (left != 0) {
                ++s0;
                R0 += insert_column(basis0,size | q);
            }
            if (right != 0) {
                ++s1;
                R1 += insert_column(basis1,size | q);
            }
            if (left != 0 || right != 0) Runion += insert_column(union_basis,size | q);
            if (left != 0 && right != 0) {
                ++k;
                e += uint(abs(left)==abs(right));
            }
        }
        const int U=int(R0+R1+2u), L=int(s0+s1), j=int(R0+R1-Runion);
        const int parent_sq=int((R+1u)*(R+1u));
        const uint proof=uint(2*L*(2*U*j-j*j) >= int(2u*k+e)*U*U);
        const uint proof_sharp=uint(2*L*(U*U-parent_sq) >= int(2u*k+e)*U*U);
        const uint nondecrease0=uint((R0+1u)*(R0+1u)*s >= uint(parent_sq)*s0);
        const uint nondecrease1=uint((R1+1u)*(R1+1u)*s >= uint(parent_sq)*s1);
        if (s!=2u*(s0+s1)-2u*k-e || (R != Runion && R!=Runion+1u)
            || ((k<=uint(j)||proof||proof_sharp)&&(!nondecrease0&&!nondecrease1))) {
            out[task]=0xFFFFFFFFu;
            return;
        }
        const uint facts[11]={s0,s1,R0,R1,k,uint(j),e,proof,proof_sharp,nondecrease0,nondecrease1};
        for (uint i=0; i<11u; ++i) signature=signature*1664525u+facts[i];
    }
    out[task]=signature;
}
'''

SEED=820829
SAMPLE_PER_DIM=512


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot import '+name)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)


def bundle(paths, only_witnesses):
    for name,path in zip(('scan_cut','metal_base'),paths):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({only_witnesses!r})')


def run(only_witnesses=False):
    import scan_cut as scan
    import metal_base as metal
    metal.MSL=SOURCE
    wall,cpu=time.perf_counter(),time.process_time()
    gpu=metal.Metal()
    results=[]
    try:
        for n in ((6,) if only_witnesses else range(5,9)):
            rng=random.Random(SEED+n)
            random_tables=[scan.random_truth(n,scan.FAMILIES[i%len(scan.FAMILIES)],rng) for i in range(SAMPLE_PER_DIM)]
            fixtures=[truth for name,n0,truth in scan.structured_fixtures() if n0==n]
            if only_witnesses:
                fixtures=[1599346001363546751,11679430202935411900,18223943850446558208]
                random_tables=[]
            truths=fixtures+random_tables
            params=array('I',[n])
            for truth in truths:
                params.extend((truth>>shift)&0xffffffff for shift in range(0,1<<n,32))
            raw,gpu_seconds,wait_seconds=gpu.compute(params,1,len(truths))
            records=array('I')
            records.frombytes(raw)
            if len(records)!=len(truths):
                raise AssertionError('GPU returned wrong length')
            for truth,digest in zip(truths,records):
                if digest==0xffffffff or digest != scan.analyze(truth,n)['signature']:
                    raise AssertionError((n,truth,digest,'CPU/GPU difference or invalid inequality'))
            results.append(dict(n=n,tables=len(truths),fixture_tables=len(fixtures),sampled_tables=len(random_tables),
                                incidences=len(truths)*((1<<n)-1),signature_sha256=hashlib.sha256(raw).hexdigest(),
                                gpu_command_seconds=gpu_seconds,commit_and_wait_seconds=wait_seconds))
    finally:
        gpu.close()
    gpu_seconds=sum(r['gpu_command_seconds'] for r in results)
    wait_seconds=sum(r['commit_and_wait_seconds'] for r in results)
    print(json.dumps(dict(status='pass',device=gpu.device_name,host=platform.node(),python=platform.python_version(),
                          seed=SEED,all_data=results,backend='Metal direct Walsh and int32 cut tests',
                          wall_seconds=time.perf_counter()-wall,process_cpu_seconds=time.process_time()-cpu,
                          gpu_command_seconds=gpu_seconds,gpu_command_fraction_of_wait=gpu_seconds/wait_seconds),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('--only-witnesses',action='store_true')
    parser.add_argument('scan_path')
    parser.add_argument('backend_path')
    args=parser.parse_args()
    if args.bundle:
        bundle((args.scan_path,args.backend_path),args.only_witnesses)
    else:
        load('scan_cut',args.scan_path)
        load('metal_base',args.backend_path)
        run(args.only_witnesses)
