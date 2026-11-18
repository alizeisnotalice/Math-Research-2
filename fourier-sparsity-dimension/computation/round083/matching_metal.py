#!/usr/bin/env python3
"""Independent Apple Metal signed sparse restrictions and rank census."""
import argparse
from array import array
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import sys
import time

SOURCE='''
#include <metal_stdlib>
using namespace metal;

int value_at(const device uint *frequencies, uint m, uint u) {
    uint lo=0,hi=m;
    while (lo<hi) {
        uint k=(lo+hi)/2;
        if (frequencies[2*k]<u) lo=k+1; else hi=k;
    }
    if (lo<m && frequencies[2*lo]==u) return as_type<int>(frequencies[2*lo+1]);
    return 0;
}

uint add_rank(thread uint *basis, uint v) {
    while (v) {
        uint i=31u-clz(v);
        if (basis[i]) { v^=basis[i]; }
        else { basis[i]=v; return 1u; }
    }
    return 0u;
}

kernel void walsh_half(const device uint *p [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint size=1u<<p[0], count=p[1];
    const device uint *freq=p+2;
    const device uint *queries=freq+2*count;
    const uint group=task/2u, pair=task%2u;
    const uint a=queries[2u*group], b=queries[2u*group+1u];
    const uint width=b?4u:2u;
    if (width==2u && pair==1u) {out[task]=0u;return;}
    const uint offsets[4]={0u,a,b,a^b};
    uint bas0[19]={0u}, bas1[19]={0u};
    uint s0=0u,s1=0u,r0=0u,r1=0u;
    for (uint i=0;i<count;++i) {
        const uint u=freq[2u*i];
        uint q=u;
        bool first=true;
        for (uint j=1u;j<width;++j) {
            uint x=u^offsets[j];
            if (x<q) q=x;
            if (x<u && value_at(freq,count,x)!=0) { first=false; break; }
        }
        if (!first) continue;
        int sum0=0,sum1=0;
        for (uint j=0;j<width;++j) {
            int value=value_at(freq,count,q^offsets[j]);
            if (pair&&(j&2u)) value=-value;
            sum0+=value;
            sum1+=(j&1u)?(-value):value;
        }
        if (sum0!=0) { ++s0; r0+=add_rank(bas0,size|q); }
        if (sum1!=0) { ++s1; r1+=add_rank(bas1,size|q); }
    }
    out[task]=s0 | (r0<<9u) | (s1<<14u) | (r1<<23u);
}
'''


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot import '+name)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)


def bundle(paths,names,depth):
    for name,path in zip(('matching_queries','metal_base'),paths):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({names!r},{depth!r})')


def run(names,depth):
    import matching_queries as match
    import metal_base as backend
    backend.MSL=SOURCE
    start=time.perf_counter()
    gpu=backend.Metal()
    records=[]
    try:
        for name in names:
            case=match.build_fixture(name)
            match.verify_fixture(case)
            n=case['n']
            ds=set(match.directions(case))
            queries=match.cuts_for(case) if depth==2 else {1:[(a,) for a in sorted(ds)]}
            histograms=[]
            for h,rows_list in queries.items():
                params=array('I',[case['n'],len(case['coeff'])])
                for point,value in sorted(case['coeff'].items()):
                    params.extend([point,value&0xffffffff])
                for rows in rows_list:
                    params.extend([rows[0],rows[1] if h==2 else 0])
                payload,gpu_seconds,commit_wait_seconds=gpu.compute(params,2,len(rows_list))
                output=array('I'); output.frombytes(payload)
                if len(output)!=2*len(rows_list):
                    raise AssertionError('incomplete GPU output')
                G=Counter()
                indices=set(random.Random(830829).sample(range(len(rows_list)),min(256,len(rows_list))))
                for i,rows in enumerate(rows_list):
                    packed=[output[2*i+j] for j in range(h)]
                    pairs=[p for k in packed for p in [(k&511,(k>>9)&31),((k>>14)&511,(k>>23)&31)]]
                    if i in indices and pairs!=match.child_pairs(case,match.restrictions(case,rows)):
                        raise AssertionError((name,rows,pairs,'signed GPU/CPU mismatch'))
                    if any(not s or not R for s,R in pairs):
                        raise AssertionError('empty child')
                    m=h+(h==2 and rows[0]^rows[1] in ds)
                    a,b=match.queries_type(rows,case['address_dim'],n)
                    G[a,b,m,tuple(sorted(pairs))]+=1
                histograms.append(dict(depth=h,n_kernels=len(rows_list),direct_cpu_cross_checks=len(indices),
                                       sha256=hashlib.sha256(payload).hexdigest(),gpu_command_seconds=gpu_seconds,
                                       commit_and_wait_seconds=commit_wait_seconds,
                                       exact_classes=[dict(address_intersection=a,data_intersection=b,active_differences=m,
                                                           child_pairs=p,n_kernels=c) for (a,b,m,p),c in sorted(G.items())]))
            records.append(dict(name=name,n=n,levels=histograms))
    finally:
        gpu.close()
    gpu_seconds=sum(level['gpu_command_seconds'] for r in records for level in r['levels'])
    wait=sum(level['commit_and_wait_seconds'] for r in records for level in r['levels'])
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),device=gpu.device_name,
                          backend='Metal independent int32 sparse signed coset restrictions',cases=records,
                          wall_seconds=time.perf_counter()-start,gpu_command_seconds=gpu_seconds,
                          gpu_command_fraction_of_wait=gpu_seconds/wait),indent=2,sort_keys=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('--depth',type=int,choices=(1,2),default=2)
    parser.add_argument('case_source')
    parser.add_argument('metal_source')
    parser.add_argument('cases',nargs='+')
    args=parser.parse_args()
    if args.bundle:
        bundle([args.case_source,args.metal_source],args.cases,args.depth)
    else:
        load('matching_queries',args.case_source)
        load('metal_base',args.metal_source)
        run(args.cases,args.depth)
