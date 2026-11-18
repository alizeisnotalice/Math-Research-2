#!/usr/bin/env python3
"""Metal audit of every connected Petersen partition for selected projected pencils."""
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

SOURCE='''
#include <metal_stdlib>
using namespace metal;

uint add(thread uint *p,uint v) {
    while (v) {
        const uint i=31u-clz(v);
        if (p[i]) v^=p[i];
        else {p[i]=v;return 1u;}
    }
    return 0u;
}

kernel void walsh_half(const device uint *p [[buffer(0)]],
                       device int *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint np=p[0];
    const device uint *parts=p+1u;
    const device uint *index=parts+2u*np;
    const device uint *record=p+index[task];
    const uint a=record[0],t=record[1],proj=(1u<<t)-1u;
    uint r[15],p_rank[15],mr=0u;
    for (uint e=0;e<15;++e) {
        uint B[20]={0u},C[20]={0u};
        uint rb=0u,rc=0u;
        for (uint d=0;d<a;++d) {
            uint row=record[2u+e*a+d];
            rb+=add(B,row);
            rc+=add(C,row&proj);
        }
        r[e]=rb;p_rank[e]=rc;mr=max(mr,rb);
    }
    const uint K=1u<<(mr/2u);
    uint d[15],sum=0u;
    for (uint e=0;e<15;++e) {
        const uint c=K<<(r[e]/2u);
        sum+=c;
        d[e]=c>>p_rank[e];
    }
    uint bg=0u;
    for (uint i=0;i<np;++i) {
        uint score=3u*K*parts[2u*i],inside=parts[2u*i+1u];
        while (inside) {
            uint j=ctz(inside);
            inside&=inside-1u;
            score+=2u*d[j];
        }
        bg=max(bg,score);
    }
    out[task]=int(2u*sum)-int(bg)-int(3u*K*t);
}
'''


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError('cannot import '+name)
    m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m)


def bundle(paths,extras):
    for name,path in zip(('petersen','metal_base'),paths):
        src=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({src!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({extras!r})')


def run(extras):
    import petersen as pet
    import metal_base as metal
    metal.MSL=SOURCE
    start=time.perf_counter(); rng=random.Random(860829)
    graph=pet.partitions()
    cases=[]
    for i in range(1536):
        kind=pet.KINDS[i%len(pet.KINDS)]
        M,center=pet.choose(kind,rng); rec=pet.normalize(M,8,rng)
        if rec is not None:
            M,_=rec
            t=rng.randrange(8,len(M[0])+1)
            cases.append((M,t,kind))
    for item in extras:
        M=item['matrices']; t=item['t']
        cases.append((M,t,item['kind']))
    params=array('I',[len(graph)])
    for k,m,_ in graph: params.extend((k,m))
    offset=len(params)+len(cases)
    for M,t,_ in cases:
        params.append(offset);offset+=2+15*len(M[0])
    for M,t,_ in cases:
        params.extend((len(M[0]),t))
        params.extend(x for Q in M for x in Q)
    gpu=metal.Metal()
    try:
        got,gt,wt=gpu.compute(params,len(cases),1)
    finally:
        gpu.close()
    out=array('i');out.frombytes(got)
    if len(out)!=len(cases):raise AssertionError('wrong graph cut output size')
    most=None; neg=eq=0
    for (M,t,k),x in zip(cases,out):
        answer=pet.estimate(M,t,True)
        y,scale=answer['minimum_margin']
        if x!=y:raise AssertionError((k,t,x,y,'GPU vs CPU whole graph partition mismatch'))
        neg+=x<0;eq+=x==0
    print(json.dumps(dict(status='pass',gpu_name=gpu.device_name,host=platform.node(),python=platform.python_version(),
                          cases=len(cases),partitions_per_pencil=len(graph),exact_partition_checks=len(graph)*len(cases),
                          negative=neg,equalities=eq,seed=860829,sha256=hashlib.sha256(got).hexdigest(),
                          gpu_command_seconds=gt,commit_and_wait_seconds=wt,wall_seconds=time.perf_counter()-start),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',action='store_true')
    p.add_argument('--pencil-file')
    p.add_argument('petersen_source')
    p.add_argument('metal_source')
    a=p.parse_args()
    extras=[] if not a.pencil_file else list(json.loads(Path(a.pencil_file).read_text())['best'].values())
    paths=(a.petersen_source,a.metal_source)
    if a.bundle:bundle(paths,extras)
    else:
        load('petersen',a.petersen_source);load('metal_base',a.metal_source);run(extras)
