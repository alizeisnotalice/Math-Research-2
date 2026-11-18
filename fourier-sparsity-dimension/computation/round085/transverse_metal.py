#!/usr/bin/env python3
"""Apple Metal exact two-budget audits of actual low-rank phase pencils."""
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

uint insert(thread uint *B,uint x) {
    while (x) {
        uint k=31u-clz(x);
        if (B[k]) x^=B[k];
        else { B[k]=x; return 1u; }
    }
    return 0u;
}

kernel void walsh_half(const device uint *p [[buffer(0)]],
                       device uint *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint case_index=task>>1u,part=task&1u;
    const device uint *record=p+p[1u+case_index];
    const uint R=record[0],a=record[1],t=record[2],s=record[3];
    const device uint *matrices=record+4u;
    const device uint *points=matrices+R*a;
    const uint projmask=(1u<<t)-1u;
    const uint common=1u<<(a/2u);
    uint weighted=0u; uint costs[64];
    for (uint i=0;i<s;++i) {
        const uint coords=points[4u*i+1u],w=points[4u*i+2u],external=points[4u*i+3u];
        uint B[11]={0u},J[11]={0u},r=0u,rank_proj=0u;
        for (uint j=0;j<a;++j) {
            uint row=0u;
            for (uint k=0;k<R;++k) if ((coords>>k)&1u) row^=matrices[k*a+j];
            r+=insert(B,row);
            rank_proj+=insert(J,row&projmask);
        }
        uint g=common<<(r/2u);
        if (!external) g-=g>>rank_proj;
        costs[i]=g;
        weighted+=g*w;
    }
    if (part==0u) { out[task]=weighted; return; }
    uint gb[7]={0u},rank=0u,subtotal=0u,mx=0u;
    for (uint j=0;j<s && rank<R;++j) {
        uint least=0xffffffffu,pos=0u;
        for (uint i=0;i<s;++i) {
            if (costs[i]<least) { least=costs[i]; pos=i; }
        }
        costs[pos]=0xffffffffu;
        if (insert(gb,points[4u*pos])) { subtotal+=least; mx=least; ++rank; }
    }
    out[task]=subtotal+mx;
}
'''


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load '+name)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)


def bundle(paths,witnesses):
    for name,path in zip(('phase_lifts','transverse','metal_base'),paths):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({witnesses!r})')


def choices(tra,seeds,witnesses):
    rng=random.Random(850829)
    for i in range(24576):
        seed=seeds[i%len(seeds)]
        kind=tra.KINDS[(i//len(seeds))%len(tra.KINDS)]
        a=rng.randint(5,min(11,2*seed['R']))
        t=rng.randint(5,a)
        gen=tra.generators(seed,kind,a,rng)
        if gen is None:
            continue
        d=rng.randrange(seed['R']+1)
        U=rng.choice(tra.subspaces(seed['R'])[d])
        yield seed,kind,*gen,a,t,U
    for w in witnesses:
        seed=next(s for s in seeds if s['name']==w['seed'])
        yield seed,'archived-'+w['kind'],w['generators'],w['generator_coords'],w['a'],w['t'],w['U']


def run(witnesses):
    import transverse as tra
    import metal_base as metal
    metal.MSL=SOURCE
    wall=time.perf_counter()
    exp=[]; meta=[]; packing=[]
    seeds=tra.compute_states()
    n_volume=0
    for seed,kind,M,coords,a,t,U in choices(tra,seeds,witnesses):
        exp.append(tra.costs(seed,M,coords,a,t,U)); meta.append((seed['name'],kind))
        x=[seed['R'],a,t,len(seed['E'])]
        x.extend(c for mat in M for c in mat)
        for e,p,w in zip(seed['E'],coords,seed['masses']):
            x.extend([e,p,w,1-((U>>e)&1)])
        packing.append(array('I',x)); n_volume+=len(seed['E'])
    params=array('I',[len(exp)])
    cursor=len(exp)+1
    for x in packing:
        params.append(cursor); cursor+=len(x)
    for x in packing:
        params.extend(x)
    g=metal.Metal()
    try:
        raw,gpu,wait=g.compute(params,2,len(exp))
    finally:
        g.close()
    result=array('I'); result.frombytes(raw)
    if len(result)!=2*len(exp):
        raise AssertionError('incomplete Metal cut output')
    cc=Counter()
    for i,row in enumerate(exp):
        if list(result[2*i:2*i+2])!=[row['weighted_bill'],row['basis_bill']]:
            raise AssertionError(('exact cut budget mismatch',i,meta[i]))
        gr0,gr1=row['greedy_ratio']; e0,e1=row['exact_ratio']
        cc['greedy_insufficient']+=gr0<gr1
        cc['exact_weighted_insufficient']+=e0<e1
        cc['exact_weighted_ties']+=e0==e1
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),device=g.device_name,
                          backend='Metal independent int32 GF2 projected ranks and affine greedy budgets',
                          cases=len(exp),bands=n_volume,counts=cc,seed=850829,sha256=hashlib.sha256(raw).hexdigest(),
                          gpu_command_seconds=gpu,commit_wait_seconds=wait,wall_seconds=time.perf_counter()-wall),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',action='store_true')
    p.add_argument('--witness-file')
    p.add_argument('phase_source')
    p.add_argument('cut_source')
    p.add_argument('metal_source')
    args=p.parse_args()
    paths=(args.phase_source,args.cut_source,args.metal_source)
    w=[] if not args.witness_file else list(json.loads(Path(args.witness_file).read_text())['minima'].values())
    if args.bundle:
        bundle(paths,w)
    else:
        for n,path in zip(('phase_lifts','transverse','metal_base'),paths):
            load(n,path)
        run(w)
