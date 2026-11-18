#!/usr/bin/env python3
"""Actual signed Fourier factors of Petersen with dense slot, on CPU and Metal."""
import argparse
from array import array
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import struct
import sys
import time

import petersen as pet

SOURCE='''
#include <metal_stdlib>
using namespace metal;

kernel void walsh_half(const device uint *p [[buffer(0)]],
                       device int *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    if (p[0]==0u) {
        const uint a=task&15u, beta=task>>4u;
        int W=1-2*int(popcount(a)&1u);
        for (uint j=0u;j<7u;++j) {
            const uint t=(beta>>(2u*j))&3u;
            int factor=0;
            for (uint x=0u;x<4u;++x) factor+=1-2*int((x==3u)^((popcount(x&t)&1u)!=0u));
            W*=factor;
        }
        out[task]=W; return;
    }
    const uint n=p[1],bits=(1u<<n),mask=bits-1u;
    const uint edge=task>>n,alpha=task&mask;
    const device uint *matrix=p+2u+edge*n;
    int W=0;
    for (uint x=0u;x<bits;++x) {
        uint q=popcount(x&alpha)&1u;
        uint h=x;
        while (h) {
            const uint j=ctz(h);
            h&=h-1u;
            q^=popcount(matrix[j]&h)&1u;
        }
        W+=1-2*int(q);
    }
    out[task]=W;
}
'''


def fwht(values):
    out=list(values); r=1
    while r<len(out):
        for p in range(0,len(out),2*r):
            for q in range(p,p+r):
                a,b=out[q],out[q+r];out[q],out[q+r]=a+b,a-b
        r*=2
    return out


def phase(M,x):
    q=0
    for i in pet.positions(x):q^=pet.parity(M[i]&x&~((1<<(i+1))-1))
    return 1-2*q


def dense14():
    H=fwht(1-2*(sum((x>>(2*i)&3)==3 for i in range(7))%2) for x in range(1<<14))
    if any(abs(v)!=128 for v in H):raise AssertionError('dense bent slot is not flat')
    return [v*(1-2*pet.parity(a)) for v in H for a in range(16)]


def powers(M):
    n=len(M[0]); transforms=[]
    for Q in M:
        W=fwht(phase(Q,x) for x in range(1<<n))
        r=pet.rank(Q)
        if sum(x!=0 for x in W)!=(1<<r) or any(abs(x)!=(1<<(n-r//2)) for x in W if x):
            raise AssertionError('true edge phase does not match exact polar rank')
        transforms.append(W)
    return transforms


def records_from(cpu_file,simul_file):
    x=json.loads(Path(cpu_file).read_text()); u=json.loads(Path(simul_file).read_text())['certificates']
    out=[]
    for key in ('all','mixed_centered'):
        w=x['best'][key]; out.append(dict(title=key,matrices=w['matrices'],t=w['t']))
    for q in u:
        if q['tag'] in ('petersen-macbook.json best nonlocal','two-line compatible dimension 2'):
            out.append(dict(title=q['tag']+' extended',matrices=q['new_matrices'],t=q['updated']['t']))
            if q['tag'].startswith('two-line'):
                out.append(dict(title=q['tag']+' before',matrices=q['old_matrices'],t=q['original']['t']))
    return out


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load '+name)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)


def source_into(a,b,inputs):
    for name,path in (('petersen',a),('metal_base',b)):
        src=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({src!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({inputs!r},True)')


def digest(values):
    return hashlib.sha256(struct.pack('<'+'i'*len(values),*values)).hexdigest()


def run(inputs,use_metal=False):
    d=dense14(); wall=time.perf_counter(); g=None
    if use_metal:
        import metal_base
        metal_base.MSL=SOURCE; g=metal_base.Metal()
    gt0=None;wait0=None
    try:
        if g is not None:
            raw,gt0,wait0=g.compute(array('I',[0]),len(d),1)
            dense_out=array('i');dense_out.frombytes(raw)
            if list(dense_out)!=d:raise AssertionError('dense slot signed CPU/GPU mismatch')
        results=[]
        for item in inputs:
            M=item['matrices']; n=len(M[0]); t=item['t']; transforms=powers(M)
            W=[v for series in transforms for v in series]
            gt=wt=None
            if g is not None:
                params=array('I',[1,n,*[u for Q in M for u in Q]])
                raw,gt,wt=g.compute(params,len(W),1)
                w=array('i');w.frombytes(raw)
                if list(w)!=W: raise AssertionError('graph edge signed CPU/GPU mismatch')
            estimate=pet.estimate(M,t,True); shiftfree=all(series[0]!=0 for series in transforms)
            inside=estimate['inside']; k=estimate['k']
            a_mask=(1<<t)-1
            m=Fraction(sum(abs(v) for e,x in enumerate(transforms) for alpha,v in enumerate(x)
                           if not (inside>>e&1 and not (alpha&a_mask))),1<<n)
            original=Fraction(*estimate['minimum_margin'])
            attained=m-Fraction(3*(k+t),2)
            if attained<original or (shiftfree and attained!=original):
                raise AssertionError('actual signed graph cut inconsistent with polar lower bound')
            # Walsh of a singleton old selector address supplies every one of 16 independent phase signs.
            whole_graph_signs=[(-v if pet.parity(a&e) else v) for e,P in enumerate(transforms)
                               for v in P if v for a in range(16)]
            results.append(dict(title=item['title'],a=n,t=t,partition_blocks=k,inside_mask=inside,
                                actual_edge_phase_transforms=len(W),graph_spectral_points=len(whole_graph_signs),
                                full_true_lift_s=len(whole_graph_signs)+len(d),full_true_R=28+n,
                                old_mass_actual_cut=[m.numerator,m.denominator],
                                actual_cut_margin=[attained.numerator,attained.denominator],
                                cut_bound_margin=[original.numerator,original.denominator],
                                all_phase_fibers_centered=shiftfree,
                                actual_edge_sha256=digest(W),actual_nonzero_graph_sha256=digest(whole_graph_signs),
                                gpu_command_seconds=gt,commit_and_wait_seconds=wt))
    finally:
        if g is not None:g.close()
    print(json.dumps(dict(status='pass',device=None if g is None else g.device_name,host=platform.node(),
                          python=platform.python_version(),dense_signed_coefficients=len(d),dense_sha256=digest(d),
                          dense_gpu_command_seconds=gt0,dense_commit_wait_seconds=wait0,
                          wall_seconds=time.perf_counter()-wall,results=results),indent=2,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',action='store_true');p.add_argument('--metal',action='store_true')
    p.add_argument('petersen_source');p.add_argument('metal_source');p.add_argument('cpu_file');p.add_argument('simul_file')
    a=p.parse_args(); z=records_from(a.cpu_file,a.simul_file)
    if a.bundle:source_into(a.petersen_source,a.metal_source,z)
    else:
        if a.metal:load('metal_base',a.metal_source)
        run(z,a.metal)
