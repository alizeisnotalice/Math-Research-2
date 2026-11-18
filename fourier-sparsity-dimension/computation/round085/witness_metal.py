#!/usr/bin/env python3
"""Direct integer Metal Walsh of the two reconstructed Boolean witness lifts."""
import argparse
from array import array
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time

SOURCE='''
#include <metal_stdlib>
using namespace metal;

kernel void walsh_half(const device uint *p [[buffer(0)]],
                       device int *out [[buffer(1)]],
                       constant uint &unused [[buffer(2)]],
                       uint task [[thread_position_in_grid]]) {
    const uint n=p[0],a=p[1],w=1u<<n;
    const uint beta=task&(w-1u),alpha=task>>n;
    const device uint *truth=p+2u;
    const device uint *shifts=truth+w;
    int W=0;
    for (uint z=0;z<(1u<<a);++z) {
        const uint shift=shifts[z]&(w-1u);
        const uint phase=(popcount(z&alpha)&1u)^(shifts[z]>>n);
        for (uint x=0;x<w;++x) {
            uint f=truth[x^shift]^phase^(popcount(beta&x)&1u);
            W+=1-2*int(f);
        }
    }
    out[task]=W;
}
'''


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load '+name)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)


def bundle(paths,witnesses,additional):
    for name,path in zip(('phase_lifts','transverse','verify_witnesses','metal_base'),paths):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({witnesses!r},{additional!r})')


def run(witnesses,additional):
    import phase_lifts as phase
    import transverse as tra
    import verify_witnesses as check
    import metal_base
    metal_base.MSL=SOURCE
    start=time.perf_counter(); results=[]
    gpu=metal_base.Metal()
    try:
        seeds={s['name']:s for s in tra.compute_states()}
        for name,item in additional.items():
            n=item['n']; size=1<<n
            W=phase.fwht(1-2*((item['truth_mask']>>x)&1) for x in range(size))
            c={u:w for u,w in enumerate(W) if w}
            state=dict(n=n,denom=size,coeff=c)
            seeds[name]=dict(name=name,n=n,R=n+1,coefficients=c,parent=state,E=[u|size for u in c])
        for label,w in sorted(witnesses.items()):
            seed=seeds[w['seed']]; a=w['a']; n=seed['n']; count=1<<(a+n)
            units=check.all_homogeneous_generator_coords(seed,w['generator_coords'])
            generators=w['generators']
            matrices=[tra.sum_matrices([generators[j] for j in tra.points(c)],a) for c in units]
            params=array('I',[n,a])
            params.extend(int(v==-1) for v in phase.truth_from_sparse(seed['parent']))
            params.extend(sum(check.quadratic(M,z)<<i for i,M in enumerate(matrices)) for z in range(1<<a))
            raw,gpu_t,wait=gpu.compute(params,count,1)
            W=array('i');W.frombytes(raw)
            _,expected,_=check.reconstruct(seed,w)
            if list(W)!=expected:
                raise AssertionError(('direct CPU and Metal Boolean Walsh differ',label))
            results.append(dict(label=label,parent_homogeneous_rank=a+n+1,n_signed_coefficients=count,
                                support_size=sum(v!=0 for v in W),sum_abs_unnormalized=sum(abs(x) for x in W),
                                sha256=hashlib.sha256(raw).hexdigest(),gpu_command_seconds=gpu_t,
                                commit_and_wait_seconds=wait))
    finally:
        gpu.close()
    print(json.dumps(dict(status='pass',device=gpu.device_name,host=platform.node(),python=platform.python_version(),
                          records=results,wall_seconds=time.perf_counter()-start),indent=2,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',action='store_true')
    p.add_argument('phase_source')
    p.add_argument('cut_source')
    p.add_argument('verify_source')
    p.add_argument('metal_source')
    p.add_argument('witness_file')
    a=p.parse_args()
    data=json.loads(Path(a.witness_file).read_text())
    w=data['minima']; extra=data.get('additional_seeds',{})
    paths=(a.phase_source,a.cut_source,a.verify_source,a.metal_source)
    if a.bundle:
        bundle(paths,w,extra)
    else:
        for name,path in zip(('phase_lifts','transverse','verify_witnesses','metal_base'),paths):
            load(name,path)
        run(w,extra)
