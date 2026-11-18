#!/usr/bin/env python3
"""Apple Metal true shifted-Boolean Walsh transforms, checked against sparse lifts."""
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
    const uint n=p[0],k=p[1],common=p[2],addr=p[3],factored=p[4];
    const uint bits=1u<<n;
    const uint u_mask=bits-1u;
    const uint beta=task&u_mask,alpha=task>>n;
    const device uint *labels=p+5;
    const device uint *tables=labels+2u*k;
    int result=0;
    for (uint u=0;u<(1u<<addr);++u) {
        uint shift=0u,flip=popcount(alpha&u)&1u;
        for (uint j=0u;j<k;++j) {
            uint gate=common?((u&1u)&((u>>(j+1u))&1u)):(((u>>(2u*j))&1u)&((u>>(2u*j+1u))&1u));
            if (gate) {shift^=labels[2u*j];flip^=labels[2u*j+1u];}
        }
        if (factored) {
            result += (1-2*int(flip^(popcount(beta&shift)&1u)))*as_type<int>(tables[beta]);
            continue;
        }
        for (uint x=0;x<bits;++x) {
            uint index=x^shift;
            uint value=(tables[index/32u]>>(index%32u))&1u;
            result+=1-2*int(value^flip^(popcount(beta&x)&1u));
        }
    }
    out[task]=result;
}
'''


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot import '+name)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)


def bundle(paths,witnesses,extra):
    for name,path in zip(('phase_lifts','metal_base'),paths):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print(f'run({witnesses!r},{extra!r})')


def check_lift(gpu,phase,old,rows,mode,factored=False):
    new,_=phase.lift(old,rows,mode,dense_cap=-1)
    n=old['n']; address=new['n']-n
    data=array('I',[n,len(rows),int(mode=='common'),address,int(factored)])
    for a,b in rows:
        data.extend((a,b))
    truth=phase.truth_from_sparse(old)
    if factored:
        data.extend(w&0xffffffff for w in phase.fwht(truth))
    else:
        for start in range(0,len(truth),32):
            data.append(sum((a==-1)<<j for j,a in enumerate(truth[start:start+32])))
    count=1<<new['n']
    raw,gt,wt=gpu.compute(data,count,1)
    direct=array('i'); direct.frombytes(raw)
    if len(direct)!=count:
        raise AssertionError('wrong Metal length')
    for i,c in enumerate(direct):
        if c*new['denom']!=new['coeff'].get(i,0)*count:
            raise AssertionError((n,mode,rows,i,c,'actual Metal Walsh mismatch'))
    return new,dict(mode=mode,active_quotient_rank=len(rows),previous_ambient=n,coefficient_count=count,
                    **phase.state_summary(new),gpu_command_seconds=gt,commit_and_wait_seconds=wt,
                    factored_address_transform=factored,
                    actual_int32_walsh_sha256=hashlib.sha256(raw).hexdigest())


def rank4_fixture(phase):
    core=7; t=3; k=3; n=t+core+k
    coeff={}
    values=[]
    for x in range(1<<n):
        address=x&7
        if address<4:
            flip=phase.parity((x>>(t+core))&(1,2,4,3)[address])
        else:
            flip=int((x>>t)&((1<<core)-1)==0)
        values.append(1-2*flip)
    coeff={u:v for u,v in enumerate(phase.fwht(values)) if v}
    return dict(n=n,denom=1<<n,coeff=coeff),[(1<<(t+core+i),0) for i in range(k)]


def run(witnesses,extra=False):
    import metal_base
    import phase_lifts as phase
    metal_base.MSL=SOURCE
    start=time.perf_counter(); gpu=metal_base.Metal(); results=[]
    try:
        for name in (() if extra else ('r49','address-two','rare-3','rare-4','nested-r49','bent-four')):
            truth,n=phase.fixture_truth(name)
            seed=phase.normalized_seed(truth,n)
            for choice in phase.greedy_select(seed,kmax=3):
                for mode in ('common','disjoint'):
                    rows=list(choice['rows']); state=seed
                    for step in (1,2,3):
                        required=state['n']+(len(rows)+1 if mode=='common' else 2*len(rows))
                        if required>15 or (step>=2 and choice['k']>2):
                            break
                        n0=state['n']
                        state,record=check_lift(gpu,phase,state,rows,mode)
                        record.update(seed=name,initial_k=choice['k'],stage=step)
                        results.append(record)
                        rows.extend((1<<j,0) for j in range(n0,state['n']))
        if extra:
            seed,rows=rank4_fixture(phase)
            new,output=check_lift(gpu,phase,seed,rows,'disjoint',factored=True)
            output.update(seed='rank-four-sparse-fallback',initial_k=3,stage=1,old=phase.state_summary(seed))
            results.append(output)
            for name in ('rare-4','rare-5'):
                old=phase.normalized_seed(*phase.fixture_truth(name))
                n=old['n']
                rows=phase.greedy_select(old,kmax=1)[0]['rows']
                first,_=phase.lift(old,rows,'disjoint',dense_cap=-1)
                a,b=rows[0]
                rebased=[(a^(1<<n),b),(1<<n,0),(1<<(n+1),0)]
                _,output=check_lift(gpu,phase,first,rebased,'disjoint')
                output.update(seed=name+'-sharp-rebase',initial_k=1,stage=2,old=phase.state_summary(first))
                results.append(output)
        for record in (() if extra else witnesses):
            seed=phase.normalized_seed(record['truth'],record['n'])
            mode='common' if record['case'].endswith('common') else 'disjoint'
            _,output=check_lift(gpu,phase,seed,record['choice']['rows'],mode)
            output.update(seed=record['case'],initial_k=record['choice']['k'],stage=1)
            results.append(output)
    finally:
        gpu.close()
    gpu_seconds=sum(r['gpu_command_seconds'] for r in results)
    wait=sum(r['commit_and_wait_seconds'] for r in results)
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),device=gpu.device_name,
                          backend='Metal integer direct Walsh summation on true shifted Boolean tables',records=results,
                          wall_seconds=time.perf_counter()-start,gpu_command_seconds=gpu_seconds,
                          gpu_command_fraction_of_wait=gpu_seconds/wait),indent=2,sort_keys=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true')
    parser.add_argument('--witness-file')
    parser.add_argument('--extra',action='store_true')
    parser.add_argument('lifts_source')
    parser.add_argument('metal_source')
    args=parser.parse_args()
    paths=(args.lifts_source,args.metal_source)
    witnesses=json.loads(Path(args.witness_file).read_text())['witnesses'] if args.witness_file else []
    if args.bundle:
        bundle(paths,witnesses,args.extra)
    else:
        load('phase_lifts',paths[0]);load('metal_base',paths[1]);run(witnesses,args.extra)
