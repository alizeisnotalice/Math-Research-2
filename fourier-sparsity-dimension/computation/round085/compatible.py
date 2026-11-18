#!/usr/bin/env python3
"""Search rank-preserving transverse extensions via compatible linear maps."""
import argparse
from collections import Counter
from fractions import Fraction
import json
import multiprocessing as mp
from pathlib import Path
import platform
import random
import time

import phase_lifts as phase
import transverse as tra


def dual_nullspace(rows,dimension):
    piv={u.bit_length()-1:u for u in phase.pivots(rows)}
    free=[]
    for j in range(dimension):
        if j in piv:
            continue
        v=1<<j
        for i,row in sorted(piv.items()):
            if phase.parity(row&v):
                v|=1<<i
        free.append(v)
    return free


def data_matrices(seed,gens,coords,a):
    return [tra.sum_matrices([gens[j] for j in tra.points(q)],a) for q in coords]


def entire_mass(seed,M):
    return Fraction(sum(w*(1<<(phase.rank(rows)//2)) for w,rows in zip(seed['masses'],M)),seed['denominator'])


def compatible_dimension(seed,gens,coords,a,M,extend=True):
    constraints=[]
    R=seed['R']
    for rows,coord in zip(M,coords):
        for h in dual_nullspace(rows,a):
            constraints.append(sum(h<<(j*a) for j in tra.points(coord)))
    K=R*a-phase.rank(constraints)
    if K<a:
        raise AssertionError('compatible linear maps lost the local subspace')
    output=dict(compatible_dimension=K,local_dimension=a,nonlocal_dimension=K-a)
    if K>a and extend:
        local=[sum(m[j]<<(i*a) for i,m in enumerate(gens)) for j in range(a)]
        if phase.rank(local)!=a:
            raise AssertionError('local maps not independently parametrized')
        v=next(v for v in dual_nullspace(constraints,R*a) if phase.rank([*local,v])>a)
        new=[]
        for i,Q in enumerate(gens):
            column=(v>>(i*a))&((1<<a)-1)
            new.append(tuple(row|(((column>>j)&1)<<a) for j,row in enumerate(Q))+(column,))
        if phase.rank(row for Q in new for row in Q)!=a+1:
            raise AssertionError('nonlocal extension does not add address rank')
        test=data_matrices(seed,new,coords,a+1)
        if [phase.rank(m) for m in M]!=[phase.rank(m) for m in test]:
            raise AssertionError('nonlocal extension changed a real rank')
        output.update(extension_vector=v,new_generators=new,certified_ranks_unchanged=True)
    return output


SEEDS=[]


def worker(task):
    i,draws,domain=task
    seed=SEEDS[i%len(SEEDS)]; R=seed['R']; L=Fraction(*seed['L'])
    rng=random.Random(850830+i*65537+104729*domain)
    ctr=Counter(); best=None; best_nonlocal=None; extensions=None; start=time.process_time()
    for j in range(draws):
        kind=tra.KINDS[(i+j)%len(tra.KINDS)]
        a=rng.randint(5,min(9,2*R))
        gen=tra.generators(seed,kind,a,rng)
        if gen is None:
            continue
        mats,coord=gen
        M=data_matrices(seed,mats,coord,a)
        absolute=entire_mass(seed,M)
        bill=absolute*L
        margin=bill-(R+a+1)
        ctr['full_address_trials']+=1
        if margin<1:
            ctr['near_one_axis']+=1
        if margin>=2 and rng.randrange(9):
            continue
        C=compatible_dimension(seed,mats,coord,a,M,extend=False)
        ctr['compatibility_rechecks']+=1
        ctr['nonlocal_found']+=C['nonlocal_dimension']>0
        ctr['possible_one_axis_counterexample']+=C['nonlocal_dimension']>0 and margin<1
        rank_gap=bill-(R+C['compatible_dimension']+1)
        ctr['compatibility_dimension_upper_bound_ties']+=rank_gap==0
        ctr['compatibility_dimension_upper_bound_failures']+=rank_gap<0
        record=dict(seed=seed['name'],a=a,kind=kind,coordinates=coord,generators=mats,L=seed['L'],
                    parent_l1=[absolute.numerator,absolute.denominator],
                    margin=[margin.numerator,margin.denominator],
                    compatibility_dimension_margin=[rank_gap.numerator,rank_gap.denominator],**C)
        if best is None or margin<Fraction(*best['margin']):
            best=record
        if C['nonlocal_dimension'] and (best_nonlocal is None or margin<Fraction(*best_nonlocal['margin'])):
            cov=compatible_dimension(seed,mats,coord,a,M)
            best_nonlocal={**record,**cov}
    return ctr,best,best_nonlocal,time.process_time()-start


def run(scan_ref,full_fixture,workers,draws,domain=0):
    global SEEDS
    start=time.perf_counter(); SEEDS=tra.compute_states()
    cert=[]
    first=scan_ref['minima']['weighted']
    ee=next(s for s in SEEDS if s['name']==first['seed'])
    C=compatible_dimension(ee,first['generators'],first['generator_coords'],first['a'],
                           data_matrices(ee,first['generators'],first['generator_coords'],first['a']))
    cert.append(dict(case='tight-rank2-pencil',a=first['a'],g=first['generators'],**C))
    fixture=full_fixture['exhaustive_rank2_generator_family'][0]['acc']
    rank_four=next(z for z in fixture if z['A_plus_B_rank']==4)
    g=[*rank_four['historical_A_B_C'],tra.zero(5),tra.zero(5)]
    Ms=data_matrices(ee,g,ee['E'],5)
    ex=compatible_dimension(ee,g,ee['E'],5,Ms)
    if ex['compatible_dimension']!=6:
        raise AssertionError('explicit nonlocal extension unexpectedly missing')
    nested=compatible_dimension(ee,ex['new_generators'],ee['E'],6,
                                data_matrices(ee,ex['new_generators'],ee['E'],6))
    cert.append(dict(case='rank4-pencil-before-extension',a=5,g=g,**ex))
    cert.append(dict(case='rank4-pencil-after-extension',a=6,**nested))
    total=Counter(); cpu=0; best=None; nonlocal_best=None
    st=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as p:
        for count,r,nl,secs in p.imap_unordered(worker,[(i,draws,domain) for i in range(4*len(SEEDS))]):
            total.update(count);cpu+=secs
            if r and (best is None or Fraction(*r['margin'])<Fraction(*best['margin'])):
                best=r
            if nl and (nonlocal_best is None or Fraction(*nl['margin'])<Fraction(*nonlocal_best['margin'])):
                nonlocal_best=nl
    elapsed=time.perf_counter()-st
    return dict(workers=workers,domain=domain,draws_per_task=draws,tasks=4*len(SEEDS),wall_seconds=time.perf_counter()-start,
                worker_wall_seconds=elapsed,mean_active_cores=cpu/elapsed,
                exact_examples=cert,random_scan=total,best_margin=best,best_nonlocal_margin=nonlocal_best)


def stream(ph,cut,ref,witness):
    for name,path in (('phase_lifts',ph),('transverse',cut)):
        source=Path(path).read_text()
        print('import types,sys')
        print(f'module=types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}]=module')
        print(f'exec(compile({source!r},{name!r},"exec"),module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__',1)[0])
    print('import argparse')
    print('p=argparse.ArgumentParser();p.add_argument("--workers",type=int,default=15);p.add_argument("--draws",type=int,default=512);p.add_argument("--domain",type=int,default=0);args=p.parse_args()')
    print('summary=run('+repr(ref)+','+repr(witness)+',args.workers,args.draws,args.domain)')
    print('print(json.dumps(dict(status="pass",host=platform.node(),python=platform.python_version(),**summary),indent=2,sort_keys=True))')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',action='store_true')
    p.add_argument('--workers',type=int,default=12)
    p.add_argument('--draws',type=int,default=512)
    p.add_argument('--domain',type=int,default=0)
    p.add_argument('phase_source')
    p.add_argument('cut_source')
    p.add_argument('ref')
    p.add_argument('witness')
    args=p.parse_args()
    ref=json.loads(Path(args.ref).read_text()); witness=json.loads(Path(args.witness).read_text())
    if args.bundle:
        stream(args.phase_source,args.cut_source,ref,witness)
    else:
        summary=run(ref,witness,args.workers,args.draws,args.domain)
        print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),**summary),indent=2,sort_keys=True))
