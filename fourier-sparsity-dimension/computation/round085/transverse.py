#!/usr/bin/env python3
"""Exact cut-mass diagnostics for quadratic translation lifts of true Boolean spectra."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
import multiprocessing as mp
from pathlib import Path
import platform
import random
import time

import phase_lifts as phase

SEED=850828
KINDS=('basis-rank2','random-rank2','rank4-and-two','holes','disjoint','generic')


def pop(x):
    return x.bit_count() if hasattr(x,'bit_count') else bin(x).count('1')


def points(mask):
    while mask:
        low=mask&-mask
        yield low.bit_length()-1
        mask^=low


def span(rows):
    result=[0]
    for u in rows:
        result.extend(u^v for v in result[:])
    return result


def xor_all(items):
    a=0
    for b in items:
        a^=b
    return a


@lru_cache(None)
def subspaces(R):
    total=(1<<(1<<R))-1
    current={1}; output=[]
    for k in range(R+1):
        output.append(tuple(sorted(current)))
        if k==R:
            break
        nxt=set()
        for H in current:
            pointsH=list(points(H)); remaining=total^H
            while remaining:
                bit=remaining&-remaining
                u=bit.bit_length()-1
                other=sum(1<<(v^u) for v in pointsH)
                remaining^=other
                nxt.add(H|other)
        current=nxt
    return output


def exact_L(state):
    R=state['n']+1
    if R>7:
        raise ValueError('exhaustive base L limited to seven homogeneous dimensions')
    classes=Counter()
    for b,co in state['coeff'].items():
        classes[abs(co)]|=1<<((1<<state['n'])|b)
    norm=sum(abs(z) for z in state['coeff'].values())
    scale=state['denom']
    best=Fraction(0); tight=[]
    for k,level in enumerate(subspaces(R)[:-1]):
        for U in level:
            mass=norm-sum(co*pop(U&S) for co,S in classes.items())
            ratio=Fraction((R-k+1)*scale,mass)
            if ratio>best:
                best=ratio; tight=[(R-k,U)]
            elif ratio==best:
                tight.append((R-k,U))
    return best,tight


def compute_states():
    states=[]
    for name in ('address-two','r49','rare-2','rare-3','rare-4','bent-four'):
        raw=phase.fixture_truth(name)
        states.append((name,phase.normalized_seed(*raw),raw))
    rng=random.Random(SEED+23)
    for n in (3,4,5,6):
        for i in range(3):
            t=phase.random_truth(n,('uniform','sparse','cubic')[i],rng)
            state=phase.normalized_seed(t,n)
            if state['n']>=2:
                states.append((f'actual-{n}-{i}',state,(t,n)))
    results=[]
    for name,state,raw in states:
        L,tight=exact_L(state)
        freq=sorted(state['coeff'])
        R=state['n']+1
        E=[b|(1<<(R-1)) for b in freq]
        W=[abs(state['coeff'][b]) for b in freq]
        cheap=sorted(range(len(W)),key=lambda i:(W[i],E[i]))
        red=[]; picked=[]
        for i in cheap:
            if phase.rank([*red,E[i]])>len(red):
                red.append(E[i]); picked.append(i)
        if len(red)!=R:
            raise AssertionError('true seed is not full in homogeneous coordinates')
        redcoords=phase.quotient_coords(red)
        results.append(dict(name=name,parent=state,n=R-1,R=R,truth_source=raw,coefficients=state['coeff'],
                            denominator=state['denom'],L=(L.numerator,L.denominator),optimal_subspaces=tight,
                            E=E,masses=W,cheap_basis_E=red,cheap_coordinates=[redcoords[e] for e in E]))
    return results


def zero(a):
    return (0,)*a


def wedge(x,y,a):
    return tuple((y if (x>>i)&1 else 0)^(x if (y>>i)&1 else 0) for i in range(a))


def sum_matrices(mats,a):
    if not mats:
        return zero(a)
    return tuple(xor_all(rows) for rows in zip(*mats))


def rand_wedge(a,rng):
    u,v=rng.randrange(1,1<<a),rng.randrange(1,1<<a)
    while u==v:
        v=rng.randrange(1,1<<a)
    return wedge(u,v,a)


def generators(seed,kind,a,rng):
    R=seed['R']
    E=seed['E']
    if kind=='basis-rank2':
        output=[rand_wedge(a,rng) if rng.randrange(8) else zero(a) for _ in range(R)]
        coord=seed['cheap_coordinates']
    elif kind=='random-rank2':
        output=[rand_wedge(a,rng) if rng.randrange(8) else zero(a) for _ in range(R)]
        coord=E
    elif kind=='rank4-and-two':
        keep=(1<<(a-2))-1
        output=[]
        for j in range(R):
            Q=rand_wedge(a,rng)
            if rng.randrange(3)==0:
                u=rng.randrange(1,1<<a)|1<<(a-1)
                v=rng.randrange(1,1<<a)|1<<(a-2)
                Q=sum_matrices([Q,wedge(u,v,a)],a)
            output.append(Q)
        coord=seed['cheap_coordinates']
    elif kind=='holes':
        # Concentrate actual candidate nonzero matrices on affine cocircuit tails.
        output=[zero(a)]*R
        q=max(2,a//2)
        for i in range(min(q,R)):
            output[i]=rand_wedge(a,rng)
        coord=seed['cheap_coordinates']
    elif kind=='disjoint':
        output=[zero(a)]*R
        for j in range(a//2):
            if j>=R:
                break
            output[j]=wedge(1<<(2*j),1<<(2*j+1),a)
        if a%2 and R>a//2:
            output[a//2]=wedge(1<<(a-1),1<<rng.randrange(a-1),a)
        coord=seed['cheap_coordinates']
    elif kind=='generic':
        output=[]
        pairs=[(i,j) for i in range(a) for j in range(i+1,a)]
        for _ in range(R):
            M=[0]*a
            for i,j in pairs:
                if rng.getrandbits(1):
                    M[i]^=1<<j; M[j]^=1<<i
            output.append(tuple(M))
        coord=E
    else:
        raise ValueError(kind)
    if phase.rank(c for M in output for c in M)!=a:
        return None
    return output,coord


def costs(seed,Ms,coords,a,t,U):
    R=seed['R']; K=1<<(a//2); mask=(1<<t)-1
    values=[]; coeffs=[]; n01=0; weighted=0; sum_p=0
    for e,w,c in zip(seed['E'],seed['masses'],coords):
        rows=sum_matrices([Ms[j] for j in points(c)],a)
        rank=phase.rank(rows)
        p=phase.rank(row&mask for row in rows)
        if rank%2 or p>rank:
            raise AssertionError('invalid polar/projected rank')
        out=(U>>e)&1==0
        term=K<<(rank//2)
        cost=term if out else term-(term>>p)
        weighted+=w*cost
        values.append((cost,e,p,rank,out))
    u=R-(pop(U).bit_length()-1)
    target=u+t+1
    rows=[]; gb=[]; rank_sum=0
    for c,e,p,r,b in sorted(values):
        if phase.rank([*rows,e])>len(rows):
            rows.append(e); gb.append((c,e,p,r,b)); rank_sum+=p
            if len(rows)==R:
                break
    basis_bill=sum(x[0] for x in gb)+max(x[0] for x in gb)
    L_num,L_den=seed['L']
    return dict(a=a,t=t,u=u,U=U,scale=K,target=target,weighted_bill=weighted,basis_bill=basis_bill,
                exact_ratio=[weighted*L_num, target*K*seed['denominator']*L_den],
                greedy_ratio=[basis_bill,target*K],
                all_ranks=sorted(Counter((r,p) for c,e,p,r,b in values).items()),
                projected_rank_sum=rank_sum,rank2_count=sum(x[3]==2 for x in values),
                basis_ranks=[x[3] for x in gb],value_list=values)


def enroll(item,seed,gen,coords,kind):
    pick={k:v for k,v in item.items() if k!='value_list'}
    pick.update(seed=seed['name'],generators=gen,generator_coords=coords,kind=kind)
    return pick


STATES=[]


def worker(task):
    n,ntasks,domain=task
    seed=STATES[n%len(STATES)]
    rng=random.Random(SEED+104729*domain+65537*n)
    hist={k:Counter() for k in KINDS}
    best={}; start=time.process_time()
    for i in range(ntasks):
        kind=KINDS[(i+n)%len(KINDS)]
        a=rng.randint(5,min(11,2*seed['R']))
        t=rng.randint(5,a)
        res=generators(seed,kind,a,rng)
        if res is None:
            hist[kind]['nonfull_address_rejected']+=1
            continue
        M,coords=res
        j=rng.randrange(4)
        if j==0:
            u,U=rng.choice(seed['optimal_subspaces'])
        elif j==1:
            U=(1<<(1<<seed['R']))-1
        else:
            d=rng.randrange(seed['R']+1)
            U=rng.choice(subspaces(seed['R'])[d])
        result=costs(seed,M,coords,a,t,U)
        e0,e1=result['exact_ratio']; g0,g1=result['greedy_ratio']
        h=hist[kind]
        h['verified_cuts']+=1
        h['with_rank2']+=bool(result['rank2_count'])
        h['greedy_insufficient']+=g0<g1
        h['weighted_insufficient']+=e0<e1
        h['weighted_ties']+=e0==e1
        for key,rat in [('weighted',Fraction(e0,e1)),('basis',Fraction(g0,g1))]:
            if key not in best or rat<best[key][0]:
                best[key]=(rat,enroll(result,seed,M,coords,kind))
        if e0<e1 and 'any_weighted_failure' not in best:
            best['any_weighted_failure']=(Fraction(e0,e1),enroll(result,seed,M,coords,kind))
    return seed['name'],hist,best,time.process_time()-start


def run(workers,draws,domain):
    global STATES
    start=time.perf_counter(); STATES=compute_states()
    tasks=[]
    for n in range(len(STATES)*4):
        tasks.append((n,draws,domain))
    hist={k:Counter() for k in KINDS}; per_seed={s['name']:Counter() for s in STATES}
    best={}; seconds=0
    wall=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        for seed,h,c,cpu in pool.imap_unordered(worker,tasks):
            for k,o in h.items():
                hist[k].update(o); per_seed[seed].update(o)
            for key,v in c.items():
                if key not in best or v[0]<best[key][0]:
                    best[key]=v
            seconds+=cpu
    elapsed=time.perf_counter()-wall
    return dict(domain=domain,workers=workers,draws_per_task=draws,ntasks=len(tasks),
                wall_seconds=time.perf_counter()-start,workers_wall_seconds=elapsed,mean_active_workers=seconds/elapsed,
                seeds=[{k:z[k] for k in ('name','n','R','L','truth_source','denominator','coefficients','optimal_subspaces')} for z in STATES],
                by_kind=hist,by_seed=per_seed,minima={k:v[1] for k,v in best.items()})


def bundle(source):
    print('import types,sys')
    print('module = types.ModuleType("phase_lifts")')
    print('sys.modules["phase_lifts"] = module')
    print('exec(compile(' + repr(Path(source).read_text()) + ',"phase_lifts","exec"),module.__dict__)')
    print(Path(__file__).read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('scan','bundle'))
    p.add_argument('--phase-source')
    p.add_argument('--workers',type=int,default=12)
    p.add_argument('--domain',type=int,default=0)
    p.add_argument('--draws',type=int,default=512)
    args=p.parse_args()
    if args.mode=='bundle':
        bundle(args.phase_source)
    else:
        result=run(args.workers,args.draws,args.domain)
        print(json.dumps(dict(status='pass',python=platform.python_version(),host=platform.node(),seed=SEED,**result),indent=2,sort_keys=True))
