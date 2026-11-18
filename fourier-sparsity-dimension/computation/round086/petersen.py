#!/usr/bin/env python3
"""True Petersen-slot quadratic pencils: exact closure cuts on the graph quotient."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
import multiprocessing as mp
import platform
import random
import time

EDGES=tuple([(i,(i+1)%5) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)]+[(i,i+5) for i in range(5)])
KINDS=('star','coordinate','bilinear','two-line','colored-skew','mixed-vertices','star-perturb')
SEED=860828


def parity(x):
    return bin(x).count('1')%2


def rank_basis(rows):
    pivots={}
    for v in rows:
        while v:
            i=v.bit_length()-1
            if i not in pivots:
                pivots[i]=v; break
            v^=pivots[i]
    return pivots


def rank(rows):
    return len(rank_basis(rows))


def xor(a,b):
    return tuple(x^y for x,y in zip(a,b))


def wedge(u,v,a):
    return tuple((v if u>>i&1 else 0)^(u if v>>i&1 else 0) for i in range(a))


def positions(v):
    while v:
        b=v&-v
        yield b.bit_length()-1
        v^=b


def closure(mask):
    labels=list(range(10))
    for i in positions(mask):
        u,v=EDGES[i]; x=labels[u]; y=labels[v]
        if x!=y:
            labels=[x if z==y else z for z in labels]
    k=len(set(labels)); inside=sum((labels[u]==labels[v])<<i for i,(u,v) in enumerate(EDGES))
    return k,inside


@lru_cache(None)
def partitions():
    vals=sorted(set(closure(i) for i in range(1<<15)))
    if any(2*(15-len(list(positions(E))))<3*k for k,E in vals if k>1):
        raise AssertionError('Petersen crossing lower bound failed')
    return tuple((k,E,tuple(positions(E))) for k,E in vals)


def normalize(raw,t,rng=None):
    # Principal compression by pivot coordinates is injective on the joint image.
    W=rank_basis(row for Q in raw for row in Q)
    if not W:
        return None
    coords=sorted(W)
    if rng is not None:
        rng.shuffle(coords)
    rows=[tuple(sum((Q[i]>>j&1)<<k for k,j in enumerate(coords)) for i in coords) for Q in raw]
    if t and len(coords)<t:
        return None
    if any(rank(A)!=rank(B) for A,B in zip(raw,rows)):
        raise AssertionError('compression changed a polar rank')
    return rows,coords


def choose(kind,rng):
    if kind=='star':
        a=11
        VV=[wedge(1,1<<(i+1),a) for i in range(10)]
        return tuple(xor(VV[i],VV[j]) for i,j in EDGES),True
    l,r=rng.randint(2,10),rng.randint(2,10)
    if kind in ('coordinate','bilinear','two-line'):
        a=l+r
        if kind=='coordinate':
            u=[1<<rng.randrange(l) for _ in range(10)]
            v=[1<<(l+rng.randrange(r)) for _ in range(10)]
        elif kind=='bilinear':
            u=[rng.randrange(1,1<<l) for _ in range(10)]
            v=[rng.randrange(1,1<<r)<<l for _ in range(10)]
        else:
            u=[1<<rng.randrange(2) for _ in range(10)]
            v=[(1<<l)|((1<<(l+i%r)) if rng.randrange(3) else 0) for i in range(10)]
        V=[wedge(x,y,a) for x,y in zip(u,v)]
        return tuple(xor(V[i],V[j]) for i,j in EDGES),True
    a=rng.randint(8,18)
    V=[]
    if kind=='colored-skew':
        first=[1<<rng.randrange(a) for _ in range(10)]
        second=[1<<rng.randrange(a) for _ in range(10)]
        V=[wedge(x,y,a) for x,y in zip(first,second)]
    elif kind=='mixed-vertices':
        for j in range(10):
            p,q=rng.randrange(1,1<<a),rng.randrange(1,1<<a)
            V.append(wedge(p,q,a))
            if rng.randrange(6)==0:
                V[-1]=xor(V[-1],wedge(rng.randrange(1,1<<a),rng.randrange(1,1<<a),a))
    elif kind=='star-perturb':
        star=1<<rng.randrange(a)
        V=[wedge(star,rng.randrange(1,1<<a),a) for _ in range(10)]
        for j in rng.sample(range(10),rng.randrange(1,5)):
            V[j]=xor(V[j],wedge(rng.randrange(1,1<<a),rng.randrange(1,1<<a),a))
    else:
        raise ValueError(kind)
    return tuple(xor(V[i],V[j]) for i,j in EDGES),False


def data(M,t):
    a=len(M[0]); mask=(1<<t)-1
    r=tuple(rank(Q) for Q in M)
    p=tuple(rank(v&mask for v in Q) for Q in M)
    if any(x%2 or y>x for x,y in zip(r,p)):
        raise AssertionError('invalid polar rank')
    K=1<<(max(r)//2)
    columns=tuple(K<<(x//2) for x in r)
    rebate=tuple(v>>y for v,y in zip(columns,p))
    return r,p,K,columns,rebate


@lru_cache(maxsize=8192)
def exact_edgebill(d,K):
    # Exact max of 2*internal rebate + 3*K*blocks, over graph-closed partitions.
    classes={}
    for i,value in enumerate(d):
        classes[value]=classes.get(value,0)|(1<<i)
    best=None
    for k,m,_ in partitions():
        score=3*K*k+2*sum(v*len(list(positions(m&B))) for v,B in classes.items())
        if best is None or (score,k,-m)>best[:3]:
            best=(score,k,-m)
    return best[0],best[1],-best[2]


def compact_edgebill(d,K):
    # Equivalent exact optimizer, optimized for the usual small number of distinct rebate weights.
    classes={}
    for i,value in enumerate(d):
        classes[value]=classes.get(value,0)|(1<<i)
    def pop(x):
        return x.bit_count() if hasattr(x,'bit_count') else bin(x).count('1')
    best=(-1,0,0)
    for k,m,_ in partitions():
        val=3*K*k+2*sum(v*pop(m&B) for v,B in classes.items())
        if (val,k,-m)>best:
            best=val,k,-m
    return best[0],best[1],-best[2]


CACHED_bill=lru_cache(maxsize=8192)(compact_edgebill)


def estimate(M,t,full=True):
    a=len(M[0]); r,p,K,c,d=data(M,t); C=sum(c)
    s_d=sorted(d); B=C-sum(d)
    lb=[2*B-3*K*(t+1)]
    for k in range(2,11):
        lower=sum(s_d[:(3*k+1)//2])
        lb.append(2*(B+lower)-3*K*(t+k))
    envelope=min(lb)
    M0=2*C-3*K*(t+10)
    best_hint=min(lb[0],M0)
    if full:
        score,k,U=CACHED_bill(d,K)
        slack=2*C-score-3*K*t
    else:
        score,k,U=None,None,None
        slack=best_hint
    return dict(a=a,t=t,ranks=r,projected=p,scale=2*K,baseline_bound=envelope,k=k,inside=U,
                minimum_margin=[slack,2*K],margin_X0=[M0,2*K],exact_partition=full)


def kernel(rows,n):
    P=rank_basis(rows)
    output=[]
    for j in range(n):
        if j in P:
            continue
        x=1<<j
        for i in sorted(P):
            if parity(P[i]&x):x|=1<<i
        output.append(x)
    return output


def compatible(M,with_extension=False):
    a=len(M[0]); constraints=[]
    for Q,(u,v) in zip(M,EDGES):
        for h in kernel(Q,a):
            constraints.append((h<<((u-1)*a) if u else 0)^(h<<((v-1)*a) if v else 0))
    dim=9*a-rank(constraints)
    if dim<a:
        raise AssertionError('local forms lost')
    record=dict(dimension=dim,local_dimension=a,nonlocal_dimension=dim-a)
    if with_extension and dim>a:
        # Local T represented by its values at e_0+e_i, using path integration of M_e.
        potentials=[None]*10;potentials[0]=(0,)*a
        while any(q is None for q in potentials):
            for Q,(u,v) in zip(M,EDGES):
                if potentials[u] is not None and potentials[v] is None:potentials[v]=xor(potentials[u],Q)
                if potentials[v] is not None and potentials[u] is None:potentials[u]=xor(potentials[v],Q)
        if any(xor(potentials[u],potentials[v])!=Q for Q,(u,v) in zip(M,EDGES)):
            raise AssertionError('matrices are not graph vertex differences')
        local=[sum(potentials[i][q]<<((i-1)*a) for i in range(1,10)) for q in range(a)]
        w=next(z for z in kernel(constraints,9*a) if rank([*local,z])>a)
        ts=[0]+[(w>>(a*i))&((1<<a)-1) for i in range(9)]
        lift=[]
        for Q,(i,j) in zip(M,EDGES):
            d=ts[i]^ts[j]
            lift.append(tuple(q|(((d>>i)&1)<<a) for i,q in enumerate(Q))+(d,))
        if any(rank(Q)!=rank(B) for Q,B in zip(M,lift)) or rank(q for T in lift for q in T)!=a+1:
            raise AssertionError('invalid advertised rank-preserving extension')
        record.update(witness_linear_map=ts,extended_matrices=lift,ranks_unchanged=True)
    return record


def min_witness(k,X,row,M,center,kind):
    key=Fraction(*row['minimum_margin'])
    record=dict(**row,matrices=M,centered=bool(center),kind=kind)
    if k not in X or key<X[k][0]:
        X[k]=(key,record)


def work(task):
    i,count,domain=task; rng=random.Random(SEED+104729*domain+7919*i)
    hist={k:Counter() for k in KINDS}; best={}; nonlocal_best=None; comp=Counter()
    start=time.process_time()
    for j in range(count):
        kind=KINDS[(i+j)%len(KINDS)]
        W,center=choose(kind,rng)
        packed=normalize(W,8,rng)
        if packed is None:
            hist[kind]['below_t8_rejected']+=1
            continue
        M,coords=packed; a=len(coords)
        t=rng.randint(8,a) if rng.randrange(3) else a
        rec=estimate(M,t,False); base=rec['baseline_bound']; denom=rec['scale']
        h=hist[kind]; h['tested']+=1
        h['mixed_r2_r4']+=2 in rec['ranks'] and 4 in rec['ranks']
        # Optimize exactly if the coarse rigorous lower envelope is vulnerable, or at a deterministic sparse sample.
        full=base<3*denom or j%41==0
        if full:
            rec=estimate(M,t,True); h['exact_partition_optimized']+=1
            val=Fraction(*rec['minimum_margin'])
            h['exact_negative']+=val<0
            h['centered_negative']+=val<0 and center
            h['exact_zero']+=val==0
            min_witness('all',best,rec,M,center,kind)
            if 2 in rec['ranks'] and 4 in rec['ranks']:
                min_witness('mixed',best,rec,M,center,kind)
                if center:
                    min_witness('mixed_centered',best,rec,M,center,kind)
        h['baseline_sufficient']+=base>=0
        margin=Fraction(*estimate(M,a,False)['margin_X0']); comp['all_zero_cut_attempts']+=1
        if margin<3 or j%53==0:
            q=compatible(M,False); comp['computed_compatibility']+=1
            comp['nonlocal']+=q['nonlocal_dimension']>0
            comp['margin_below_one']+=margin<1
            comp['can_make_negative_by_one_axis']+=q['nonlocal_dimension']>0 and margin<Fraction(3,2)
            # New transverse axis adds 3/2 in the original mass units for L(G)=2/3.
            if q['nonlocal_dimension'] and (nonlocal_best is None or margin<nonlocal_best[0]):
                nonlocal_best=(margin,dict(margin=[margin.numerator,margin.denominator],old_matrices=M,
                                          a=a,kind=kind,**compatible(M,True)))
    return hist,best,comp,nonlocal_best,time.process_time()-start


def run(workers,tasks,draws,domain):
    begin=time.perf_counter(); partitions()
    hist={k:Counter() for k in KINDS}; comp=Counter(); best={}; nonlocal_best=None; cpu=0
    processing=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as p:
        for h,c,q,r,secs in p.imap_unordered(work,[(i,draws,domain) for i in range(tasks)]):
            for k,v in h.items():hist[k].update(v)
            comp.update(q);cpu+=secs
            for k,(f,x) in c.items():
                if k not in best or f<best[k][0]:best[k]=f,x
            if r and (nonlocal_best is None or r[0]<nonlocal_best[0]):nonlocal_best=r
    elapsed=time.perf_counter()-processing
    return dict(worker_threads=workers,tasks=tasks,draws_per_task=draws,domain=domain,
                total_wall_seconds=time.perf_counter()-begin,workers_wall_seconds=elapsed,
                worker_cpu_seconds=cpu,mean_active_workers=cpu/elapsed,partitions=len(partitions()),
                by_kind=hist,compatibility=comp,best={k:x[1] for k,x in best.items()},
                best_nonlocal=nonlocal_best[1] if nonlocal_best else None)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers',type=int,default=12);p.add_argument('--tasks',type=int,default=48)
    p.add_argument('--draws',type=int,default=512);p.add_argument('--domain',type=int,default=0)
    args=p.parse_args()
    data=run(args.workers,args.tasks,args.draws,args.domain)
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),seed=SEED,**data),
                     indent=2,sort_keys=True))
