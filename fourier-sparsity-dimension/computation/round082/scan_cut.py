#!/usr/bin/env python3
"""Exact one-query inequalities and realized P comparisons on Boolean spectra."""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import multiprocessing as mp
import os
import platform
import random
import struct
import time

MASK32 = (1 << 32) - 1
FAMILIES = ('uniform', 'odd', 'sparse', 'quadratic', 'cubic', 'quartic', 'selector')
SEED = 820828


def parity(x):
    return bin(x).count('1') % 2


def fwht(values):
    out = list(values)
    stride = 1
    while stride < len(out):
        for start in range(0,len(out),2*stride):
            for i in range(start,start+stride):
                a,b=out[i],out[i+stride]
                out[i],out[i+stride]=a+b,a-b
        stride *= 2
    return out


def spectrum(truth,n):
    return fwht(1-2*((truth>>u)&1) for u in range(1<<n))


def rank(columns,n):
    pivots=[0]*(n+1)
    r=0
    for column in columns:
        while column:
            i=column.bit_length()-1
            p=pivots[i]
            if p:
                column^=p
            else:
                pivots[i]=column
                r+=1
                break
    return r


def affine_rank(frequencies,n):
    return rank(((1<<n)|u for u in frequencies),n)


def span(gens):
    elements={0}
    for g in gens:
        if g not in elements:
            elements.update(g^w for w in list(elements))
    return elements


def eval_direction(W,n,query,parent_R,parent_s,details=False):
    pivot=1<<(query.bit_length()-1)
    S0=[]
    S1=[]
    union=[]
    k=e=0
    C0={}
    C1={}
    for q in representatives(n,pivot):
        a,b=W[q],W[q^query]
        ap,am=a+b,a-b
        if ap:
            S0.append(q)
            if details:
                C0[q]=ap//2
        if am:
            S1.append(q)
            if details:
                C1[q]=am//2
        if ap or am:
            union.append(q)
        if ap and am:
            k+=1
            e+=abs(ap)==abs(am)
    s0,s1=len(S0),len(S1)
    R0,R1=affine_rank(S0,n),affine_rank(S1,n)
    R_union=affine_rank(union,n)
    j=R0+R1-R_union
    defect=parent_R-R_union
    U=R0+R1+2
    L=s0+s1
    signed_margin=2*L*(2*U*j-j*j)-(2*k+e)*U*U
    actual_margin=2*L*(U*U-(parent_R+1)**2)-(2*k+e)*U*U
    left_delta=(R0+1)**2*parent_s-(parent_R+1)**2*s0
    right_delta=(R1+1)**2*parent_s-(parent_R+1)**2*s1
    if defect not in (0,1) or parent_s != 2*L-2*k-e or not s0 or not s1:
        raise AssertionError((query,'inconsistent cut identity'))
    if (signed_margin>=0 or actual_margin>=0 or k<=j) and max(left_delta,right_delta)<0:
        raise AssertionError((query,'false positive cut proof'))
    facts=(s0,s1,R0,R1,k,j,e,int(signed_margin>=0),int(actual_margin>=0),int(left_delta>=0),int(right_delta>=0))
    if not details:
        return facts,max(left_delta,right_delta),signed_margin
    return dict(query=query,s0=s0,s1=s1,R0=R0,R1=R1,k=k,j=j,e=e,projection_defect=defect,
                signed_margin=signed_margin,sharpened_margin=actual_margin,child_vs_parent_numerators=[left_delta,right_delta],
                support0=S0,support1=S1,walsh0=C0,walsh1=C1)


@lru_cache(maxsize=None)
def representatives(n,pivot):
    return tuple(q for q in range(1<<n) if not q&pivot)


def analyze(truth,n,details=False):
    W=spectrum(truth,n)
    S=[u for u,a in enumerate(W) if a]
    R=affine_rank(S,n)
    difference_space=span(u^S[0] for u in S)
    counts=Counter()
    signed_max=None
    actual_max=None
    signature=0
    all_cuts=[]
    for q in range(1,1<<n):
        facts,actual_delta,margin=eval_direction(W,n,q,R,len(S))
        s0,s1,R0,R1,k,j,e,proof,sharp,nondecrease0,nondecrease1=facts
        for x in facts:
            signature=(signature*1664525+x)&MASK32
        counts['all_directions']+=1
        if q in difference_space:
            counts['effective_directions']+=1
            counts['easy_cuts']+=k<=j
            counts['sufficient_cuts']+=proof
            counts['rank_corrected_cuts']+=sharp
            counts['nondecreasing_cuts']+=nondecrease0 or nondecrease1
            counts['strict_gain_cuts']+=actual_delta>0
            signed_max=margin if signed_max is None else max(signed_max,margin)
            actual_max=actual_delta if actual_max is None else max(actual_max,actual_delta)
        if details:
            case=eval_direction(W,n,q,R,len(S),True)
            case['effective']=q in difference_space
            all_cuts.append(case)
    small=dict(n=n,truth_mask=truth,s=len(S),R=R,signature=signature,counts=dict(counts),
               best_sufficient_margin=signed_max,best_possible_child_delta=actual_max)
    if details:
        small.update(walsh=dict((u,W[u]) for u in S), all_cuts=all_cuts)
    return small


@lru_cache(maxsize=None)
def monomial_masks(n,degree):
    return [sum(1<<x for x in range(1<<n) if x&monomial==monomial)
            for monomial in range(1<<n) if bin(monomial).count('1')<=degree]


def random_truth(n,family,rng):
    size=1<<n
    if family=='uniform':
        return rng.getrandbits(size)
    if family=='odd':
        tail=rng.getrandbits(size-1)
        return tail|((1^parity(tail))<<(size-1))
    if family=='sparse':
        w=rng.randint(2,min(size//2,17))
        return sum(1<<point for point in rng.sample(range(size),w))
    if family in ('quadratic','cubic','quartic'):
        degree=('quadratic','cubic','quartic').index(family)+2
        truth=0
        for m in monomial_masks(n,degree):
            if rng.getrandbits(1):
                truth^=m
        return truth
    if family=='selector':
        left=rng.randint(1,n-2)
        right=n-1-left
        a,b=rng.getrandbits(1<<left),rng.getrandbits(1<<right)
        truth=0
        for x in range(size):
            val=(a>>(x&((1<<left)-1)))&1 if not x>>(n-1) else (b>>((x>>left)&((1<<right)-1)))&1
            truth|=val<<x
        return truth
    raise AssertionError(family)


def fixture_truth(n,fn):
    return sum(int(fn(x))<<x for x in range(1<<n))


def r49(x):
    match=((0,1),(2,4),(3,7),(5,6))
    i=next(i for i,line in enumerate(match) if x&7 in line)
    return x>>(3+i)&1


def weak_fold_1(x):
    match=((0,1),(2,7),(4,3),(6,5))
    i=next(i for i,line in enumerate(match) if x&7 in line)
    return x>>(3+i)&1


def quadratic_bent(x,n):
    return sum(((x>>i)&1)*((x>>(i+1))&1) for i in range(0,n-1,2))%2


def structured_fixtures():
    fixtures=[]
    for n in range(5,9):
        fixtures.append((f'quadratic-rank-{n//2*2}-n{n}',n,fixture_truth(n,lambda x:quadratic_bent(x,n))))
        fixtures.append((f'odd-single-flip-n{n}',n,1))
        fixtures.append((f'cubic-selector-n{n}',n,fixture_truth(n,lambda x:(x&7)==7 or ((x>>(n-1))&1 and (x&3)==0))))
    fixtures.extend([
        ('address-one-tensor-bent',5,fixture_truth(5,lambda x:((x>>(1+(x&1)))&1)^quadratic_bent(x>>3,2))),
        ('address-two',6,fixture_truth(6,lambda x:(x>>(2+(x&3)))&1)),
        ('address-two-tensor-bent',8,fixture_truth(8,lambda x:((x>>(2+(x&3)))&1)^quadratic_bent(x>>6,2))),
        ('r49-height-one',7,fixture_truth(7,r49)),
        ('r49-height-one-modulated',8,fixture_truth(8,lambda x:r49(x)^(x>>7))),
        ('weak-folding-matching',7,fixture_truth(7,weak_fold_1)),
    ])
    return fixtures


def worker(task):
    block,n,count,domain=task
    rng=random.Random(SEED+1000003*domain+7919*block+101*n)
    start=time.process_time()
    summary={f:Counter() for f in FAMILIES}
    best={}
    digest=hashlib.sha256()
    for i in range(count):
        family=FAMILIES[(i+block)%len(FAMILIES)]
        truth=random_truth(n,family,rng)
        record=analyze(truth,n)
        digest.update(struct.pack('<I',record['signature']))
        c=summary[family]
        metrics=record['counts']
        c['tables']+=1
        c['all_directions']+=metrics['all_directions']
        c['effective_directions']+=metrics.get('effective_directions',0)
        if metrics.get('effective_directions',0)==0:
            c['characters']+=1
            continue
        checks={
            'no_easy_cut':not metrics.get('easy_cuts'),
            'all_sufficient_fail':not metrics.get('sufficient_cuts'),
            'all_rank_corrected_fail':not metrics.get('rank_corrected_cuts'),
            'no_strict_gain':not metrics.get('strict_gain_cuts'),
            'all_children_drop':not metrics.get('nondecreasing_cuts'),
        }
        for key,yes in checks.items():
            if yes:
                c[key]+=1
                pair=(family,key)
                quality=((record['R']+1)**2*65536//record['s'],-truth)
                if pair not in best or quality>best[pair][0]:
                    best[pair]=(quality,truth,record['s'],record['R'])
    return dict(id=block,n=n,worker_cpu_seconds=time.process_time()-start,
                metrics={f:dict(c) for f,c in summary.items()},sha256=digest.hexdigest(),
                certificate_candidates=[dict(family=f,test=k,quality=q,n=n,truth=t,s=s,R=r)
                                        for (f,k),(q,t,s,r) in sorted(best.items())])


def parallel(workers,sizes,domain):
    tasks=[]
    block=0
    for n,total in zip(range(5,9),sizes):
        ntask=math_tasks(workers,total)
        for i in range(ntask):
            count=total//ntask+(i<total%ntask)
            tasks.append((block,n,count,domain))
            block+=1
    # Mix dimensions through the queue so all worker cores get large tasks from the start.
    random.Random(SEED).shuffle(tasks)
    started=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        parts=list(pool.imap_unordered(worker,tasks))
    wall=time.perf_counter()-started
    totals={n:{f:Counter() for f in FAMILIES} for n in range(5,9)}
    selected={}
    for p in sorted(parts,key=lambda p:p['id']):
        for f,v in p['metrics'].items():
            totals[p['n']][f].update(v)
        for rec in p['certificate_candidates']:
            key=rec['n'],rec['family'],rec['test']
            if key not in selected or rec['quality']>selected[key]['quality']:
                selected[key]=rec
    cpu=sum(p['worker_cpu_seconds'] for p in parts)
    return dict(workers=workers,n_batches=len(parts),domain=domain,totals=totals,wall_seconds=wall,
                sum_worker_cpu_seconds=cpu,mean_workers_active=cpu/wall,
                batch_ledger=[{k:v for k,v in p.items() if k!='certificate_candidates'} for p in sorted(parts,key=lambda p:p['id'])],
                certificates=[rec for key,rec in sorted(selected.items())])


def math_tasks(workers,total):
    return max(1,min(4*workers,total))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('structured','parallel','witness'))
    parser.add_argument('--workers',type=int,default=os.cpu_count())
    parser.add_argument('--sizes',type=int,nargs=4,default=[10240,8192,6144,4096])
    parser.add_argument('--domain',type=int,default=0)
    parser.add_argument('--truth',type=int)
    parser.add_argument('--bits',type=int)
    args=parser.parse_args()
    if args.mode=='structured':
        stats=[]
        for name,n,truth in structured_fixtures():
            record=analyze(truth,n)
            if not record['counts'].get('nondecreasing_cuts') or not record['counts'].get('sufficient_cuts'):
                record=analyze(truth,n,True)
            stats.append(dict(name=name,**record))
    elif args.mode=='witness':
        stats=analyze(args.truth,args.bits,True)
    else:
        stats=parallel(args.workers,args.sizes,args.domain)
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),seed=SEED,
                          mode=args.mode,records=stats),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
