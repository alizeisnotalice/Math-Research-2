#!/usr/bin/env python3
"""Exact finite and recursive spectra for affine-selected quadratic shifts."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
import multiprocessing as mp
import os
import platform
import random
import time

SEED=840828
FAMILIES=('uniform','sparse','cubic','quartic','rare','selector')


def parity(x):
    return bin(x).count('1')%2


def pop(x):
    return x.bit_count() if hasattr(x,'bit_count') else bin(x).count('1')


def fwht(x):
    out=list(x); w=1
    while w<len(out):
        for a in range(0,len(out),2*w):
            for j in range(a,a+w):
                v,u=out[j],out[j+w]
                out[j],out[j+w]=v+u,v-u
        w*=2
    return out


def pivots(values):
    d={}
    for v in values:
        while v:
            pivot=v.bit_length()-1
            if pivot in d:
                v^=d[pivot]
            else:
                d[pivot]=v
                break
    return tuple(d.values())


def rank(values):
    return len(pivots(values))


def quotient_coords(basis):
    lookup={0:0}
    for i,row in enumerate(basis):
        lookup.update({v^row:u^(1<<i) for v,u in list(lookup.items())})
    return lookup


def normalized_seed(truth,n):
    W=fwht(1-2*((truth>>x)&1) for x in range(1<<n))
    support=[u for u,w in enumerate(W) if w]
    first=support[0]
    rows=pivots(u^first for u in support)
    coords=quotient_coords(rows)
    return dict(n=len(rows),denom=1<<n,coeff={coords[u^first]:W[u] for u in support})


def truth_from_sparse(state):
    n=state['n']
    dense=fwht(state['coeff'].get(i,0) for i in range(1<<n))
    if any(abs(x)!=state['denom'] for x in dense):
        raise AssertionError('sparse state does not give a Boolean truth table')
    return [v//state['denom'] for v in dense]


def code_masks(S,n,limit,seed):
    axis=[sum((u>>j&1)<<i for i,u in enumerate(S)) for j in range(n)]
    if n<=limit:
        masks=[0]
        for col in axis:
            masks.extend(m^col for m in masks[:])
        forms=list(range(1,1<<n))
        return forms,masks,True
    rng=random.Random(seed+n)
    forms=sorted({1<<j for j in range(n)}|{rng.randrange(1,1<<n) for _ in range(4096)})
    return forms,{f:sum(parity(f&u)<<i for i,u in enumerate(S)) for f in forms},False


def evaluate_affines(rows,point):
    return sum((parity(a&point)^b)<<i for i,(a,b) in enumerate(rows))


def select_label_basis(S,rows,n):
    k=len(rows)
    anchors=[u for u in S if evaluate_affines(rows,u)==0]
    if not anchors:
        raise AssertionError('zero flat is not incident to support')
    fixed=anchors[0]
    origin_rows=[]; relative=[]
    for u in S:
        v=evaluate_affines(rows,u)
        if rank([*relative,v])>len(relative):
            origin_rows.append(u)
            relative.append(v)
            if len(relative)==k:
                break
    code=quotient_coords(relative)
    # Dual rows are adjusted so the reported selected true frequencies have unit phases.
    adjusted=[]
    for i in range(k):
        dual=sum(((code[1<<j]>>i)&1)<<j for j in range(k))
        linear=constant=0
        for j,(a,b) in enumerate(rows):
            if dual>>j&1:
                linear^=a;constant^=b
        adjusted.append((linear,constant))
    if [evaluate_affines(adjusted,u) for u in origin_rows]!=[1<<i for i in range(k)]:
        raise AssertionError('false affine-basis phases')
    # Complete the selected true labels to an affine basis of the entire spectrum.
    complete=[fixed,*origin_rows]
    for u in S:
        dirs=[x^fixed for x in complete]
        if rank([*dirs,u^fixed])>len(complete)-1:
            complete.append(u)
            if len(complete)==n+1:
                break
    if len(complete)!=n+1:
        raise AssertionError('incomplete true spectral affine basis')
    return adjusted,complete


def greedy_select(state,kmax=5,complete_limit=13):
    S=sorted(state['coeff']); n=state['n']; N=len(S)
    if not n:
        return []
    all_points=(1<<N)-1
    forms,masks,complete=code_masks(S,n,complete_limit,SEED)
    rows=[]; linear_space={0}; Z=all_points; out=[]
    for k in range(1,min(n,kmax)+1):
        best=None
        for a in forms:
            if a in linear_space:
                continue
            m=masks[a]
            for const in (0,1):
                vanish=Z&(m if const else ~m)
                w=pop(vanish)
                if w and (best is None or (w,-a,-const)>best[:3]):
                    best=(w,-a,-const,vanish)
        if best is None:
            break
        s0,nega,negc,Z=best
        a,c=-nega,-negc
        rows.append((a,c)); linear_space.update({u^a for u in linear_space})
        adjusted,aff_basis=select_label_basis(S,rows,n)
        weights=Counter(pop(evaluate_affines(adjusted,u)) for u in S)
        common_s=N+3*(N-s0)
        disjoint_s=sum(c*(1<<(2*w)) for w,c in weights.items())
        out.append(dict(k=k,zero=s0,active=N-s0,rows=adjusted,selected_affine_basis=aff_basis,
                        disjoint_weight_histogram=dict(weights),
                        common_p2=frac((n+k+3)**2,common_s),disjoint_p2=frac((n+2*k+2)**2,disjoint_s),
                        all_candidates_per_greedy_step=complete))
    return out


def frac(x,y):
    f=Fraction(x,y)
    return [f.numerator,f.denominator]


def state_summary(state):
    n=state['n']; s=len(state['coeff']); D=n+2
    return dict(v=n,s=s,D=D,P_squared=frac(D*D,s),denom=state['denom'])


def phase_data(state,rows,mode):
    n=state['n']; k=len(rows)
    labels={u:evaluate_affines(rows,u) for u in state['coeff']}
    active=sum(s!=0 for s in labels.values())
    if not active:
        raise AssertionError('trivial selected affine quotient')
    top=1 if mode=='common' else max(pop(x) for x in labels.values())
    addr=k+1 if mode=='common' else 2*k
    expected=len(labels)+3*active if mode=='common' else sum(1<<(2*pop(s)) for s in labels.values())
    return labels, active, top, addr, expected


def branches(bits,mode):
    if mode=='common':
        if not bits:
            return [(0,1)]
        base=bits<<1
        return [(0,1),(1,1),(base,1),(base|1,-1)]
    result=[(0,1)]
    for j in range(bits.bit_length()):
        if bits>>j&1:
            offset=2*j
            result=[(v|(u<<offset),sgn*(1-2*(u==3))) for v,sgn in result for u in range(4)]
    return result


def lift(state,rows,mode,dense_cap=16):
    n=state['n']; k=len(rows)
    labels,active,top,a,expected=phase_data(state,rows,mode)
    cf={}
    for u,c in state['coeff'].items():
        value=labels[u]
        power=bool(value) if mode=='common' else pop(value)
        for addr,p in branches(value,mode):
            cf[u|(addr<<n)]=c*p*(1<<(top-power))
    v=n+a
    new=dict(n=v,coeff=cf,denom=state['denom']<<top)
    if len(cf)!=expected or rank((1<<v)|u for u in cf)!=v+1 or sum(c*c for c in cf.values())!=new['denom']**2:
        raise AssertionError('sparse quadratic lift identities failed')
    dense_checked=False
    if v<=dense_cap:
        old=truth_from_sparse(state)
        true_values=[]
        for z in range(1<<a):
            shift=flip=0
            for j,(p,c) in enumerate(rows):
                gate=(z&1)*(z>>(j+1)&1) if mode=='common' else (z>>(2*j)&1)*(z>>(2*j+1)&1)
                if gate:
                    shift^=p;flip^=c
            true_values.extend((1-2*flip)*old[x^shift] for x in range(1<<n))
        W=fwht(true_values)
        if any(c*new['denom']!=cf.get(u,0)*(1<<v) for u,c in enumerate(W)):
            raise AssertionError('independent actual Boolean shift transform disagrees')
        dense_checked=True
    return new,dense_checked


def follow_recursive(seed,opening,mode,steps=5,cap=4096):
    state=seed
    rows=[tuple(z) for z in opening['rows']]
    split=Counter(pop(evaluate_affines(rows,u)) for u in seed['coeff'])
    base=state_summary(seed); k0=len(rows); m0=opening['active']
    trace=[dict(stage=0,dense_checked=True,**base)]
    x=1
    for t in range(1,steps+1):
        if mode=='common':
            s=base['s']-m0+4**t*m0
            D=base['D']+(2**t-1)*(k0+1)
            m=4**t*m0
            k=2**t*(k0+1)-1
        else:
            x=x*(1+x)**2
            s=sum(c*x**w for w,c in split.items())
            D=base['D']+k0*(3**t-1)
            m=s-split.get(0,0)
            k=k0*3**t
        item=dict(stage=t,s=s,D=D,v=D-2,P_squared=frac(D*D,s),active_tail=m,tail_codim=k,
                  explicit_sparse_checked=False,dense_checked=False)
        if state is not None and s<=cap:
            old_n=state['n']
            state,dense=lift(state,rows,mode)
            if len(state['coeff'])!=s or state['n']+2!=D or sum(evaluate_affines(rows,u)==0 and u>>old_n==0 for u in state['coeff'])!=base['s']-m0:
                raise AssertionError('same-tail closed recurrence disagrees')
            rows.extend((1<<j,0) for j in range(old_n,state['n']))
            if len(rows)!=k or sum(evaluate_affines(rows,u)!=0 for u in state['coeff'])!=m:
                raise AssertionError('actual lifted surviving tail has wrong size or quotient')
            item['explicit_sparse_checked']=True
            item['dense_checked']=dense
        else:
            state=None
        trace.append(item)
    return trace


def greedy_adaptive(seed,steps=5,cap=4096):
    state=seed
    output=[dict(stage=0,**state_summary(state))]
    for t in range(1,steps+1):
        options=greedy_select(state,kmax=6)
        if not options:
            break
        best=max(options,key=lambda record:Fraction(*record['common_p2']))
        if len(state['coeff'])+3*best['active']>cap or state['n']+best['k']+1>22:
            output.append(dict(stage=t,next_test_only=best,reason='explicit cap'))
            break
        old_case=state_summary(state)
        state,dense=lift(state,best['rows'],'common')
        output.append(dict(stage=t,choice=best,strict_increase=Fraction(*best['common_p2'])>Fraction(*old_case['P_squared']),
                           dense_checked=dense,**state_summary(state)))
    return output


def fixture_truth(name):
    if name=='r49':
        n=6; matches=((0,1),(2,4),(3,7),(5,6))
        def f(x):
            i=next(i for i,e in enumerate(matches) if (x&7) in e)
            return 0 if i==0 else x>>(2+i)&1
    elif name=='address-two':
        n=5
        def f(x):
            i=x&3
            return 0 if i==0 else x>>(i+1)&1
    elif name.startswith('rare-'):
        leaf=int(name.split('-')[-1]); n=leaf+2
        def f(x):
            return x>>leaf&1 if x>>(leaf+1) else int(x&((1<<leaf)-1)==0)
    elif name=='bent-four':
        n=4
        def f(x):
            return (x&1)*(x>>1&1)^(x>>2&1)*(x>>3&1)
    elif name=='nested-r49':
        source,v=fixture_truth('r49'); n=v+2
        def f(x):
            return x>>v&1 if x>>(v+1) else source>>(x&((1<<v)-1))&1
    else:
        raise KeyError(name)
    return sum(f(x)<<x for x in range(1<<n)),n


@lru_cache(maxsize=None)
def monomial_masks(n,degree):
    return [sum(1<<x for x in range(1<<n) if x&m==m) for m in range(1<<n) if pop(m)<=degree]


def random_truth(n,family,rng):
    size=1<<n
    if family=='uniform':
        return rng.getrandbits(size)
    if family=='sparse':
        return sum(1<<p for p in rng.sample(range(size),rng.randrange(2,min(size//2,18))))
    if family in ('cubic','quartic'):
        truth=0
        for mask in monomial_masks(n,3+(family=='quartic')):
            if rng.getrandbits(1):
                truth^=mask
        return truth
    if family in ('rare','selector'):
        # All branches are actual truth tables. 'rare' reserves one independent character.
        if family=='rare':
            leaf=rng.getrandbits(1<<(n-2)); v=(1<<(1<<(n-2)))-1
            return leaf|(leaf<<(1<<(n-2)))|(v<<(3*(1<<(n-2))))
        left=rng.getrandbits(1<<(n-1)); right=rng.getrandbits(1<<(n-1))
        return left|(right<<(1<<(n-1)))
    raise KeyError(family)


def scan_worker(task):
    block,n,n_tables,domain=task
    rng=random.Random(SEED+100003*domain+7919*block+101*n)
    start=time.process_time()
    results={f:Counter() for f in FAMILIES}; chosen={}
    for i in range(n_tables):
        f=FAMILIES[(i+block)%len(FAMILIES)]
        truth=random_truth(n,f,rng)
        state=normalized_seed(truth,n)
        summary=state_summary(state)
        old=Fraction(*summary['P_squared'])
        c=results[f]; c['samples']+=1
        c['nonfull']+=summary['s'] < (1<<state['n'])
        results_f=greedy_select(state,kmax=5)
        for mode in ('common','disjoint'):
            if not results_f:
                continue
            key=mode+'_p2'
            choice=max(results_f,key=lambda z:Fraction(*z[key]))
            new=Fraction(*choice[key])
            c[mode+'_increase']+=new>old
            c[mode+'_increase_rank_at_least_two']+=any(Fraction(*x[key])>old and x['k']>=2 for x in results_f)
            c[mode+'_increase_rank_at_least_three']+=any(Fraction(*x[key])>old and x['k']>=3 for x in results_f)
            if new>old and ((f,mode) not in chosen or new/old>Fraction(*chosen[f,mode]['gain_squared'])):
                q=(new/old)
                chosen[f,mode]=dict(seed_truth=truth,n=n,family=f,mode=mode,seed=summary,choice=choice,
                                    gain_squared=[q.numerator,q.denominator])
    return n,results,list(chosen.values()),time.process_time()-start


def sample_run(workers,domain,samples):
    start=time.perf_counter(); tasks=[]; block=0
    for n in range(5,9):
        chunks=min(samples,4*workers)
        for i in range(chunks):
            tasks.append((block,n,samples//chunks+(i<samples%chunks),domain));block+=1
    random.Random(SEED).shuffle(tasks)
    results={n:{f:Counter() for f in FAMILIES} for n in range(5,9)}
    chosen={}; cpu=0
    with mp.get_context('fork').Pool(workers) as pool:
        for n,hist,records,time_s in pool.imap_unordered(scan_worker,tasks):
            cpu+=time_s
            for f,d in hist.items():
                results[n][f].update(d)
            for r in records:
                key=n,r['family'],r['mode']
                if key not in chosen or Fraction(*r['gain_squared'])>Fraction(*chosen[key]['gain_squared']):
                    chosen[key]=r
    elapsed=time.perf_counter()-start
    return dict(samples_per_dimension=samples,workers=workers,domain=domain,wall_seconds=elapsed,
                aggregate_cpu_seconds=cpu,mean_active_cores=cpu/elapsed,by_family=results,
                witnesses=[r for k,r in sorted(chosen.items())])


def run_fixtures():
    start=time.perf_counter(); output=[]
    for name in ('address-two','r49','rare-2','rare-3','rare-4','rare-5','nested-r49','bent-four'):
        truth,n=fixture_truth(name)
        state=normalized_seed(truth,n)
        features=greedy_select(state)
        group=[]
        for choice in features:
            group.append(dict(k=choice['k'],common=follow_recursive(state,choice,'common'),
                              disjoint=follow_recursive(state,choice,'disjoint')))
        output.append(dict(name=name,seed_truth=truth,seed_ambient_dimension=n,seed=state_summary(state),
                           affine_selected=features,iterates=group,adaptive=greedy_adaptive(state)))
    return dict(wall_seconds=time.perf_counter()-start,fixtures=output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('fixtures','samples'))
    p.add_argument('--workers',type=int,default=os.cpu_count())
    p.add_argument('--domain',type=int,default=0)
    p.add_argument('--samples',type=int,default=2048)
    args=p.parse_args()
    result=run_fixtures() if args.mode=='fixtures' else sample_run(args.workers,args.domain,args.samples)
    print(json.dumps(dict(status='pass',python=platform.python_version(),host=platform.node(),mode=args.mode,
                          seed=SEED,**result),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
