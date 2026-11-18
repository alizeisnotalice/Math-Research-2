#!/usr/bin/env python3
"""Exact all one-/two-query cut census on small guarded matching selectors."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
import multiprocessing as mp
import os
import platform
import time

MATCHES={
    'r49':((0,1),(2,4),(3,7),(5,6)),
    'repeated':((0,1),(2,3),(4,6),(5,7)),
}
MODES=('character','global-bent','private-bent','private-guarded-bent','private-guarded-odd')


def parity(x):
    return bin(x).count('1')%2


def bit_count(x):
    return x.bit_count() if hasattr(x,'bit_count') else bin(x).count('1')


def basis(values,n):
    pivots=[0]*(n+1)
    for u in values:
        while u:
            i=u.bit_length()-1
            if pivots[i]:
                u^=pivots[i]
            else:
                pivots[i]=u
                break
    return tuple(u for u in pivots if u)


def span(rows):
    output=[0]
    for r in rows:
        output.extend(u^r for u in output[:])
    return output


def dimension(values,n):
    return len(basis(values,n))


def affine_directions(freq,n):
    first=next(iter(freq))
    return basis((u^first for u in freq),n)


def affine_rank(freq,n):
    return 1+len(affine_directions(freq,n))


def gf_binomial(d,r):
    if not 0<=r<=d:
        return 0
    numerator=denominator=1
    for j in range(r):
        numerator*=(1<<d)-(1<<j)
        denominator*=(1<<r)-(1<<j)
    return numerator//denominator


def type_count(n,ad,h,a,b):
    # Choose the coordinate intersections, then the graph of an isomorphism
    # between residual address and data subspaces.
    m=h-a-b
    bd=n-ad
    gl=1
    for j in range(m):
        gl*=(1<<m)-(1<<j)
    return gf_binomial(ad,a)*gf_binomial(bd,b)*gf_binomial(ad-a,m)*gf_binomial(bd-b,m)*gl


def fwht(values):
    out=list(values)
    w=1
    while w<len(out):
        for pos in range(0,len(out),2*w):
            for i in range(pos,pos+w):
                a,b=out[i],out[i+w]
                out[i],out[i+w]=a+b,a-b
        w*=2
    return out


def seed_phase(z,g):
    if g==2:
        return 1-2*((z&1)*((z>>1)&1))
    if g==3:
        return -1 if z==0 else 1
    return 1


def build_fixture(name):
    match,mode=name.split(':')
    cells=MATCHES[match]
    ad=3
    b=N=4
    guarded=mode in ('character','global-bent','private-guarded-bent','private-guarded-odd')
    g=3 if mode=='private-guarded-odd' else (0 if mode=='character' else 2)
    guards=[0]+[1<<(ad+j) for j in range(N-1)] if guarded else [0]*N
    data_start=ad+(N-1 if guarded else 0)
    shared=mode=='global-bent'
    effective_n=data_start+(g if shared else g*N)
    shift=[data_start if shared else data_start+i*g for i in range(N)]
    Wg=fwht(seed_phase(z,g) for z in range(1<<g))
    coeff=Counter()
    denom=N*(1<<g)
    for i,(x,y) in enumerate(cells):
        for a in range(1<<ad):
            if parity(a&(x^y)):
                continue
            address_sign=1-2*parity(a&x)
            for z,seed_coeff in enumerate(Wg):
                if seed_coeff:
                    coeff[a|guards[i]|(z<<shift[i])]+=address_sign*seed_coeff
    coeff={u:v for u,v in sorted(coeff.items()) if v}
    if sum(v*v for v in coeff.values())!=denom**2:
        raise AssertionError(('Walsh energy',name))
    R=affine_rank(coeff,effective_n)
    if R!=effective_n+1:
        raise AssertionError(('intrinsic coordinates are not full',name,R,effective_n))
    return dict(name=name,n=effective_n,address_dim=ad,coeff=coeff,denom=denom,R=R,guards=guards,
                seed_bits=g,seed_shifts=shift,match=cells,mode=mode,shared=shared)


def actual_truth(case):
    label={x:i for i,line in enumerate(case['match']) for x in line}
    seed_mask=(1<<case['seed_bits'])-1
    for x in range(1<<case['n']):
        i=label[x&7]
        sign=1-2*parity(case['guards'][i]&x)
        phase=seed_phase((x>>case['seed_shifts'][i])&seed_mask,case['seed_bits'])
        yield sign*phase


def verify_fixture(case):
    n=case['n']
    raw=fwht(actual_truth(case))
    for u,value in enumerate(raw):
        if value*case['denom']!=case['coeff'].get(u,0)*(1<<n):
            raise AssertionError((case['name'],u,'dense true-Walsh mismatch'))


def directions(case):
    return sorted({x^y for x,y in combinations(case['coeff'],2)})


def queries_type(rows,ad,n):
    h=len(rows)
    ma=(1<<ad)-1
    return h-dimension((r>>ad for r in rows),n),h-dimension((r&ma for r in rows),n)


def restrictions(case,rows):
    h=len(rows)
    translations=[0,rows[0]] if h==1 else [0,rows[0],rows[1],rows[0]^rows[1]]
    bins={}
    for frequency,value in case['coeff'].items():
        q,t=min((frequency^shift,code) for code,shift in enumerate(translations))
        if q not in bins:
            bins[q]=[0]*(1<<h)
        bins[q][t]=value
    all_children=[{} for _ in translations]
    for q,small in bins.items():
        for out,value in enumerate(fwht(small)):
            if value:
                all_children[out][q]=value
    if any(sum(v*v for v in f.values())!=case['denom']**2 for f in all_children):
        raise AssertionError((case['name'],rows,'real child Parseval'))
    return all_children


def child_pairs(case,child):
    return [(len(f),affine_rank(f,case['n'])) for f in child]


def pair_score(p):
    s,R=p
    return Fraction((R+1)**2,s)


def bits_mask(n,points):
    bits=bytearray(((1<<n)+7)//8)
    for u in points:
        bits[u>>3]|=1<<(u&7)
    return int.from_bytes(bits,'little')


@lru_cache(maxsize=None)
def universe(n,pivot):
    return bits_mask(n,(u for u in range(1,1<<n) if not u&pivot))


def summarize_type(hist,witnesses,a,b,m,pairs,rows,count=1):
    if not count:
        return
    key=(a,b,m,tuple(sorted(pairs)))
    hist[key]+=count
    if rows and key not in witnesses:
        witnesses[key]=tuple(rows)


def one_diff_planes(case,q,Ds,children,hist,witnesses):
    n,ad=case['n'],case['address_dim']
    ma=(1<<ad)-1
    pivot=1<<(q.bit_length()-1)
    eq=lambda w:min(w,w^q)
    forbidden=bits_mask(n,(eq(w) for w in Ds))
    eligible=universe(n,pivot)&~forbidden
    rowspace_A=[eq(1<<i) for i in range(ad)]
    rowspace_B=[eq(1<<i) for i in range(ad,n)]
    maskA=bits_mask(n,span(basis(rowspace_A,n)))
    maskB=bits_mask(n,span(basis(rowspace_B,n)))
    subspaces=[bits_mask(n,span(affine_directions(F,n))) for F in children]
    s_and_r=child_pairs(case,children)
    for flags in range(16):
        selected=eligible
        for i,mask in enumerate([maskA,maskB,*subspaces]):
            selected&=mask if flags&(1<<i) else ~mask
        count=bit_count(selected)
        if count:
            a=(q<=ma)+(flags&1!=0)
            b=(q&ma==0)+(flags&2!=0)
            pairs=[(s,R-int(flags&(4<<i)!=0)) for i,(s,R) in enumerate(s_and_r)]*2
            w=(selected&-selected).bit_length()-1
            summarize_type(hist,witnesses,a,b,1,pairs,[q,w],count)


POOL_CASES={}


def group_work(task):
    name,first,last,depth=task
    start=time.process_time()
    case,ordered,Ds=POOL_CASES[name]
    n,ad=case['n'],case['address_dim']
    h1=Counter(); h2=Counter(); w1={}; w2={}
    for i in range(first,last):
        q=ordered[i]
        one=restrictions(case,(q,))
        pairs=child_pairs(case,one)
        a,b=queries_type((q,),ad,n)
        summarize_type(h1,w1,a,b,1,pairs,(q,))
        if depth<=1:
            continue
        one_diff_planes(case,q,Ds,one,h2,w2)
        for r in ordered[i+1:]:
            s=q^r
            triple=s in Ds
            if triple and s<r:
                continue
            H=(q,r)
            a,b=queries_type(H,ad,n)
            results=restrictions(case,H)
            summarize_type(h2,w2,a,b,2+triple,child_pairs(case,results),H)
    return name,h1,h2,w1,w2,time.process_time()-start


def finalize(case,hist,witnesses,h):
    n,ad=case['n'],case['address_dim']
    source=[(1<<i) for i in range(n)]
    groups=Counter()
    for (a,b,m,p),count in hist.items():
        groups[a,b]+=count
    for a in range(h+1):
        for b in range(h+1-a):
            count=type_count(n,ad,h,a,b)-groups[a,b]
            if count<0:
                raise AssertionError(('negative no-collision multiplicity',case['name'],h,a,b))
            fixed=[(len(case['coeff']),case['R']-h)]*(1<<h)
            summarize_type(hist,witnesses,a,b,0,fixed,(),count)
    if sum(hist.values())!=gf_binomial(n,h):
        raise AssertionError(('incomplete query counts',case['name'],h))
    parent=pair_score((len(case['coeff']),case['R']))
    by_type=defaultdict(Counter)
    best={}
    overall=None
    for (a,b,m,pairs),c in sorted(hist.items()):
        value=max(map(pair_score,pairs))
        local_type=(a,b,m)
        arow=by_type[local_type]
        arow['n_kernels']+=c
        if value>=parent:
            arow['nondecreasing_kernels']+=c
        if value>parent:
            arow['increasing_kernels']+=c
        test=dict(p2=[value.numerator,value.denominator],child_pairs=pairs,
                  H_basis=witnesses.get((a,b,m,pairs)))
        if local_type not in best or value>Fraction(*best[local_type]['p2']):
            best[local_type]=test
        if overall is None or value>Fraction(*overall['p2']):
            overall=dict(a=a,b=b,active_differences=m,**test)
    return dict(depth=h,total=len_hist(hist),overall_best=overall,
                n_recovering=sum(v.get('nondecreasing_kernels',0) for v in by_type.values()),
                n_improving=sum(v.get('increasing_kernels',0) for v in by_type.values()),
                types=[dict(address_intersection=a,data_intersection=b,active_differences=m,
                            **v,best=best[a,b,m]) for (a,b,m),v in sorted(by_type.items())],
                exact_classes=[dict(address_intersection=a,data_intersection=b,active_differences=m,child_pairs=p,
                                    n_kernels=c) for (a,b,m,p),c in sorted(hist.items())])


def len_hist(hist):
    return sum(hist.values())


def cuts_for(case):
    ordered=directions(case)
    Ds=set(ordered)
    cuts={1:[(q,) for q in ordered],2:[]}
    for i,q in enumerate(ordered):
        for r in ordered[i+1:]:
            s=q^r
            if s in Ds and s<r:
                continue
            cuts[2].append((q,r))
    return cuts


def run_cases(names,workers,depth):
    start=time.perf_counter()
    tasks=[]
    for name in names:
        case=build_fixture(name)
        verify_fixture(case)
        Ds=directions(case)
        POOL_CASES[name]=(case,Ds,set(Ds))
        chunk=max(1,(len(Ds)+12*workers-1)//(12*workers))
        tasks.extend((name,i,min(i+chunk,len(Ds)),depth) for i in range(0,len(Ds),chunk))
    # Long early pair batches first. Alternate fixtures by assigning each the same ranking metric.
    tasks.sort(key=lambda row:len(POOL_CASES[row[0]][1])-row[1],reverse=True)
    values={name:(Counter(),Counter(),{}, {}) for name in names}
    cpu=0
    with mp.get_context('fork').Pool(workers) as pool:
        for name,h1,h2,w1,w2,cpu_seconds in pool.imap_unordered(group_work,tasks,chunksize=1):
            a,b,x,y=values[name]
            a.update(h1); b.update(h2)
            for target,source in ((x,w1),(y,w2)):
                for key,w in source.items():
                    if key not in target or w<target[key]:
                        target[key]=w
            cpu+=cpu_seconds
    records=[]
    for name in names:
        case,ds,_=POOL_CASES[name]
        h1,h2,w1,w2=values[name]
        levels=[finalize(case,h1,w1,1)]
        if depth==2:
            levels.append(finalize(case,h2,w2,2))
        p2=pair_score((len(case['coeff']),case['R']))
        records.append(dict(**case,n_difference_directions=len(ds),P_squared=[p2.numerator,p2.denominator],levels=levels))
    elapsed=time.perf_counter()-start
    return dict(cases=records,wall_seconds=elapsed,sum_worker_cpu_seconds=cpu,
                mean_active_worker_cores=cpu/elapsed,workers=workers,tasks=len(tasks))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=os.cpu_count())
    parser.add_argument('--depth',type=int,choices=(1,2),default=2)
    parser.add_argument('case',nargs='+')
    args=parser.parse_args()
    output=run_cases(args.case,args.workers,args.depth)
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),**output),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    main()
