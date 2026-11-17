#!/usr/bin/env python3
"""Exact k=2,n=4 enumeration and bounded random affine-flat searches.

A codimension-k affine flat is stored by a canonical basis of its normal
space and the k right-hand-side bits on that basis.  All intersection tests
are exact bitset operations on the ambient finite cube.
"""
import argparse, hashlib, json, random, sys, time
from dataclasses import dataclass
from pathlib import Path

def canonical(vectors, n):
    piv=[0]*n
    for x in vectors:
        for b in range(n-1,-1,-1):
            if x>>b&1:
                if piv[b]: x ^= piv[b]
                else:
                    # Reduce the new leading row by all existing lower pivots
                    # before installing it; this makes the returned RREF basis
                    # canonical for a subspace, rather than only echelon.
                    for c in range(b-1,-1,-1):
                        if x>>c&1 and piv[c]: x ^= piv[c]
                    piv[b]=x
                    for c in range(n):
                        if c!=b and piv[c]>>b&1:piv[c]^=x
                    break
    return tuple(x for x in piv if x)
def rank(vectors,n):return len(canonical(vectors,n))
def normal_spaces(n,k):
    spaces=set()
    # Grow spans rather than ordered generator tuples, then set-deduplicate.
    levels={()}
    for _ in range(k):
        nxt=set()
        for s in levels:
            for x in range(1,1<<n):
                if rank(s+(x,),n)==len(s)+1:nxt.add(canonical(s+(x,),n))
        levels=nxt
    return sorted(levels)
def chars(n):
    full=(1<<(1<<n))-1
    return [sum(((x&a).bit_count()&1)<<x for x in range(1<<n)) for a in range(1<<n)],full
def flat_mask(basis, offset, char, full):
    mask=full
    for j,a in enumerate(basis):mask &= char[a] if offset>>j&1 else full^char[a]
    return mask
@dataclass(frozen=True)
class Flat: basis:tuple; offset:int; mask:int
def exhaustive_k2_n4():
    n,k=4,2; char,full=chars(n); spaces=normal_spaces(n,k)
    flats=[Flat(s,o,flat_mask(s,o,char,full)) for s in spaces for o in range(1<<k)]
    counts=[0]*5; maxrank=0; witness=None
    def dfs(start,chosen,used,normals):
        nonlocal maxrank,witness
        if chosen:
            r=rank(normals,n);counts[len(chosen)]+=1
            if r>maxrank:maxrank=r;witness=chosen[:]
            if r>3: raise AssertionError('counterexample')
        if len(chosen)==4:return
        for i in range(start,len(flats)):
            f=flats[i]
            if not used&f.mask:dfs(i+1,chosen+[f],used|f.mask,normals+list(f.basis))
    dfs(0,[],0,[])
    return {'n':n,'k':k,'normal_spaces':len(spaces),'flats':len(flats),'family_counts_by_m':counts[1:],
            'max_normal_rank':maxrank,'maxrank_witness':encode(witness),'counterexample_found':False}
def encode(flats):return [{'normal_basis':list(f.basis),'offset_bits':f.offset,'point_mask_hex':hex(f.mask)} for f in (flats or [])]
def random_basis(n,k,rng):
    out=[]
    while len(out)<k:
        x=rng.randrange(1,1<<n)
        if rank(out+[x],n)>len(out):out.append(x)
    return canonical(out,n)
def bounded_greedy(n,k,seed,trials,target):
    rng=random.Random(seed);char,full=chars(n); best=[]; counter=None
    for t in range(trials):
        chosen=[]; used=0; normals=[]
        for step in range(1<<k):
            options=[]
            # Sample a bounded choice pool; all candidates are exact flats.
            for _ in range(48):
                b=random_basis(n,k,rng)
                for off in range(1<<k):
                    m=flat_mask(b,off,char,full)
                    if not used&m: options.append(Flat(b,off,m))
            if not options:break
            rng.shuffle(options)
            options.sort(key=lambda f:rank(normals+list(f.basis),n),reverse=True)
            f=options[0];chosen.append(f);used|=f.mask;normals+=list(f.basis)
            r=rank(normals,n)
            if r>=target:
                counter=chosen[:];break
        if len(chosen)>len(best):best=chosen[:]
        if counter:return {'n':n,'k':k,'seed':seed,'trials':trials,'target_rank':target,'counterexample_found':True,
                           'normal_rank':rank([a for f in counter for a in f.basis],n),'flats':encode(counter),'trial_found':t}
    return {'n':n,'k':k,'seed':seed,'trials':trials,'target_rank':target,'counterexample_found':False,
            'best_m':len(best),'best_normal_rank':rank([a for f in best for a in f.basis],n),'best_flats':encode(best)}
def punctured_partition_lift(n,k,seed,trials,target):
    """Delete three cells of a standard full partition, then seek free normals."""
    rng=random.Random(seed);char,full=chars(n); base=tuple(1<<j for j in range(k)); best=[]; counter=None
    for t in range(trials):
        kept=rng.sample(range(1<<k),5)
        chosen=[Flat(base,o,flat_mask(base,o,char,full)) for o in kept];used=0
        for f in chosen:used|=f.mask
        normals=[a for f in chosen for a in f.basis]
        while len(chosen)<1<<k:
            options=[]
            for _ in range(96):
                b=random_basis(n,k,rng)
                for off in range(1<<k):
                    m=flat_mask(b,off,char,full)
                    if not used&m:options.append(Flat(b,off,m))
            if not options:break
            options.sort(key=lambda f:rank(normals+list(f.basis),n),reverse=True)
            f=options[0];chosen.append(f);used|=f.mask;normals+=list(f.basis)
            if rank(normals,n)>=target:counter=chosen[:];break
        if rank(normals,n)>rank([a for f in best for a in f.basis],n):best=chosen[:]
        if counter:break
    chosen=counter or best
    return {'n':n,'k':k,'seed':seed,'trials':trials,'method':'punctured_standard_partition_then_free_normal_greedy',
            'target_rank':target,'counterexample_found':counter is not None,'best_m':len(chosen),
            'best_normal_rank':rank([a for f in chosen for a in f.basis],n),'best_flats':encode(chosen)}
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--trials',type=int,default=20000);a=ap.parse_args()
beg=time.monotonic(); exact=exhaustive_k2_n4(); searches=[]
for n in (8,9):
    for seed in (2026111801,2026111802):searches.append(bounded_greedy(n,3,seed,a.trials,8))
searches.append(punctured_partition_lift(8,3,2026111803,max(1,a.trials//4),8))
result={'status':'passed','scope':'exact k2n4 plus bounded random k3n8n9 counterexample search; negative random results are not proof',
        'exact_k2_n4':exact,'bounded_k3_searches':searches,'elapsed_seconds':time.monotonic()-beg,
        'python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','k2_maxrank':exact['max_normal_rank'],'k3_counterexamples':sum(x['counterexample_found'] for x in searches),'elapsed_seconds':result['elapsed_seconds']}))
