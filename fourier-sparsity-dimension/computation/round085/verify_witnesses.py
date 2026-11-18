#!/usr/bin/env python3
"""Actual Boolean transform and complete data cut audit for the archived pencils."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import struct
import time

import phase_lifts as phase
import transverse as tra


def quadratic(M,x):
    return tra.xor_all(phase.parity(M[i]&x&~((1<<(i+1))-1)) for i in tra.points(x))


def all_homogeneous_generator_coords(seed,coords):
    points=[]; gamma=[]
    for e,c in zip(seed['E'],coords):
        if phase.rank([*points,e])>len(points):
            points.append(e);gamma.append(c)
    E=phase.quotient_coords(points)
    return [tra.xor_all(gamma[j] for j in tra.points(E[1<<k])) for k in range(seed['R'])]


def reconstruct(seed,w):
    a=w['a']; n=seed['n']; R=seed['R']
    generators=w['generators']; coords=w['generator_coords']
    unit_coords=all_homogeneous_generator_coords(seed,coords)
    affine_M=[tra.sum_matrices([generators[j] for j in tra.points(c)],a) for c in unit_coords]
    signs=phase.truth_from_sparse(seed['parent'])
    W_base=phase.fwht(signs)
    actual=[]
    for z in range(1<<a):
        shift=sum(quadratic(M,z)<<j for j,M in enumerate(affine_M[:n]))
        flip=quadratic(affine_M[n],z)
        actual.extend((1-2*flip)*signs[x^shift] for x in range(1<<n))
    W=phase.fwht(actual)
    support={q:c for q,c in enumerate(W) if c}
    if phase.rank((1<<(n+a))|q for q in support)!=R+a:
        raise AssertionError('actual true parent homogeneous rank')
    for beta,e,tag in zip(sorted(seed['coefficients']),seed['E'],coords):
        M=tra.sum_matrices([generators[j] for j in tra.points(tag)],a)
        P=phase.fwht(1-2*quadratic(M,z) for z in range(1<<a))
        if sum(c!=0 for c in P)!=(1<<phase.rank(M)):
            raise AssertionError('phase support disagrees with its polar rank')
        for alpha,c in enumerate(P):
            if W_base[beta]*c!=W[beta|(alpha<<n)]:
                raise AssertionError('actual true Walsh differs from phase multiplication')
    return actual,W,support


def all_data_cuts(seed,w):
    dims=[len(cuts) for cuts in tra.subspaces(seed['R'])]
    by=Counter(); minimum=None; low=None
    for level in tra.subspaces(seed['R']):
        for U in level:
            row=tra.costs(seed,w['generators'],w['generator_coords'],w['a'],w['t'],U)
            ratio=Fraction(*row['exact_ratio'])
            greedy=Fraction(*row['greedy_ratio'])
            by['weighted_insufficient']+=ratio<1
            by['weighted_ties']+=ratio==1
            by['greedy_insufficient']+=greedy<1
            if minimum is None or ratio<minimum:
                minimum=ratio; low=row
    return dict(data_subspaces=sum(dims),counts=by,worst_weighted_ratio=[minimum.numerator,minimum.denominator],
                worst_U=low['U'],worst_u=low['u'])


def audit(seed,w,label):
    actual,W,support=reconstruct(seed,w)
    a=w['a']; n=seed['n']; size=1<<(a+n)
    mass=Fraction(sum(abs(c) for c in W),size)
    L=Fraction(*seed['L']); Xzero=Fraction(seed['R']+a+1,1)/mass
    ranks=Counter()
    unit=all_homogeneous_generator_coords(seed,w['generator_coords'])
    nonzero_image=None
    for beta in range(1<<n):
        e=(1<<n)|beta
        gamma=tra.xor_all(unit[j] for j in tra.points(e))
        M=tra.sum_matrices([w['generators'][i] for i in tra.points(gamma)],a)
        p=phase.pivots(M); ranks[len(p)]+=1
        if seed['coefficients'].get(beta,0) and p:
            image=set(tra.span(p))-{0}
            nonzero_image=image if nonzero_image is None else nonzero_image&image
    census=all_data_cuts(seed,w)
    if label=='weighted' and a==5 and census['counts']['weighted_insufficient']==0 and w['U']==1 and Xzero==L:
        theorem_data=dict(exact_L=list(seed['L']),uses_R84_theorem_for_t_le4=True)
    else:
        theorem_data=None
    return dict(witness_label=label,record=w,seed_L=seed['L'],seed_truth=seed['truth_source'],
                parent_s=len(support),parent_R=a+seed['R'],parent_l1=[mass.numerator,mass.denominator],
                actual_X0_ratio=[Xzero.numerator,Xzero.denominator],all_data_cuts_for_same_J=census,
                hull_polar_rank_histogram=ranks,has_common_nonzero_image_vector=bool(nonzero_image),
                exact_parent_L_certificate=theorem_data,
                parent_signed_walsh_sha256=hashlib.sha256(struct.pack('<'+'i'*len(W),*W)).hexdigest(),
                signed_nonzero_parent_Walsh=support)


def plane_matrices(a):
    return sorted({tra.wedge(u,v,a) for u,v in combinations(range(1,1<<a),2)})


def enumerate_fixed_family(seed,a):
    # Actual rare-2 support: 0,A,B,A+B,C phases, each occurring twice.
    A=tra.wedge(1,2,a); matrices=plane_matrices(a)
    if len(matrices)!=(((1<<a)-1)*((1<<a)-2)//6):
        raise AssertionError('rank-two atlas is incomplete')
    R=seed['R']; affine_space=(1<<(1<<R))-1
    if set(seed['coefficients'])!={0,1,2,3,4,8,9,10,11,12}:
        raise AssertionError('rare-2 support/coordinate assumption')
    by=Counter(); representatives={}
    for B in matrices:
        for C in matrices:
            if phase.rank([*A,*B,*C]) != a:
                continue
            d=phase.rank(tra.sum_matrices([A,B],a)); by[d]+=1
            representatives.setdefault(d,(A,B,C,tra.zero(a),tra.zero(a)))
    out=[]
    for d,c in sorted(by.items()):
        w=dict(generators=representatives[d],generator_coords=seed['E'],a=a,t=a)
        test=all_data_cuts(seed,w)
        out.append(dict(A_plus_B_rank=d,pencils=c,**test,
                        historical_A_B_C=representatives[d][:3]))
    return dict(a=a,fixed_first_wedge=A,raw_BC_pair_count=len(matrices)**2,fullspan_pencils=sum(by.values()),acc=out,
                covered_pencil_U_pairs=sum(z['pencils']*z['data_subspaces'] for z in out))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('witness_file')
    args=p.parse_args()
    start=time.perf_counter(); seeds={s['name']:s for s in tra.compute_states()}
    doc=json.loads(Path(args.witness_file).read_text())
    evidence=[]
    for label in ('basis','weighted'):
        w=doc['minima'][label]
        evidence.append(audit(seeds[w['seed']],w,label))
    records=[enumerate_fixed_family(seeds['rare-2'],a) for a in (5,6)]
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),wall_seconds=time.perf_counter()-start,
                          evidence=evidence,exhaustive_rank2_generator_family=records),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
