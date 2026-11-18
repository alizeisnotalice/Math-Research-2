#!/usr/bin/env python3
"""High-transverse equality lifts from complete decision-tree seeds."""
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import platform
import random
import time

import phase_lifts as phase
import transverse as tra
import verify_witnesses as verify


def leaf_paths(N):
    h=N.bit_length()-1
    paths=[]
    for i in range(N):
        node=0; anc=[]; assignment=0
        for position in range(h-1,-1,-1):
            bit=(i>>position)&1
            anc.append(1<<node)
            assignment|=bit<<node
            node=2*node+1+bit
        paths.append((anc,assignment))
    return paths


def seed_truth(N,n,x):
    node=0
    while node<N-1:
        node=2*node+1+((x>>node)&1)
    i=node-(N-1)
    if i==0:
        return 1
    return 1-2*((x>>(N-2+i))&1)


def seed_case(N):
    a=N-1; n=2*N-2; R=n+1
    coeff={}
    for i,(paths,assignment) in enumerate(leaf_paths(N)):
        g=0 if i==0 else 1<<(a+i-1)
        for alpha in tra.span(paths):
            coeff[alpha^g]=1-2*phase.parity(alpha&assignment)
    if len(coeff)!=N*N or phase.rank((1<<n)|r for r in coeff)!=R:
        raise AssertionError('seed support and rank checks')
    if n<=15:
        W=phase.fwht(seed_truth(N,n,x) for x in range(1<<n))
        if any(v*N!=coeff.get(u,0)*(1<<n) for u,v in enumerate(W)):
            raise AssertionError('tree seed coefficients differ from actual Fourier transform')
    return dict(name='tree-'+str(N),n=n,R=R,E=[(1<<n)|u for u in sorted(coeff)],masses=[1]*len(coeff),
                denominator=N,coefficients=coeff,L=(2,1),parent=dict(n=n,denom=N,coeff=coeff),truth_source=None)


def matrix_lift(seed,N):
    a=2*N; n=seed['n']; real_input_addr=N-1
    M=[tra.zero(a) for _ in range(n+1)]
    for j in range(N//2):
        node=N//2-1+j
        M[node]=tra.wedge(1<<(2*j),1<<(2*j+1),a)
    fixed=tra.wedge(1,1<<N,a)
    M[n]=fixed
    for j in range(1,N):
        M[real_input_addr+j-1]=tra.sum_matrices([fixed,tra.wedge(1<<j,1<<(N+j),a)],a)
    return dict(seed=seed['name'],kind='tree-high-transverse-equality',a=a,t=a,U=1,generators=M,generator_coords=seed['E'])


def verify_circuit(coeff):
    conv=Counter()
    for (u,c),(v,d) in combinations(coeff.items(),2):
        conv[u^v]+=c*d
    if any(conv.values()):
        raise AssertionError('sparse true Boolean convolution failed')


def verify_tree(N):
    seed=seed_case(N); W=matrix_lift(seed,N)
    a=W['a']; n=seed['n']; addr=N-1
    coeff={}; planes=Counter(); common=None
    for i,(path,asst) in enumerate(leaf_paths(N)):
        g=0 if i==0 else 1<<(addr+i-1)
        for alpha in tra.span(path):
            u=alpha^g
            v=1<<(N+i)
            if alpha&path[-1]:
                v^=1<<(i^1)
            x=1<<i
            polar=tra.sum_matrices([W['generators'][j] for j in tra.points((1<<n)|u)],a)
            if polar!=tra.wedge(x,v,a):
                raise AssertionError('actual rank-two phase differs')
            plane=(0,*sorted((x,v,x^v)))
            planes[plane]+=1
            weight=seed['coefficients'][u]
            for q in plane:
                coeff[u|(q<<n)]=weight*(1-2*(q==(x^v)))
            image=set(plane)-{0}
            common=image if common is None else common&image
    if len(planes)!=2*N or any(v!=N//2 for v in planes.values()):
        raise AssertionError('image plane multiplicity changed')
    if len(coeff)!=4*N*N or phase.rank((1<<(a+n))|u for u in coeff)!=4*N-1 or common:
        raise AssertionError('wrong exact true lift geometry')
    verify_circuit(coeff)
    hull_ghost=(1<<n)|(1<<(addr+0))|(1<<(addr+1))  # combination g0+g1+g2
    H=tra.sum_matrices([W['generators'][j] for j in tra.points(hull_ghost)],a)
    if phase.rank(H)!=6:
        raise AssertionError('missing hull label does not have rank six')
    sample_count=0; dense=None
    if N==4:
        base,W_true,direct=verify.reconstruct(seed,W)
        if any(v*(2*N)!=coeff.get(u,0)*(1<<(n+a)) for u,v in enumerate(W_true)):
            raise AssertionError('exact signed sparse vs actual dense child failed')
        dense=len(W_true)
    else:
        rng=random.Random(850831+N)
        for i in range(512):
            x=rng.getrandbits(n); z=rng.getrandbits(a); y=x|(z<<n)
            Q=[verify.quadratic(M,z) for M in W['generators']]
            shift=sum(b<<i for i,b in enumerate(Q[:n]))
            actual=seed_truth(N,n,x^shift)*(1-2*Q[n])
            estimate=sum(v*(1-2*phase.parity(u&y)) for u,v in coeff.items())
            if actual*(2*N)!=estimate:
                raise AssertionError('actual pointwise translated function differs')
            sample_count+=1
    return dict(N=N,seed_support=N*N,seed_rank=2*N-1,lift_s=4*N*N,lift_rank=4*N-1,t=a,
                seed_l1=N,lift_l1=2*N,L_seed=2,L_lift=2,whole_hull_ghost_rank=6,image_planes=len(planes),
                signed_sparse_convolution_pass=True,actual_dense_signed_coefficients=dense,
                sampled_actual_values=sample_count,plane_atom_old_mass=[1,2]),seed,W


def rare_six():
    seed=next(x for x in tra.compute_states() if x['name']=='rare-2')
    a=6; c,e1,e2,e3,v,w=[1<<j for j in range(a)]
    base=tra.wedge(c,e1,a)
    generators=(tra.wedge(c,e2,a),tra.wedge(c,e3,a),tra.sum_matrices([tra.wedge(v,w,a),base],a),tra.wedge(c,v,a),base)
    spec=dict(seed='rare-2',kind='low-rank-actual-high-rank-hull',a=6,t=6,U=1,
              generators=generators,generator_coords=seed['E'])
    cert=verify.audit(seed,spec,'rare-six')
    return cert,spec


if __name__=='__main__':
    start=time.perf_counter(); group=[]
    for N in (4,8,16):
        r,seed,W=verify_tree(N); group.append(r)
        if N==4:
            export_seed=seed; export=W
    cert,rare=rare_six()
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),wall_seconds=time.perf_counter()-start,
                          complete_trees=group,t6_rare=cert),indent=2,sort_keys=True))
    src=json.loads(Path('signed-witnesses-input.json').read_text())
    src['minima']['rare-six']=rare
    src['minima']['tree-four']=export
    m=sum((seed_truth(4,6,x)==-1)<<x for x in range(64))
    src['additional_seeds']={export_seed['name']:dict(n=6,truth_mask=m)}
    Path('signed-witnesses-input.json').write_text(json.dumps(src,indent=2,sort_keys=True)+'\n')
