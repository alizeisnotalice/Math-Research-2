#!/usr/bin/env python3
"""Exact feasible simultaneous rank-preserving extensions of Petersen pencils."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import platform
import random
import time

import petersen as pet


def dot(x,y):return pet.parity(x&y)


def edge_cycles():
    return pet.kernel([sum(((u==x or v==x)<<i) for i,(u,v) in enumerate(pet.EDGES)) for x in range(10)],15)


def edge_values(T,a):
    p=[0]+[(T>>(i*a))&((1<<a)-1) for i in range(9)]
    return [p[u]^p[v] for u,v in pet.EDGES]


def reduce(x,P):
    while x:
        h=x.bit_length()-1
        if h not in P:return x
        x^=P[h]
    return 0


def preimage(B,q):
    rows={}
    for i,u in enumerate(B):
        t=1<<i
        while u:
            h=u.bit_length()-1
            if h not in rows:
                rows[h]=(u,t); break
            v,w=rows[h];u^=v;t^=w
    u=q;t=0
    while u:
        v,w=rows[u.bit_length()-1];u^=v;t^=w
    return t


def potentials(M):
    a=len(M[0]);nodes=[None]*10;nodes[0]=(0,)*a
    while any(v is None for v in nodes):
        for Q,(u,v) in zip(M,pet.EDGES):
            if nodes[u] is not None and nodes[v] is None:nodes[v]=pet.xor(nodes[u],Q)
            if nodes[v] is not None and nodes[u] is None:nodes[u]=pet.xor(nodes[v],Q)
    if any(pet.xor(nodes[u],nodes[v])!=tuple(Q) for Q,(u,v) in zip(M,pet.EDGES)):
        raise AssertionError('edge pencil failed cycle consistency')
    return nodes


def quotient_basis(M):
    a=len(M[0]); rows=[]
    for Q,(u,v) in zip(M,pet.EDGES):
        for x in pet.kernel(Q,a):
            rows.append((x<<((u-1)*a) if u else 0)^(x<<((v-1)*a) if v else 0))
    nodes=potentials(M)
    local=[sum(nodes[i][j]<<((i-1)*a) for i in range(1,10)) for j in range(a)]
    if pet.rank(local)!=a:
        raise AssertionError('local rank does not span address')
    selected=[];K=pet.kernel(rows,9*a);current=a
    for v in K:
        if pet.rank([*local,*selected,v])>current:
            selected.append(v);current+=1
    if len(selected)!=len(K)-a:
        raise AssertionError('nonlocal quotient rank wrong')
    return selected


def six_forms(M,nonlocal_basis):
    a=len(M[0]);n=len(nonlocal_basis)
    init=[edge_values(x,a) for x in nonlocal_basis]
    lifted=[[preimage(Q,v) for Q,v in zip(M,T)] for T in init]
    cycles=edge_cycles()
    forms=[[0]*n for _ in range(n)]
    for i,T in enumerate(init):
        for j,S in enumerate(lifted):
            edge_cov=sum(dot(v,w)<<q for q,(v,w) in enumerate(zip(T,S)))
            forms[i][j]=sum(dot(edge_cov,C)<<c for c,C in enumerate(cycles))
    if any(forms[i][i] or forms[i][j]!=forms[j][i] for i in range(n) for j in range(n)):
        raise AssertionError('quotient cycle forms are not alternating')
    return forms


def value(forms,left,right):
    x=0
    for i in pet.positions(left):
        for j in pet.positions(right):
            x^=forms[i][j]
    return x


def largest_isotropic(forms,cap=7):
    n=len(forms)
    if not any(any(row) for row in forms):return [1<<i for i in range(n)],True,1
    if n>cap:
        found=[]
        for v in range(1,1<<min(n,12)):
            if pet.rank([*found,v])>len(found) and all(value(forms,u,v)==0 for u in found):found.append(v)
        return found,False,0
    levels={1:()};maxbasis=(); attempted=0
    while levels:
        nxt={}
        for mask,bs in levels.items():
            for q in range(1,1<<n):
                if mask>>q&1 or any(value(forms,u,q) for u in bs):continue
                attempted+=1
                outside=sum(1<<(x^q) for x in pet.positions(mask))
                newer=mask|outside
                nxt.setdefault(newer,(*bs,q))
        maxbasis=next(iter(levels.values()))
        levels=nxt
    return maxbasis,True,attempted


def explicit_extension(M,basis,bs):
    a=len(M[0]); q=len(bs)
    lifted=[]
    for b in bs:
        T=0
        for j in pet.positions(b):T^=basis[j]
        lifted.append(edge_values(T,a))
    output=[]
    for e,Q in enumerate(M):
        C=[z[e] for z in lifted]
        pre=[preimage(Q,col) for col in C]
        old=[y|sum((v>>i&1)<<(a+j) for j,v in enumerate(C)) for i,y in enumerate(Q)]
        new=[v|sum(dot(v,h)<<(a+j) for j,h in enumerate(pre)) for v in C]
        newQ=old+new
        if pet.rank(newQ)!=pet.rank(Q):raise AssertionError('extension changes actual polar ranks')
        output.append(newQ)
    potentials(output)
    if pet.rank(row for Q in output for row in Q)!=a+q:
        raise AssertionError('extension failed exactly predicted effective new rank')
    return output


def certify(tag,M):
    a=len(M[0]); basis=quotient_basis(M); forms=six_forms(M,basis)
    top,complete,computations=largest_isotropic(forms)
    more=explicit_extension(M,basis,top)
    parent=pet.estimate(M,a,True); child=pet.estimate(more,len(more[0]),True)
    return dict(tag=tag,original=parent,old_matrices=M,nonlocal_dimension=len(basis),quotient_six_forms=forms,
                maximum_feasible_axes=len(top) if complete else None,certified_axes=len(top),coordinates_of_axes=top,
                is_exhaustive=complete,isotropic_partial_extensions=computations,
                updated=child,new_matrices=more)


def scan():
    records=[]
    for f in ('petersen-mini.json','petersen-macbook.json'):
        d=json.load(open(f)); M=d['best_nonlocal']['old_matrices'];records.append(certify(f+' best nonlocal',M))
        records.append(certify(f+' common-star',d['best']['all']['matrices']))
        records.append(certify(f+' centered-mixed',d['best']['mixed_centered']['matrices']))
    rng=random.Random(860831)
    hist=Counter(); mins={}; selected=[]
    for i in range(6144):
        kind=pet.KINDS[i%len(pet.KINDS)]
        raw,center=pet.choose(kind,rng)
        case=pet.normalize(raw,8,rng)
        if case is None:continue
        M,_=case; a=len(M[0]); q=pet.compatible(M,False)['nonlocal_dimension']
        hist['pencils']+=1;hist['nonlocal']+=bool(q)
        if not q:continue
        B=Fraction(*pet.estimate(M,a,False)['margin_X0'])
        prob=B-Fraction(3*q,2)
        hist['upper_bound_reaches_negative']+=prob<0
        if prob<=1 or q>=2 and not any(z[0]==q for z in selected):
            selected.append((q,kind,M,B,prob))
    for q,kind,M,B,b in selected:
        if q not in mins or b<mins[q][0]:mins[q]=b,kind,M
    for q,(_,kind,M) in sorted(mins.items()):
        r=certify(kind+' compatible dimension '+str(q),M); records.append(r)
        hist['exact_isotropic_cases']+=1
        hist['negative_actual_cut']+=Fraction(*r['updated']['margin_X0'])<0
    return dict(seed=860831,search=hist,certificates=records)


if __name__=='__main__':
    st=time.perf_counter(); r=scan()
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),
                          wall_seconds=time.perf_counter()-st,**r),indent=2,sort_keys=True))
