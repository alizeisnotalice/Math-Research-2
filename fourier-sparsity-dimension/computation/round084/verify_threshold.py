"""Exact two-layer spectral folding and threshold-path replay."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from random import Random
import json
import platform


def rank(xs):
    b={}
    for v in xs:
        while v and v.bit_length() in b:
            v ^= b[v.bit_length()]
        if v:
            b[v.bit_length()]=v
    return len(b)


def rows(H,A):
    b={}
    for v,s in zip(H,A):
        while v and v.bit_length() in b:
            t,a=b[v.bit_length()]
            v ^= t; s ^= a
        if v:
            b[v.bit_length()]=v,s
    return sorted(b.values(),reverse=True)


def folding(S,H,A):
    B=rows(H,A)
    out=defaultdict(int)
    for x,c in S.items():
        sign=0
        for v,b in B:
            if x & (1<<(v.bit_length()-1)):
                x ^= v; sign ^= b
        out[x]+=c*(1-2*sign)
    return {x:c for x,c in out.items() if c}


def orth(edge):
    a,b=edge
    bs=[]
    for k in range(1,16):
        if ((a^b)&k).bit_count()%2==0 and rank(bs+[k])>len(bs):
            bs.append(k)
    assert len(bs)==3
    return bs,[(a&k).bit_count()%2 for k in bs]


def run():
    rng=Random(8407)
    edges=list(range(16))
    rng.shuffle(edges)
    edges=list(zip(edges[::2],edges[1::2]))
    bands=[]
    for a,b in edges:
        bands.append({i:1-2*((i&a).bit_count()%2) for i in range(16)
                      if not((i&(a^b)).bit_count()%2)})
    S={}; inner_addr=[]; guards=[]; inner_ranges=[]; offset=4
    for i,outer in enumerate(bands):
        inner_addr.append(offset)
        start=offset
        offset+=4
        gs=[]
        for j,inner in enumerate(bands):
            # Each bent seed has an independent guard and two genuine variables.
            g=1<<(offset+2); gs.append(g)
            for v,c in outer.items():
                for w,a in inner.items():
                    for z in range(4):
                        S[v|(w<<inner_addr[-1])|g|(z<<offset)]=c*a*(1 if z!=3 else -1)
            offset+=3
        guards.append(gs)
        inner_ranges.append(((1<<offset)-1)^((1<<start)-1))
    marker=1<<offset
    assert rank(marker|v for v in S)==228
    assert len(S)==16384
    parent=lambda sp:dict(s=len(sp),D=rank(marker|x for x in sp)+1,
                           P_squared=str(Fraction((rank(marker|x for x in sp)+1)**2,len(sp))))
    seedR=3; inners=256; innerr=28
    checked=0
    def calc(H,A):
        nonlocal checked
        sp=folding(S,H,A)
        assert sum(c*c for c in sp.values())==128**2
        constraints=[(0,0)]
        for x,b in zip(H,A):
            constraints += [(v^x,s^b) for v,s in constraints]
        pure=[(x,b) for x,b in constraints if x<16]
        addr=[u for u in range(16) if all((u&x).bit_count()%2==b for x,b in pure)]
        f=sum(a in addr and b in addr for a,b in edges)
        t=sum((a in addr)^(b in addr) for a,b in edges)
        D=rank(marker|x for x in sp)+1
        assert D<=(f+t)*innerr+len(addr).bit_length()
        assert len(sp)>=inners*((f+t)**2-2*len(H)*len(addr))
        checked+=1
        return sp
    H,A=orth(edges[0])
    outer=calc(H,A)
    assert parent(outer)['P_squared']==str(Fraction(innerr+1,16)**2)
    h2,a2=orth(edges[0])
    combined=calc(H+[x<<inner_addr[0] for x in h2],A+a2)
    assert parent(combined)['P_squared']=='4'
    pure_recover=None
    for p,q in combinations(range(1,16),2):
        for v in product((0,1),repeat=2):
            cell=[u for u in range(16) if ((u&p).bit_count()%2,(u&q).bit_count()%2)==v]
            if not any(a in cell and b in cell for a,b in edges):
                sp=calc([p,q],v)
                pure_recover=parent(sp)
                break
        if pure_recover: break
    assert Fraction(parent(S)['P_squared'])<Fraction(pure_recover['P_squared'])<4
    effective=list(S)
    for _ in range(35):
        dim=rng.choice((1,2,3,4))
        H=[]
        while len(H)<dim:
            i=rng.randrange(8); j=rng.randrange(8)
            v=rng.choice([rng.randrange(1,16),rng.randrange(1,16)<<inner_addr[i],
                          rng.choice(effective)^rng.choice(effective),
                          guards[i][rng.randrange(8)]^guards[j][rng.randrange(8)]])
            if rank(H+[v])>len(H): H.append(v)
        calc(H,[rng.randrange(2) for _ in H])
    return dict(seed=8407,matching=edges,query_answer_pairs=checked,
                root=parent(S),two_outer_queries=pure_recover,
                outer_edge_three_queries=parent(outer),both_edges_six_queries=parent(combined))


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,example=run()),indent=2))
