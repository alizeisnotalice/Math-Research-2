"""Exact reassembly of compressed Walsh spectra for the R74 gluing lemma."""
import collections
import itertools
import json
import math
import platform
import random


def dot(a, b):
    return (a & b).bit_count() % 2


def rank(vectors):
    basis = {}
    for a in vectors:
        while a and a.bit_length() in basis:
            a ^= basis[a.bit_length()]
        if a:
            basis[a.bit_length()] = a
    return len(basis)


def walsh(f):
    return {a:v for a in range(len(f))
            if (v:=sum(s*(-1)**dot(a,x) for x,s in enumerate(f)))}


def compress(spectrum, c, j, matrix, shift):
    out = collections.Counter()
    for a,v in spectrum.items():
        k = a & ((1<<c)-1)
        image = sum(dot(k,row)<<i for i,row in enumerate(matrix))
        out[(a>>c)<<j | image] += v*(-1)**dot(k,shift)
    return {a:v for a,v in out.items() if v}


def trial(children, j, rng):
    M=1<<j
    g=len(children); adr=g.bit_length()-1
    best=[]; Q=R=0; trials=0
    for f,c in children:
        si=walsh(f)
        q=sum(1 for a in si if a>>c)
        r=rank(a>>c for a in si)
        Q+=q; R+=r
        best_score=(-1,-1,-1); average_num=0; count=0
        for rows in itertools.product(range(1<<c),repeat=j):
            for shift in range(1<<c):
                compressed=compress(si,c,j,rows,shift)
                actual=rank(a>>j for a in compressed)
                count+=1; average_num+=actual
                score=actual,sum(a<M for a in compressed),len(compressed)
                if score>best_score:
                    best_score=score
                    winner=compressed,rows,shift
        assert 2*M*(r*count-average_num)<=(q-r)*count
        best.append(winner); trials+=count
    whole=collections.Counter(); offset=j+adr
    offsets=[]
    for i,((f,c),(compressed,rows,shift)) in enumerate(zip(children,best)):
        offsets.append(offset)
        for a,v in compressed.items():
            common=a&(M-1); exterior=a>>j
            data=(exterior<<offset)|(common<<adr)
            for z in range(g):
                whole[data|z]+=v*(-1)**dot(z,i)
        offset+=(len(f).bit_length()-1-c)
    whole={a:v for a,v in whole.items() if v}
    assert len(whole)<=g*(Q+M)
    assert 2*M*(R-rank(whole))<=Q-R
    for t in range(1<<offset) if offset<=7 else (rng.randrange(1<<offset) for _ in range(8)):
        i=t&(g-1); z=(t>>adr)&(M-1); f,c=children[i]
        rows,shift=best[i][1:]
        raw=shift
        for a,row in enumerate(rows):
            if z&(1<<a): raw^=row
        raw|=((t>>offsets[i])&((len(f)>>c)-1))<<c
        assert sum(v*(-1)**dot(a,t) for a,v in whole.items())==len(f)*g*f[raw]
    return dict(maps_with_shifts=trials,s=len(whole),d=rank(whole),
                Q=Q,R=R,common_grid_bound=g*(Q+M))


def main():
    rng=random.Random(74073)
    cases=140; maps=0
    for _ in range(cases):
        b=rng.choice((2,4,8))
        prepared=[([rng.choice((-1,1)) for _ in range(8)],rng.randint(0,3))
                  for _ in range(b)]
        sample=rng.sample(prepared,b//2)
        maps+=trial(sample,rng.randint(0,2),rng)['maps_with_shifts']
    f=[(-1)**((x>>1)&1) if (x&1) else 1 for x in range(4)]
    simple=trial([(f,1)]*4,1,rng)
    print(json.dumps(dict(python=platform.python_version(),seed=74073,
        random_regluings=cases,linear_maps_and_shifts=maps,
        common_band_example=simple,
        checks="expected actual exterior rank, exact sparsity bound, assembled Boolean evaluations"),
        indent=2))


if __name__=='__main__':
    main()
