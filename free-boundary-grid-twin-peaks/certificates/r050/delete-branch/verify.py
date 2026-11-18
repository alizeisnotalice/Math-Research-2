from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib
import json
import platform

P = Path(__file__).resolve().parent
CERT = P.parents[1]
NORMALIZERS={3:143808,4:779299174080,5:7063261481362958280000}

def n3_punctured():
    # Independent enumeration by punctured words, using the inverse formula.
    n=3
    V=frozenset(range(n*n))
    adj=[frozenset(u for u in V if abs(u%n-v%n)+abs(u//n-v//n)==1) for v in V]

    @cache
    def words(vertices, prefix):
        if prefix==vertices:
            return 1
        return sum(words(vertices,prefix|{v}) for v in vertices-prefix
                   if adj[v]&prefix)

    def grew(vertices, order, seen):
        if seen==vertices:
            yield order
        else:
            for v in sorted(vertices-seen):
                if not seen or adj[v]&seen:
                    yield from grew(vertices,order+(v,),seen|{v})

    out=Counter()
    number=0
    for b in range(1,8):
        for branch in combinations(V,b):
            B=frozenset(branch)
            roots={q:words(B,frozenset({q})) for q in B}
            if not any(roots.values()):
                continue
            boundary=frozenset.union(*(adj[q] for q in B)) - B
            for sigma in grew(V-B,(),frozenset()):
                s=next(i for i,v in enumerate(sigma,1) if v in boundary)
                if s<2:
                    continue
                m=s+b-2
                ways=comb(m,b)
                assert sum(comb(m-(k-1),b-1) for k in range(2,s+1))==ways
                p=sigma[0]
                for q,A in roots.items():
                    if A:
                        dist=abs(q%n-p%n)+abs(q//n-p//n)
                        assert dist>=2
                        out[s+b,b,dist]+=A*ways
                number+=1
    return out,number

def check(n):
    path=P/('n3.tsv' if n==3 else f'n{n}-b8.tsv')
    original={}
    removed={}
    by_b=defaultdict(list)
    for line in path.read_text().splitlines():
        i,b,d,f,t=map(int,line.split())
        assert 0<=t<=f
        original[i,b,d]=f
        if t:
            removed[i,b,d]=t
        by_b[b].append((d,f,t))
    prev=Counter()
    d= CERT/'r049/low-branch'/f'n{n}-contact.tsv'
    for line in d.read_text().splitlines():
        i,b,h,r,D,c=map(int,line.split())
        prev[i,b,r]+=c
    assert dict(prev)==original
    used=sum(original.values())
    legal=sum(removed.values())
    answer=dict(n=n,cells=len(original),covered_words=used,deletable_words=legal,
                covered_fraction=str(Fraction(used,NORMALIZERS[n])),
                deletable_fraction_in_size_scope=str(Fraction(legal,used)))
    if n==3:
        assert legal==103264
        independent,number=n3_punctured()
        assert dict(independent)==removed
        answer['checked_punctured_pairs']=number
    results=[]
    for b, rows in sorted(by_b.items()):
        a=sum(f for _,f,_ in rows)
        z=sum(t for _,_,t in rows)
        d2a=sum(f for d,f,t in rows if d==2)
        d2z=sum(t for d,f,t in rows if d==2)
        results.append(dict(b=b,all=a,deletable=z,fraction=str(Fraction(z,a)),
                            distance2_fraction=str(Fraction(d2z,d2a))))
    answer['by_b']=results
    answer['data_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    return answer

if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),checks=[check(n) for n in (3,4,5)]),indent=2))
