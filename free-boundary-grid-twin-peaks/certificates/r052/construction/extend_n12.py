import os
"""Disjoint extensions of the archived n12 quadrant numerator."""
from argparse import ArgumentParser
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
import hashlib
import json
from math import comb,factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PROJECT=Path(os.environ.get('A022_RESEARCH_ROOT', Path(__file__).resolve().parents[3]))
R50=PROJECT/'certificates/r050/n11-gluing'
ROOT6=PROJECT/'computation/results/independent-peaks-n6-r004.json'
R51=PROJECT/'certificates/r051/n12/n12-conditional-final.json'
sys.path.insert(0,str(R50))
from merge_n11 import prefix_rank_counts


def read(path):return json.loads(path.read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def dec(x):
    with localcontext() as c:
        c.prec=57
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def triplet(a,b,c):return comb(a+b+c,a)*comb(b+c,b)


def chain_shuffles(M):
    return [[sum(comb(a-i+j-1,j-1)*triplet(M-a,M-j,M) for a in range(i,M+1))
             if i and j else 0 for j in range(M+1)] for i in range(M+1)]


def gate_chain(M):
    W=chain_shuffles(M)
    exact=[[0]*(M+1) for _ in range(M+1)]
    positions=set(range(3*M))
    for first in combinations(range(3*M),M):
        remaining=positions.difference(first)
        for middle in combinations(sorted(remaining),M):
            last_begin=min(remaining.difference(middle))
            for i,s in enumerate(first,1):
                for j,t in enumerate(middle,1):
                    exact[i][j]+=s<middle[0] and t<last_begin
    assert exact==W
    return triplet(M,M,M)


def safe_shuffle(M,ell_l,ell_r):
    """ell_l restricts first L before the ell_l-th R. ell_r is its mirror."""
    A=4*M
    B=2*M
    answer=comb(A,B)
    for ell in (ell_l,ell_r):
        if ell:
            answer-=comb(A-ell,B)
    return answer


def gate_pair(M):
    # Explicit positions of the two pair words, including simultaneous gates.
    B=2*M
    for r in range(B+1):
        for s in range(B+1):
            count=0
            for A in combinations(range(2*B),B):
                C=tuple(sorted(set(range(2*B)).difference(A)))
                count+=(not r or A[0]<C[r-1]) and (not s or C[0]<A[s-1])
            assert count==safe_shuffle(M,r,s)
    return comb(2*B,B)


class Book:
    def __init__(self):
        self.m=6
        self.M=self.m*self.m
        self.paths={'root':ROOT6,'plan':R50/'plan-K7.json','completion':R50/'results-K7/summary.json'}
        p,raw=(read(self.paths[key]) for key in ('plan','completion'))
        assert digest(self.paths['plan'])==raw['plan_sha256']
        self.root=list(map(int,read(ROOT6)['single_peak_root_counts']))
        self.ranks=prefix_rank_counts({'rows':[r for r in p['rows'] if r['w']==r['h']==6]},
                                     {(6,6,int(t['seen'])):int(t['completion']) for t in raw['records'] if t['w']==t['h']==6},
                                     {(6,6):self.root})
        for x in range(self.m):
            for y in range(self.m):
                for c in range(self.m):
                    b=self.rank(x,y,c)
                    assert 0<sum(b)<=self.root[y*self.m+x]
                    assert not any(b[:abs(x-c)+self.m-y])

    def rank(self,x,y,c):
        if x==self.m-1:
            x=0;c=self.m-1-c
        return self.ranks[self.m,self.m,y*self.m+x,c]

    def pair(self):
        m,M=self.m,self.M
        return {(x,y):sum(a*self.root[c]*comb(2*M-j,M) for c in range(m)
                          for j,a in enumerate(self.rank(x,y,c)) if a)
                for x in range(m) for y in range(m)}

    def chain(self):
        m,M=self.m,self.M
        G=chain_shuffles(M)
        # For a fixed connector column c of the first winner, the first
        # loser's root is (c,m-1). Transpose its block to count the ranks
        # of its right-edge connectors. The new loser is off the upper
        # seam, so y=m-1 is excluded from its possible local root sites.
        middle={c:[0]*(M+1) for c in range(m)}
        for c,weighted in middle.items():
            for y in range(m-1):
                r=self.root[y*m+m-1]
                for j,a in enumerate(self.rank(m-1,c,y)):
                    weighted[j]+=a*r
        weights={c:[sum(G[i][j]*r for j,r in enumerate(hist) if r)
                    for i in range(M+1)] for c,hist in middle.items()}
        return {(x,y):sum(a*weights[c][i] for c in range(m) for i,a in enumerate(self.rank(x,y,c)) if a)
                for x in range(m-1) for y in range(m)}


def frame(u,v,m):
    """Two adjacent winning blocks, transformed to a vertical pairing."""
    n=2*m
    x,y,z,w=u%n,u//n,v%n,v//n
    if x//m==z//m and y//m!=w//m:
        x,y,z,w=y,x,w,z
    if x//m==z//m or y//m!=w//m:
        return None
    return ((min(x,n-1-x),min(y,n-1-y)),(min(z,n-1-z),min(w,n-1-w)))


def counts(book,old):
    m,M=book.m,book.M
    n=2*m
    pair,chain=book.pair(),book.chain()
    P=comb(4*M,2*M)
    single=comb(4*M,M)
    safe=lambda xy:max(xy)<m-1
    old_rows={tuple(r['r']):r for r in old['strata']}
    added={}
    checked=0
    for u,v in combinations(range(n*n),2):
        x,y,z,w=u%n,u//n,v%n,v//n
        diag=x//m!=z//m and y//m!=w//m
        distance=abs(x-z)+abs(y-w)
        r=frame(u,v,m)
        seam=linear=corner=0
        if r is not None:
            xy,zw=r
            l,r0=xy[0]==m-1,zw[0]==m-1
            if l or r0:
                ell_l=m-zw[0]+abs(xy[1]-zw[1])
                ell_r=m-xy[0]+abs(xy[1]-zw[1])
                if distance>1:
                    seam=pair[xy]*pair[zw]*safe_shuffle(M,ell_l if l else 0,ell_r if r0 else 0)
                    # The old surviving roots violate its only possible axis;
                    # neither diagonal dominant winner nor transpose applies.
                    assert (u,v) not in old_rows
            else:
                archived=old_rows[(u,v)]
                assert int(archived['vertical'] if y//m==w//m else archived['horizontal'])==P*pair[xy]*pair[zw]
                checked+=1
                if safe(zw):linear+=single*chain[xy]*book.root[zw[1]*m+zw[0]]
                if safe(xy):linear+=single*chain[zw]*book.root[xy[1]*m+xy[0]]
        elif diag:
            p=(min(x,n-1-x),min(y,n-1-y));q=(min(z,n-1-z),min(w,n-1-w))
            if (p==(m-1,m-1) and safe(q)) or (q==(m-1,m-1) and safe(p)):
                other=q if p==(m-1,m-1) else p
                corner=single*triplet(M-1,M,M)*book.root[0]**3*book.root[other[1]*m+other[0]]
                assert (u,v) not in old_rows
        if seam or linear or corner:
            assert sum(c>0 for c in (seam,linear,corner))==1
            added[u,v]={'peaks':[u,v],'distance':distance,'central_seam_winner':str(seam),
                        'suppression_chain':str(linear),'diagonal_corner_winner':str(corner),'additive_count':str(seam+linear+corner)}
    return added,checked


def main():
    p=ArgumentParser()
    p.add_argument('output',type=Path)
    a=p.parse_args()
    old=read(R51)
    gates={'chain_three_color_complete_words':[gate_chain(M) for M in (1,2,3,4)],
           'paired_two_color_complete_words':[gate_pair(M) for M in (1,2,3)]}
    book=Book()
    added, crosscheck=counts(book,old)
    upper=int(old['computed_T_upper'])
    rows=[]
    fields=('central_seam_winner','suppression_chain','diagonal_corner_winner','additive_count')
    for t in old['thresholds']:
        d=t['strict_threshold']
        total={k:sum(int(v[k]) for v in added.values() if v['distance']>d) for k in fields}
        assert sum(total[k] for k in fields[:3])==total['additive_count']
        before=int(t['sum_parts']['total'])
        after=before+total['additive_count']
        z=Fraction(after,upper)
        rows.append({'strict_threshold':d,'new_counts':{k:str(v) for k,v in total.items()},'old_R51':str(before),
                     'combined':str(after),'conditional_lower_r51_T4_exact':str(z),
                     'conditional_lower_r51_T4_decimal':dec(z),'numerator_gain':dec(Fraction(after,before))})
    out={'scope':'Three mutually disjoint finite additions to the r51 n12 probability numerator',
         'source_R51_sha256':digest(R51),'source_R50_sha256':{kind:digest(path) for kind,path in book.paths.items()},
         'maximum_rss_process_budget_mib':128,'new_independent_set_queries':0,
         'old_axis_single_overlap_crosschecked_peak_pairs':crosscheck,'gates':gates,'old_R51_T4_A2_upper':str(upper),
         'thresholds':rows,'new_strata':list(added.values())}
    with a.output.open('x') as f:
        json.dump(out,f,indent=2);f.write('\n')
    for r in rows:print(r['strict_threshold'],r['conditional_lower_r51_T4_decimal'],r['numerator_gain'])

if __name__=='__main__':main()
