"""Exact two-axis union of rank-truncated odd-quadrant gluings."""
from argparse import ArgumentParser
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import hashlib
import json
from math import comb,factorial,prod
from pathlib import Path
import os

from build_prefixes import HERE, adjacency, bits, canonical, make_plan, neighbors

PROJECT=Path(os.environ.get('A022_RESEARCH_ROOT', Path(__file__).resolve().parents[3]))
N11_ARCHIVE=PROJECT/'certificates/r049/odd-quadrants/n11-numerator.json'
DENOM=PROJECT/'certificates/r049/n11-dense/verification.json'
RANK5=PROJECT/'certificates/r048/n10-gluing/ranks-m5.json'
ROOT6=PROJECT/'computation/results/independent-peaks-n6-r004.json'
N=11
W={0:5,1:6}


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def as_decimal(r):
    with localcontext() as ctx:
        ctx.prec=57
        return str(Decimal(r.numerator)/Decimal(r.denominator))


def root_tables():
    base=read(N11_ARCHIVE)
    small=read(RANK5)
    rect=list(map(int,base['rectangle_root_counts']))
    root={(5,5):list(map(int,small['rooted_counts'])),(5,6):rect,
          (6,5):[rect[x*5+y] for y in range(5) for x in range(6)],
          (6,6):list(map(int,read(ROOT6)['single_peak_root_counts']))}
    assert all(len(r)==w*h for (w,h),r in root.items())
    return root


def prefix_rank_counts(plan,completions,root):
    ranks={}
    for row in plan['rows']:
        w,h,u,c,j=(row[k] for k in ('w','h','root','connector_x','local_rank'))
        keys=(min(w,h),max(w,h))
        ans=0
        for S,wgt in row['prefixes']:
            if j==1:
                ans+=root[w,h][u]*wgt
            else:
                ans+=completions[*keys,int(S)]*wgt
        if ans:
            rank=ranks.setdefault((w,h,u,c),[0]*(w*h+1))
            rank[j]=ans
    for (w,h,u,c),rank in ranks.items():
        assert 0 < sum(rank) <= root[w,h][u]
    return ranks


def full_square5_ranks():
    book=read(RANK5)
    result={(5,5,row['root'],c):[0,*map(int,bin)] for row in book['roots']
            for c,bin in enumerate(row['rank_counts'])}
    for y in range(5):
        for c in range(5):
            result[5,5,5*y+3,c]=result[5,5,5*y+1,4-c]
    return result


def small_replays():
    plan=read(HERE/'plan-small.json')
    archive=read(HERE/'small-results/summary.json')
    assert digest(HERE/'plan-small.json')==archive['plan_sha256']
    adj=adjacency(3,3)
    full=(1<<9)-1
    @lru_cache(maxsize=None)
    def suffix(T):
        if T==full:
            return 1
        return sum(suffix(T|(1<<v)) for v in bits(neighbors(adj,T)&~T))
    known={(3,3,int(r['seen'])):int(r['completion']) for r in archive['records']}
    assert all(c==suffix(T) for (_,_,T),c in known.items())
    roots={(3,3):[suffix(1<<i) for i in range(9)]}
    ranks=prefix_rank_counts(plan,known,roots)
    brute=defaultdict(lambda:[0]*10)
    for word in permutations(range(9)):
        u=word[0]
        if u%3==2:
            continue
        seen=1<<u
        for v in word[1:]:
            if not adj[v]&seen:
                break
            seen|=1<<v
        else:
            for j,v in enumerate(word,1):
                if v//3==2 and (j<=4 or j==abs(v%3-u%3)+abs(v//3-u//3)+1):
                    brute[3,3,u,v%3][j]+=1
    assert dict(brute)==ranks
    return {'n3_queried_prefixes':len(known),'n3_exhaustive_orders':factorial(9),
            'n3_rank_vectors_checked':len(ranks)}


def split_merge_factor(M_a,M_b,L_a,L_b):
    total=M_a+M_b+L_a+L_b
    fact=[factorial(j) for j in range(total+1)]
    arr=[[0]*(M_b+2) for _ in range(M_a+2)]
    for a in range(M_a,-1,-1):
        for b in range(M_b,-1,-1):
            fact_rest=fact[total-a-b-1]
            other=fact[M_a-a]*fact[M_b-b]
            p=fact_rest//(other*fact[L_a-1]*fact[L_b])+fact_rest//(other*fact[L_a]*fact[L_b-1])
            arr[a][b]=comb(a+b,a)*p+arr[a+1][b]+arr[a][b+1]-arr[a+1][b+1]
    return arr


def interleaving_gate(a,b,c,d):
    total=a+b+c+d
    positions=set(range(total))
    counts=[[0]*(b+1) for _ in range(a+1)]
    for A in combinations(range(total),a):
        left=positions-set(A)
        for B in combinations(sorted(left),b):
            pivot=min(left-set(B))
            i=sum(x<pivot for x in A)
            j=sum(y<pivot for y in B)
            counts[i][j]+=comb(c+d,c)
    gate=split_merge_factor(a,b,c,d)
    for i in range(a+1):
        for j in range(b+1):
            assert gate[i][j]==sum(counts[x][y] for x in range(i,a+1) for y in range(j,b+1))


def form_pairs(ranks,root):
    pair_counts={}
    for w,h in root:
        other_h=N-h
        M=w*h
        L=w*other_h
        for y in range(h):
            for x in range(w-1):
                u=y*w+x
                pair_counts[w,h,u]=sum(a*root[w,other_h][c]*comb(M+L-j,L)
                                       for c in range(w) for j,a in enumerate(ranks[w,h,u,c]) if a)
    return pair_counts


def coordinates(v):
    x,y=v%N,v//N
    col,row=int(x>=5),int(y>=5)
    return col,row,(N-1-x if col else x),(N-1-y if row else y)


def union_strata(ranks,root):
    pair=form_pairs(ranks,root)
    all_merge=comb(121,55)
    out=[]
    for u,v in combinations(range(121),2):
        c,r,x,y=coordinates(u)
        c2,r2,x2,y2=coordinates(v)
        w,h,ww,hh=W[c],W[r],W[c2],W[r2]
        vertical=c!=c2 and x<w-1 and x2<ww-1
        horizontal=r!=r2 and y<h-1 and y2<hh-1
        if not (vertical or horizontal):
            continue
        col=all_merge*pair[w,h,y*w+x]*pair[ww,hh,y2*ww+x2] if vertical else 0
        row=all_merge*pair[h,w,x*h+y]*pair[hh,ww,x2*hh+y2] if horizontal else 0
        both=0
        if vertical and horizontal:
            a,b=w*h,ww*hh
            ell,ell2=w*hh,ww*h
            C=overlap_cached(a,b,ell,ell2)
            R=root[w,hh][0]*root[ww,h][0]
            U=ranks[w,h,y*w+x,w-1]
            V=ranks[ww,hh,y2*ww+x2,ww-1]
            assert U==ranks[h,w,x*h+y,h-1] and V==ranks[hh,ww,x2*hh+y2,hh-1]
            both=R*sum(U[i]*V[j]*C[i][j] for i in range(1,a+1) if U[i] for j in range(1,b+1) if V[j])
            assert 0<both<=min(row,col)
        dist=abs(u%N-v%N)+abs(u//N-v//N)
        out.append({'peaks':[u,v],'manhattan':dist,'vertical':str(col),'horizontal':str(row),
                    'intersection':str(both),'union':str(col+row-both)})
    return out


@lru_cache(maxsize=None)
def overlap_cached(a,b,ell,ell2):return split_merge_factor(a,b,ell,ell2)


def main():
    ap=ArgumentParser()
    ap.add_argument('plan',type=Path)
    ap.add_argument('completions',type=Path)
    ap.add_argument('out',type=Path)
    op=ap.parse_args()
    gate=small_replays()
    for shape in ((1,1,1,1),(2,3,1,2),(4,3,2,3)):
        interleaving_gate(*shape)
    plan=read(op.plan)
    archive=read(op.completions)
    assert digest(op.plan)==archive['plan_sha256']
    completions={(q['w'],q['h'],int(q['seen'])):int(q['completion']) for q in archive['records']}
    assert len(completions)==len(plan['tasks'])
    root=root_tables()
    ranks=prefix_rank_counts(plan,completions,root)
    ranks.update(full_square5_ranks())
    rows=union_strata(ranks,root)
    upper=int(read(DENOM)['A2_upper'])
    incumbent=int(read(N11_ARCHIVE)['constructed_words'])
    terms=[]
    for thresh in (Fraction(11,3),Fraction(11,2),Fraction(6),Fraction(7)):
        chosen=[r for r in rows if r['manhattan']>thresh]
        a,b,c,o=(sum(int(r[key]) for r in chosen) for key in ('vertical','horizontal','intersection','union'))
        assert a==b and o==a+b-c and o<=upper
        if thresh==Fraction(11,2):
            assert a>=incumbent
        ans=Fraction(o,upper)
        terms.append({'threshold':str(thresh),'root_pairs':len(chosen),'vertical':str(a),'horizontal':str(b),'intersection':str(c),
                      'union':str(o),'conditional_lower_exact':str(ans),'conditional_lower_decimal':as_decimal(ans),
                      'factor_over_r49_numerator':as_decimal(Fraction(o,incumbent))})
    output={'n':11,'scope':'Exact odd-quadrant union with rank-truncated rectangle and large square words',
            'all_permitted_shapes':plan['shapes'],'rank5_source_sha256':digest(RANK5),'root6_source_sha256':digest(ROOT6),
            'r49_numerator_source_sha256':digest(N11_ARCHIVE),'plan_source_sha256':digest(op.plan),
            'completion_source_sha256':digest(op.completions),'A2_upper_r49_sha256':digest(DENOM),'A2_upper_r49':str(upper),
            'incumbent_word_count':str(incumbent),'gates':gate,
            'thresholds':terms,'strata':rows}
    with op.out.open('x') as f:
        json.dump(output,f,indent=2);f.write('\n')
    for row in terms:print(row['threshold'],row['conditional_lower_decimal'],row['factor_over_r49_numerator'])


if __name__=='__main__':main()
