"""Disjoint odd-quadrant family in which one winner suppresses both losers."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb,factorial
import json
from pathlib import Path

from build_prefixes import HERE, adjacency, bits, canonical, neighbors
from merge_n11 import PROJECT, W, as_decimal, coordinates, digest, read, root_tables

PREVIOUS=PROJECT/'certificates/r048/n10-gluing'


def double_truncated(w,h,K,root_below,root_right,lookup):
    area=w*h
    adj=adjacency(w,h)
    bottom=((1<<(w-1))-1)<<(w*(h-1))
    right=sum(1<<(w*y+w-1) for y in range(h-1))
    out={}
    @lru_cache(maxsize=None)
    def sum_below(S):return sum(root_below[v%w] for v in bits(S))
    @lru_cache(maxsize=None)
    def sum_right(S):return sum(root_right[v//w] for v in bits(S))
    for y in range(h-1):
        for x in range(w-1):
            u=w*y+x
            counts=[0]*(area+1)
            current={1<<u:1}
            for j in range(2,K+1):
                nxt=defaultdict(int)
                for S,number in current.items():
                    for v in bits(neighbors(adj,S)&~S):
                        bit=1<<v
                        T=S|bit
                        if j<K:nxt[T]+=number
                        if bit&bottom and S&right:
                            weight=root_below[v%w]*sum_right(S&right)
                        elif bit&right and S&bottom:
                            weight=root_right[v//w]*sum_below(S&bottom)
                        else:
                            continue
                        counts[j]+=number*lookup(w,h,T)*weight
                current=nxt
            out[u]=counts
    return out


def double_full_square(m,root_below,root_right):
    book=read(PREVIOUS/f'double-m{m}.json')
    out={}
    for row in book['roots']:
        hist=[0]*(m*m+1)
        for a,ac in enumerate(row['last_connector_rank_counts']):
            for b,abc in enumerate(ac):
                weight=root_below[a]*root_right[b]
                for j,value in enumerate(abc,1):hist[j]+=int(value)*weight
        out[row['root']]=hist
    return out


def gates():
    small=read(HERE/'small-results/summary.json')
    counts={(q['w'],q['h'],int(q['seen'])):int(q['completion']) for q in small['records']}
    def lookup(w,h,T):return counts[canonical(w,h,T)]
    a,b=[2,3],[5,7]
    rank=double_truncated(3,3,4,a,b,lookup)
    old=double_full_square(3,a,b)
    for u,vec in rank.items():assert vec[:5]==old[u][:5],u
    return len(rank)


def multinomial(a,b,c):return factorial(a+b+c)//(factorial(a)*factorial(b)*factorial(c))


def main():
    p=ArgumentParser()
    p.add_argument('axis_ledger',type=Path)
    p.add_argument('completion_ledger',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    gate=gates()
    q=read(a.axis_ledger)
    archive=read(a.completion_ledger)
    assert digest(a.completion_ledger)==q['completion_source_sha256']
    counts={(r['w'],r['h'],int(r['seen'])):int(r['completion']) for r in archive['records']}
    @lru_cache(maxsize=None)
    def lookup(w,h,T):return counts[canonical(w,h,T)]
    root=root_tables()
    value={}
    dims={(w,h):K for w,h,K in q['all_permitted_shapes']}
    for (w,h) in root:
        right=11-w;below=11-h
        Rdown=root[w,below][:w-1]
        Rright=[root[right,h][y*right] for y in range(h-1)]
        H=double_full_square(5,Rdown,Rright) if w==h==5 else double_truncated(w,h,dims[w,h],Rdown,Rright,lookup)
        for u,ranks in H.items():
            value[w,h,u]=sum(r*multinomial(w*h-j,w*below,right*h) for j,r in enumerate(ranks) if r)
    extra=[]
    for u,v in combinations(range(121),2):
        c,r,x,y=coordinates(u);cc,rr,xx,yy=coordinates(v)
        if c==cc or r==rr:continue
        w,h,ww,hh=W[c],W[r],W[cc],W[rr]
        if x>=w-1 or y>=h-1 or xx>=ww-1 or yy>=hh-1:continue
        U=y*w+x;V=yy*ww+xx
        new=(comb(121,ww*hh)*value[w,h,U]*root[ww,hh][V]
             +comb(121,w*h)*root[w,h][U]*value[ww,hh,V])
        if new:
            d=abs(u%11-v%11)+abs(u//11-v//11)
            extra.append({'peaks':[u,v],'manhattan':d,'disjoint_words':str(new)})
    upper=int(q['A2_upper_r49'])
    terms=[]
    for row in q['thresholds']:
        th=Fraction(row['threshold'])
        extra_count=sum(int(r['disjoint_words']) for r in extra if r['manhattan']>th)
        total=extra_count+int(row['union'])
        assert total<=upper
        ratio=Fraction(total,upper)
        terms.append({'threshold':str(th),'disjoint_extra':str(extra_count),'axis_union':row['union'],'combined':str(total),
                      'conditional_lower_exact':str(ratio),'conditional_lower_decimal':as_decimal(ratio),
                      'factor_over_r49_numerator':as_decimal(Fraction(total,int(q['incumbent_word_count'])))})
    out={'n':11,'scope':'Single dominant among diagonal winners, disjoint from the two pairings',
         'n3_weighted_double_vectors_matched':gate,'source_double5_sha256':digest(PREVIOUS/'double-m5.json'),
         'axis_ledger_sha256':digest(a.axis_ledger),'A2_upper_r49':str(upper),'additional_root_pairs':len(extra),
         'unique_existing_completion_queries_accessed':lookup.cache_info().misses,'additional_ie_queries':0,
         'thresholds':terms,'additional_strata':extra}
    with a.output.open('x') as f:
        json.dump(out,f,indent=2);f.write('\n')
    for t in terms:print(t['threshold'],t['conditional_lower_decimal'],t['factor_over_r49_numerator'])

if __name__=='__main__':main()
