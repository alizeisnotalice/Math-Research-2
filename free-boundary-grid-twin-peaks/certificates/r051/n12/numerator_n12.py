"""Reuse exact size-six connector completions for the size-twelve numerator."""
from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, factorial, prod
import json
from pathlib import Path
import sys

from count_rectangles import HERE, PROJECT, dec, digest

R50=PROJECT/'certificates/r050/n11-gluing'
ROOT6=PROJECT/'computation/results/independent-peaks-n6-r004.json'
sys.path.insert(0,str(R50))
from add_dominant import double_truncated
from build_prefixes import canonical
from merge_n11 import prefix_rank_counts, split_merge_factor


def read(p):return json.loads(p.read_text())


def main():
    parser=ArgumentParser()
    parser.add_argument('denominator',type=Path)
    parser.add_argument('output',type=Path)
    a=parser.parse_args()
    N=12
    m=6
    M=m*m
    raw_plan=R50/'plan-K7.json'
    raw_completions=R50/'results-K7/summary.json'
    source=read(raw_plan)
    complete=read(raw_completions)
    assert digest(raw_plan)==complete['plan_sha256']
    lookup_table={(6,6,int(r['seen'])):int(r['completion']) for r in complete['records'] if r['w']==r['h']==6}
    rows=[r for r in source['rows'] if r['w']==r['h']==6]
    C={r['local_rank'] for r in rows}
    cutoff=next(k for w,h,k in source['shapes'] if w==h==6)
    root=list(map(int,read(ROOT6)['single_peak_root_counts']))
    assert len(root)==M and all(root[y*m+x]==root[x*m+y] for x in range(m) for y in range(m))
    ranks=prefix_rank_counts({'rows':rows},lookup_table,{(6,6):root})
    @lru_cache(maxsize=None)
    def lookup(w,h,T):return lookup_table[canonical(w,h,T)]
    weighted_double=double_truncated(m,m,cutoff,root[:m-1],[root[m*y] for y in range(m-1)],lookup)
    pair={u:sum(val*root[c]*comb(2*M-j,M) for c in range(m)
                for j,val in enumerate(ranks[m,m,u,c]) if val)
          for u in range(M) if u%m<m-1}
    dom={u:sum(val*factorial(3*M-j)//prod((factorial(M-j),factorial(M),factorial(M)))
               for j,val in enumerate(row) if val)
         for u,row in weighted_double.items()}
    # The division must be integer at each factor, independent of the rank count.
    for j in range(1,M+1):
        assert factorial(3*M-j)% (factorial(M-j)*factorial(M)**2)==0
    factor=split_merge_factor(M,M,M,M)
    pair_shuffle=comb(4*M,2*M)
    triple_shuffle=comb(4*M,M)
    strata=[]
    for u,v in combinations(range(N*N),2):
        x,y,z,w=u%N,u//N,v%N,v//N
        xc,yc,zc,wc=(min(k,N-1-k) for k in (x,y,z,w))
        vert=x//m!=z//m and xc<m-1 and zc<m-1
        horiz=y//m!=w//m and yc<m-1 and wc<m-1
        if not vert and not horiz: continue
        a0,b0=yc*m+xc,wc*m+zc
        V=pair_shuffle*pair[a0]*pair[b0] if vert else 0
        H=pair_shuffle*pair[xc*m+yc]*pair[zc*m+wc] if horiz else 0
        J=D=0
        if vert and horiz:
            r1,r2=ranks[m,m,a0,m-1],ranks[m,m,b0,m-1]
            assert r1==ranks[m,m,xc*m+yc,m-1] and r2==ranks[m,m,zc*m+wc,m-1]
            J=root[0]**2*sum(p*q*factor[i][j] for i,p in enumerate(r1) if p for j,q in enumerate(r2) if q)
            assert 0<J<=min(H,V)
            D=triple_shuffle*(dom[a0]*root[b0]+root[a0]*dom[b0])
        d=abs(x-z)+abs(y-w)
        strata.append({'r':[u,v],'distance':d,'vertical':str(V),'horizontal':str(H),'overlap':str(J),
                       'dominant_disjoint':str(D),'total':str(V+H-J+D)})
    upper=read(a.denominator)['computed_results']
    T=int(upper['three_block_T']['top_three_upper'])
    Q=int(upper['quadrants']['top_three_upper'])
    aligned=(sum(root[x]**2 for x in range(m//2)))**2*factorial(4*M)//factorial(M)**4
    terms=[]
    for th in (4,6,7):
        sub=[r for r in strata if r['distance']>th]
        totals={k:sum(int(r[k]) for r in sub) for k in ('vertical','horizontal','overlap','dominant_disjoint','total')}
        assert totals['total']==totals['vertical']+totals['horizontal']-totals['overlap']+totals['dominant_disjoint']
        assert totals['vertical']==totals['horizontal']
        if th in (4,6):assert totals['vertical']>=aligned
        val=totals['total']
        terms.append({'strict_threshold':th,'sum_parts':{k:str(v) for k,v in totals.items()},'chosen_peak_pairs':len(sub),
                      'conditional_lower_T_exact':str(Fraction(val,T)),'conditional_lower_T':dec(Fraction(val,T)),
                      'conditional_lower_four_blocks':dec(Fraction(val,Q)),'numerator_over_aligned':dec(Fraction(val,aligned)),
                      'unconditional_numerator':dec(Fraction(val,factorial(144)))})
    out={'n':12,'scope':'Strict finite original-conditional lower, reused size-six prefix queries only',
         'root_vector_source_sha256':digest(ROOT6),'K7_plan_sha256':digest(raw_plan),
         'K7_completions_sha256':digest(raw_completions),'denominator_sha256':digest(a.denominator),
         'new_independent_set_queries':0,'eligible_local_ranks':sorted(C),'num_cached_prefixes':len(lookup_table),
         'aligned_numerator':str(aligned),'aligned_conditional_lower_T':dec(Fraction(aligned,T)),
         'computed_T_upper':str(T),'computed_four_blocks_upper':str(Q),'conditional_improvement_from_T_vs_four':dec(Fraction(Q,T)),
         'thresholds':terms,'strata':strata}
    with a.output.open('x') as f:
        json.dump(out,f,indent=2);f.write('\n')
    for r in terms:print(r['strict_threshold'],r['conditional_lower_T'],r['numerator_over_aligned'])

if __name__=='__main__':main()
