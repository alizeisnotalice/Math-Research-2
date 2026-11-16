#!/usr/bin/env python3
"""Exact new controls: core coverage and admissible residual credit inside J."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dist(p,q):return max(abs(p[0]-q[0]),abs(p[1]-q[1]))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();d=json.loads(a.input.read_text());n=d['n'];h=F(d['spacing']);w=F(d['width']);oldmass=F(d['old_atom_mass'])
    assert n==5 and h==F(6,25) and w==F(24,25) and oldmass==F(1,25)
    S=list(product([i*h for i in range(n)],repeat=2));Q=[tuple(map(F,q['position'])) for q in d['outside']];m=[oldmass]*len(S)+[F(q['mass']) for q in d['outside']];P=S+Q
    assert Q==[(F(-1,20),F(12,25)),(F(6,5),F(12,25))] and m[-2:]==[F(1,50)]*2
    axes=[sorted({p[k] for p in P}) for k in (0,1)];rows={};rects=0
    for (l,r),(b,t) in product(combinations_with_replacement(axes[0],2),combinations_with_replacement(axes[1],2)):
        rects+=1;ids=tuple(i for i,p in enumerate(P) if l<=p[0]<=r and b<=p[1]<=t)
        if not ids or ids in rows:continue
        bb=(min(P[i][0] for i in ids),max(P[i][0] for i in ids),min(P[i][1] for i in ids),max(P[i][1] for i in ids));mass=sum((m[i] for i in ids),F(0));span=max(bb[1]-bb[0],bb[3]-bb[2]);rows[ids]=(bb,mass,span**2,mass>span**2)
    eligible=[(ids,r) for ids,r in rows.items() if r[3]];L=max(r[1] for ids,r in eligible);maxima=[ids for ids,r in eligible if r[1]==L]
    assert L==1 and maxima==[tuple(range(25))]
    adjacent=[set() for _ in P]
    for ids,r in eligible:
        for i in ids[1:]:adjacent[ids[0]].add(i);adjacent[i].add(ids[0])
    reached={0};stack=[0]
    while stack:
        fresh=adjacent[stack.pop()]-reached;reached|=fresh;stack.extend(fresh)
    assert len(reached)==len(P)
    total=sum(m,F(0));span=max(axes[0][-1]-axes[0][0],axes[1][-1]-axes[1][0]);gap=span**2-F(35,33)*total
    assert gap==F(1213,2640)>0
    # Direct analytic maximum-mass proof controls; no need for2^27 subsets.
    assert F(24,25)+F(2,50)==1 and dist(Q[0],Q[1])**2>1
    assert (w-Q[0][0])**2>F(51,50) and Q[1][0]**2>F(51,50) and (Q[1][0]-Q[0][0])**2>F(26,25)
    e=F(1,25);J=(-e,F(1),-e,F(1));assert all(not(J[0]<=q[0]<=J[1] and J[2]<=q[1]<=J[3]) for q in Q)
    assert dist(Q[0],Q[1])**2>F(1,25) # residual pair has empty receiver
    den=1<<d['sqrt_bits'];k=isqrt(2*den*den);rl,ru=F(k,den),F(k+1,den);assert rl**2<=2<ru**2
    creditlo,credithi=(10-ru)/1000,(10-rl)/1000;assert creditlo>F(8,1000)
    covers=[]
    for k in d['coverage_control_n']:
        assert k>=5 and 2*k*k>(k+2)**2 and k>2
        covers.append(dict(n=k,strict_pair_gap_squared_comparison_margin=2*k*k-(k+2)**2,two_by_two_side_over_gap=k-1))
    holes=[]
    for item in d['old_halo_holes']:
        nn=item['n'];hh=F(nn+1,nn*nn);pp=list(product([i*hh for i in range(nn)],repeat=2));x=tuple(map(F,item['receiver']));dd=Counter(dist(q,x) for q in pp);count=0;slacks=[]
        for R,c in sorted(dd.items()):
            count+=c;slack=4*R*R-F(count,nn*nn);assert slack>0;slacks.append(slack)
        assert -F(1,nn*nn)<=x[0]<=1 and -F(1,nn*nn)<=x[1]<=1
        holes.append(dict(n=nn,receiver=list(map(str,x)),closed_prefix_events=len(dd),strict_minimum_nonhalo_prefix_slack=str(min(slacks))))
    result=dict(status='PASS_NEW_EXACT_CORE_CREDIT_AND_GRID_COVERAGE_CONTROLS',source_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,actual_positive_atoms=len(P),direct_coordinate_rectangles=rects,complete_bbox_closures=len(rows),eligible_closures=len(eligible),maximum_mass=str(L),unique_maximum_S=True,connected=True,nonterminal_gap=str(gap),residual_positive_atoms=2,residual_eligible_labels='two singletons only',outside_measure_of_J='0',residual_halo_intersection_J_exact='(10-sqrt(2))/1000',residual_halo_intersection_J_lower=str(creditlo),residual_halo_intersection_J_upper=str(credithi),coverage_controls=covers,old_internal_halo_holes=holes,elapsed_seconds=time.perf_counter()-start,scope='New exact admissible27atom example proving H(nu) meets J despite nu(J)=0, plus strict filling inequalities for n>=5 and new point holes for n3,4. No complete Delta/outer-gap area calculated and no counterexample to3.5 claimed.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],sources=len(P),closures=len(rows),interior_residual_credit='(10-sqrt(2))/1000',seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
