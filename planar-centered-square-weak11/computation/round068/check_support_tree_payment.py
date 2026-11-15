#!/usr/bin/env python3
"""Bounded exact R68 support-tree/chain checks; no W/Wprime or LP evaluation."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,hashlib,json,os,sys,time
import sympy as sp

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def span2(p):return max(max(x for x,y in p)-min(x for x,y in p),max(y for x,y in p)-min(y for x,y in p))**2
def rectangle_union_area(rects):
    # Independent rational coordinate-cell integration, not tree addition.
    xx=sorted({r[k]for r in rects for k in (0,1)});yy=sorted({r[k]for r in rects for k in (2,3)})
    area=F(0);cells=0
    for xl,xr in zip(xx,xx[1:]):
        for yl,yr in zip(yy,yy[1:]):
            x=(xl+xr)/2;y=(yl+yr)/2
            if any(l<x<r and b<y<t for l,r,b,t in rects):area+=(xr-xl)*(yr-yl);cells+=1
    return area,cells

def tree(p):
    edges=sorted((abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1])-abs((p[i][0]-p[j][0])*(p[i][1]-p[j][1])),i,j)for i,j in combinations(range(len(p)),2))
    parent=list(range(len(p)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    chosen=[]
    for c,i,j in edges:
        a,b=root(i),root(j)
        if a!=b:parent[a]=b;chosen.append((i,j,c))
    assert len(chosen)==len(p)-1
    return sum((x[2]for x in chosen),F(0)),chosen

def chain_after_reflection(p):
    hits=[]
    for sx in (-1,1):
        for sy in (-1,1):
            q=sorted((sx*x,sy*y)for x,y in p)
            if all(q[i][1]<=q[i+1][1]for i in range(len(q)-1)):hits.append((sx,sy,q))
    return hits

def analyze(c):
    p=[tuple(map(F,z))for z in c['points']];mass=list(map(F,c['masses']));n=len(p);k=c['S_count'];S=(1<<k)-1
    assert len(mass)==n and min(mass)>0 and len(set(p))==n and sum(mass[:k])==1
    eligible=[]
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1];m=sum(mass[i]for i in ids)
        if m>span2([p[i]for i in ids]):eligible.append((mask,m))
    L=max(m for mask,m in eligible);maxmasks=[mask for mask,m in eligible if m==L]
    assert L==1 and maxmasks==[S]
    # Full eligible hypergraph connectivity via union of its hyperedges.
    seen=1
    while True:
        old=seen
        for mask,m in eligible:
            if mask&seen:seen|=mask
        if seen==old:break
    assert seen==(1<<n)-1
    fullmass=sum(mass);terminal_margin=span2(p)-F(35,33)*fullmass;assert terminal_margin>0
    s=p[:k];w=max(x for x,y in s)-min(x for x,y in s);h=max(y for x,y in s)-min(y for x,y in s)
    assert w>F(3,4)and h>F(3,4)and span2(s)<1
    rects=[(x-F(1,2),x+F(1,2),y-F(1,2),y+F(1,2))for x,y in s]
    area,cells=rectangle_union_area(rects);cost,edges=tree(s);bound=1+cost;assert area<=bound
    assert str(bound)==c['expected_tree_bound']
    chains=chain_after_reflection(s);assert bool(chains)==c['expected_chain']
    chain_formula=None
    if chains:
        q=chains[0][2];product=sum((q[i+1][0]-q[i][0])*(q[i+1][1]-q[i][1])for i in range(k-1))
        chain_formula=1+w+h-product;assert area==chain_formula and area<3
    if 'expected_Z_area'in c:assert str(area)==c['expected_Z_area']
    if 'expected_tree_cost'in c:assert str(cost)==c['expected_tree_cost']
    out=dict(id=c['id'],positive_sources=n,subset_count=(1<<n)-1,strict_eligible_subsets=len(eligible),S_mask=S,S_rank=k,unique_maximum_masks=maxmasks,L='1',S_width=str(w),S_height=str(h),connected=True,nonterminal_margin_exact=str(terminal_margin),support_union_area_exact=str(area),compressed_covered_cells=cells,MST_cost_exact=str(cost),tree_upper_exact=str(bound),tree_edges=[dict(i=i,j=j,cost=str(v))for i,j,v in edges],certifies_tree_5over2_extra=cost<=F(5,2),certifies_Wprime_le_3point5_by_tree=bound<=F(7,2),chain_after_axis_reflection=bool(chains),chain_exact_formula=str(chain_formula)if chain_formula is not None else None,coarse_bbox_area=str((1+w)*(1+h)),Wprime_evaluated=False)
    if c['id']=='two-row-tree-method-barrier':
        assert area>F(7,2) and cost>F(5,2)
        out['support_union_excess_over_7over2']=str(area-F(7,2));out['MST_excess_over_5over2']=str(cost-F(5,2));out['scope']='Does not refute Wprime, whose value is not computed; refutes universality of the support-only/MST sufficient test.'
    return out

def main():
    assert os.path.ismount('/Volumes/data');start=time.perf_counter()
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    data=json.loads(args.input.read_text());s,x,y=sp.symbols('s x y',positive=True)
    identities=[sp.expand(s*s-(s-x)*(s-y)-(s*(x+y)-x*y))]
    dx=sp.symbols('dx1:4',nonnegative=True);dy=sp.symbols('dy1:4',nonnegative=True)
    identities.append(sp.expand(s*s+sum(s*s-(s-a)*(s-b)for a,b in zip(dx,dy))-(s*s+s*sum(dx+dy)-sum(a*b for a,b in zip(dx,dy)))))
    assert all(z==0 for z in identities)
    result=dict(process_pid=os.getpid(),round=68,status='PASS_BOUNDED_EXACT_SUPPORT_TREE_AND_CHAIN_CHECKS',python=sys.version,executable=sys.executable,sympy=sp.__version__,script_sha256=sha(Path(__file__)),input_sha256=sha(args.input),symbolic_zero_identities=list(map(str,identities)),examples=[analyze(c)for c in data['examples']],elapsed_seconds=time.perf_counter()-start,scope='Exact rational eligibility/connectivity/domain/tree/rectangle-union checks and symbolic algebra only. General area lemma is proved analytically. No Wprime measurement, optimization, floating LP, MATLAB or general3.5 selection certificate.')
    assert str(args.out.resolve()).startswith('/Volumes/data/');args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],python=result['python'],sympy=result['sympy'],examples=[{k:r[k]for k in ['id','S_rank','support_union_area_exact','MST_cost_exact','tree_upper_exact','certifies_Wprime_le_3point5_by_tree']}for r in result['examples']],seconds=result['elapsed_seconds']),indent=2))
if __name__=='__main__':main()
