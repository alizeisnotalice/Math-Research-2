#!/usr/bin/env python3
"""Independent127-subset/net patch audit: inclusion-exclusion and all125 trees.

No author imports, full-halo experiment, Prim/Kruskal, or coordinate-cell area.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product,combinations
import argparse,hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def intersection(rects):return (max(r[0]for r in rects),min(r[1]for r in rects),max(r[2]for r in rects),min(r[3]for r in rects))
def area(r):return max(F(0),r[1]-r[0])*max(F(0),r[3]-r[2])
def square(p,r):return (p[0]-r,p[0]+r,p[1]-r,p[1]+r)
def uncovered(J,covers):
    covered=F(0);nonzero=0
    for k in range(1,len(covers)+1):
        for selected in combinations(covers,k):
            a=area(intersection([J]+list(selected)))
            if a:nonzero+=1
            covered+=(1 if k%2 else-1)*a
    return area(J)-covered,nonzero

def certify(data):
    if 'sources'in data:
        p=[tuple(map(F,q['point']))for q in data['sources']];m=[F(q['mass'])for q in data['sources']]
    else:p=[tuple(map(F,q))for q in data['points']];m=list(map(F,data.get('initial_masses',data.get('masses'))))
    S=data['S_mask'];n=len(p)
    assert p==[(F(0),F(0)),(F(1,2),F(0)),(F(9,10),F(0)),(F(0),F(9,10)),(F(9,10),F(9,10)),(F(-4,5),F(0)),(F(13,10),F(0))]
    assert m==[F(499,1000)]*2+[F(1,1500)]*3+[F(249,500)]*2 and S==31
    rows=[];closed={};parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1];mass=sum(m[i]for i in ids);b=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids));d2=max(b[1]-b[0],b[3]-b[2])**2
        closure=sum(1<<i for i,(x,y)in enumerate(p)if b[0]<=x<=b[1]and b[2]<=y<=b[3]);row=dict(mask=mask,mass=str(mass),span_squared=str(d2),eligible=d2<mass,closure=closure);rows.append(row)
        if closure==mask:closed[mask]=row
        if d2<mass:
            for i in ids[1:]:parent[find(i)]=find(ids[0])
    L=max(F(r['mass'])for r in rows if r['eligible']);winners=[r['mask']for r in rows if r['eligible']and F(r['mass'])==L]
    assert L==1 and winners==[31]and len({find(i)for i in range(n)})==1 and len(closed)==30 and sum(r['eligible']for r in closed.values())==17
    def cost(i,j):
        x,y=abs(p[i][0]-p[j][0]),abs(p[i][1]-p[j][1]);return x+y-x*y
    tree_costs=[]
    for seq in product(range(5),repeat=3):
        degree=[1]*5
        for x in seq:degree[x]+=1
        c=F(0)
        for x in seq:
            leaf=next(i for i,d in enumerate(degree)if d==1);c+=cost(leaf,x);degree[leaf]-=1;degree[x]-=1
        i,j=[i for i,d in enumerate(degree)if d==1];c+=cost(i,j);tree_costs.append(c)
    assert len(tree_costs)==125 and min(tree_costs)==F(27,10)>F(5,2)
    E=(F(0),F(1,2),F(-1,2),F(1,2));rh=F(177,500);rt=F(23,1000);rq=F(353,1000)
    Iplus=[E,square(p[0],rh),square(p[1],rh)]+[square(p[i],rt)for i in (2,3,4)];Rplus=[square(p[i],rq)for i in (5,6)]
    assert 4*rh*rh>F(501,1000)and 4*rt*rt>F(1,500)and 4*rq*rq>F(249,500)
    cover_cases=[]
    for mask in range(1,S+1):
        mass=sum(m[i]for i in range(5)if mask>>i&1);h=(mask&3).bit_count()
        if h==2:assert mass<=1;kind='two_heavy_in_E'
        elif h==1:assert mass<=F(501,1000);kind='one_heavy_in_177over500_square'
        else:assert mass<=F(1,500);kind='tail_in_23over1000_squares'
        cover_cases.append(dict(mask=mask,mass=str(mass),cover_case=kind))
    assert (p[6][0]-p[5][0])**2>m[5]+m[6]
    for a in Iplus:
        for b in Rplus:
            q=intersection([a,b]);assert q[0]>q[1]or q[2]>q[3]
    r=F(499,1000);J=[(F(-499,1000),F(-301,1000),-r,r),(F(801,1000),F(999,1000),-r,r)]
    assert J[0][1]<J[1][0]and 4*r*r<F(997,1000)
    for box,ids in zip(J,((0,5),(1,6))):
        for i in ids:assert max(abs(p[i][0]-box[0]),abs(p[i][0]-box[1]),abs(p[i][1]-box[2]),abs(p[i][1]-box[3]))<=r
    pieces=[uncovered(box,Iplus+Rplus)for box in J];patch=sum(a for a,k in pieces)
    assert patch==F(12231,50000)>F(1,5)
    # Independent exact scalar/rectangle check of the heavy-core application.
    H=[0,1];V=[2,3,4];M=sum(m[i]for i in H);tail=1-M;rc=F(43,50)
    near=[i for i in (5,6)if min(max(abs(p[i][0]-p[j][0]),abs(p[i][1]-p[j][1]))**2 for j in V)<tail];beta=sum((m[i]for i in near),F(0))
    assert tail==F(1,500)and beta==0 and 4*rc*rc<=3*M
    C=(F(-1,2),F(7,5),F(-1,2),F(7,5));K=intersection([square(p[i],rc)for i in H]);ZH=(F(-1,2),F(1),F(-1,2),F(1,2))
    ck=area(intersection([C,K]));zk=area(intersection([ZH,K]));gamma=area(C)-ck+zk;bound=gamma+9*(tail+beta)
    assert ck==F(1037,625)and zk==F(61,50)and gamma==F(7927,2500)and bound==F(1993,625)<F(7,2)
    margin=F(21,10)**2-F(35,33)*sum(m);assert margin==F(7567,3300)>0
    return dict(status='PASS_R69_INDEPENDENT_NET_MIXED_PATCH',all_nonempty_subsets=rows,closed_bboxes=len(closed),eligible_closed_bboxes=sum(r['eligible']for r in closed.values()),maximum_masks=winners,L='1',connected=True,nonterminal_margin=str(margin),all125_support_trees_minimum=str(min(tree_costs)),internal_cover_cases=cover_cases,Iplus=list(map(lambda r:list(map(str,r)),Iplus)),Rplus=list(map(lambda r:list(map(str,r)),Rplus)),all_internal_residual_cover_pairs_strictly_disjoint=True,old_credit_exact='0',fixed_witness_radius=str(r),fixed_witness_area=str(4*r*r),fixed_witness_mass='997/1000',patch_area_by_inclusion_exclusion=str(patch),patch_terms_nonzero=[k for a,k in pieces],net_increment_lower=str(patch),strict_excess_over_one_fifth=str(patch-F(1,5)),heavy_core_application=dict(H_mask=3,M=str(M),tail_mass=str(tail),beta=str(beta),r=str(rc),C=list(map(str,C)),K=list(map(str,K)),C_intersect_K_area=str(ck),Z_H_intersect_K_area=str(zk),Gamma=str(gamma),W_upper_by_general_lemma=str(bound)),scope='Exact127-subset/source-cover/wholepatch plus2×255 inclusion-exclusion terms and125support trees. No full halo area. Refutes net P-credit<=L/5, not general3.5; separately proved heavy-core class pays this input.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not sys.flags.optimize and Path('/Volumes/data').is_mount()
    out=a.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    t=time.perf_counter();v=certify(json.loads(a.input.read_text()));v.update(entry_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,pid=os.getpid(),elapsed_seconds=time.perf_counter()-t);out.write_text(json.dumps(v,indent=2)+'\n');print(json.dumps({k:x for k,x in v.items()if k not in ('all_nonempty_subsets','internal_cover_cases','Iplus','Rplus')},indent=2))
if __name__=='__main__':main()
