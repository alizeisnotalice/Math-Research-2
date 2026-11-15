#!/usr/bin/env python3
"""Exact63-subset/whole-patch audit; deliberately no full halo-area calculation."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def certify(data):
    p=[tuple(map(F,q))for q in data['points']];m=list(map(F,data.get('initial_masses',data.get('masses'))));S=data['S_mask']
    assert p==[(F(3,10),F(3,10)),(F(6,5),F(6,5)),(F(3,10),F(6,5)),(F(6,5),F(3,10)),(F(-1,50),F(3,10)),(F(7,50),F(7,50))]
    assert m==[F(1,10),F(449,500),F(1,1000),F(1,1000),F(1,10),F(1,1000)]and S==15
    n=len(p);rows=[];closed={};parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1];mass=sum(m[i]for i in ids)
        bounds=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids))
        d2=max(bounds[1]-bounds[0],bounds[3]-bounds[2])**2;closure=sum(1<<i for i,(x,y)in enumerate(p)if bounds[0]<=x<=bounds[1]and bounds[2]<=y<=bounds[3])
        row=dict(mask=mask,mass=str(mass),span_squared=str(d2),eligible=d2<mass,closure=closure);rows.append(row)
        if closure==mask:closed[mask]=row
        if d2<mass:
            for i in ids[1:]:parent[find(i)]=find(ids[0])
    L=max(F(v['mass'])for v in rows if v['eligible']);winners=[v['mask']for v in rows if v['eligible']and F(v['mass'])==L]
    assert L==1 and winners==[15]and len({find(i)for i in range(n)})==1
    assert (F(6,5)-F(3,10))**2>F(9,16)
    d2=F(61,50)**2;margin=d2-F(35,33)*sum(m);assert margin>0
    h=F(1,10000);c=F(7,50);P=(c-h,c+h,c-h,c+h)
    def far(point):return max(abs(point[0]-P[0]),abs(point[0]-P[1]),abs(point[1]-P[2]),abs(point[1]-P[3]))
    def near(point):return max(P[0]-point[0],F(0),point[0]-P[1],P[2]-point[1],point[1]-P[3])
    N=49;B=N&~S;nr=32;r=F(161,1000);mn=sum(m[i]for i in range(n)if N>>i&1);bm=sum(m[i]for i in range(n)if B>>i&1)
    assert N in closed and closed[N]['eligible']and mn==F(201,1000)and bm==F(101,1000)
    assert all(far(p[i])<r for i in (0,4,5))and 4*r*r<mn
    assert all(4*near(p[i])**2>sum(m)for i in (1,2,3))
    assert 4*near(p[4])**2>bm and 4*far(p[5])**2<m[5]
    # H(S) has no witness here: p0 alone fails, any other S atom cannot
    # be captured even if all S mass is credited.
    assert 4*near(p[0])**2>m[0]and all(4*near(p[i])**2>L for i in (1,2,3))
    area=(P[1]-P[0])*(P[3]-P[2]);assert area==F(1,25000000)
    costs=[abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1])-abs(p[i][0]-p[j][0])*abs(p[i][1]-p[j][1])for i in range(4)for j in range(i)]
    assert min(costs)==F(9,10)and (1+F(9,10))**2==F(361,100)>F(7,2)
    assert sum(abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1])-abs(p[i][0]-p[j][0])*abs(p[i][1]-p[j][1])for i,j in ((0,2),(2,1),(1,3)))==F(27,10)>F(5,2)
    return dict(status='PASS_R68_REACHABLE_FULL_SOURCE_OLD_HALO_PATCH',sources=n,all_nonempty_subsets=rows,closed_bbox_count=len(closed),eligible_closures=sum(r['eligible']for r in closed.values()),L=str(L),maximum_masks=winners,S_rank=4,both_axes_fat=True,connected=True,nonterminal_margin=str(margin),patch=list(map(str,P)),area=str(area),support_MST_cost='27/10',support_tree_bound='37/10',support_union_area='361/100',actual_full_source_mask=N,actual_residual_source_mask=nr,full_source_minus_S_mask=B,common_witness_radius=str(r),common_witness_area=str(4*r*r),full_source_mass=str(mn),full_source_residual_mass=str(bm),minimum_q_capture_area=str(4*near(p[4])**2),minimum_p0_capture_area=str(4*near(p[0])**2),minimum_far_S_capture_area=str(min(4*near(p[i])**2 for i in (1,2,3))),patch_in_Wprime_of_actual_N_and_Hremaining=True,patch_disjoint_from_H_of_S=True,scope='Exact rational whole-patch proof, not sampled receivers or full halo area. Actual full N reachability alone cannot remove old-halo overpayment; actual residual N is a proper subset of Nfull minus S.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert sys.flags.optimize==0 and Path('/Volumes/data').is_mount();out=args.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();r=certify(json.loads(args.input.read_text()));r.update(entry_sha256=sha(Path(__file__)),input_sha256=sha(args.input),python=sys.version,pid=os.getpid(),elapsed_seconds=time.perf_counter()-start);out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k!='all_nonempty_subsets'},indent=2))
if __name__=='__main__':main()
