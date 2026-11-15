#!/usr/bin/env python3
"""Independent exact whole-patch obstruction audit using all internal subsets.

Constrained-center minimum covering radii replace the author's staged
column exclusion. No complete halo area or producer import.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(data):
    p=[tuple(map(F,q))for q in data['points']];m=list(map(F,data['initial_masses']));S=data['S_mask'];n=len(p)
    assert p==[(F(0),F(0)),(F(99,100),F(99,100)),(F(0),F(99,100)),(F(99,100),F(0)),(F(-3,10),F(0))]
    assert m==[F(2499,5000),F(1,2),F(1,10000),F(1,10000),F(1,2)]and S==15
    rows=[];bbox={};parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1];mass=sum(m[i]for i in ids);b=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids));bbox[mask]=(b,mass)
        d2=max(b[1]-b[0],b[3]-b[2])**2;eligible=d2<mass;rows.append(dict(mask=mask,mass=str(mass),span_squared=str(d2),eligible=eligible))
        if eligible:
            for i in ids[1:]:parent[find(i)]=find(ids[0])
    L=max(F(v['mass'])for v in rows if v['eligible']);maximum=[v['mask']for v in rows if v['eligible']and F(v['mass'])==L]
    assert L==1 and maximum==[15]and len({find(i)for i in range(n)})==1
    patches=[tuple(map(F,q))for q in data['patches']];r=F(data['common_half_side']);assert r==F(4999,10000)and len(patches)==2
    # Exact constrained-center minimax radius over a closed center rectangle.
    # Intersecting each admissible center interval with [lo,hi] gives this max.
    def minradius(b,P):return max(F(0),(b[1]-b[0])/2,(b[3]-b[2])/2,b[1]-P[1],P[0]-b[0],b[3]-P[3],P[2]-b[2])
    area=F(0);certs=[];witness_mass=m[0]+m[4];assert 4*r*r<witness_mass
    for P in patches:
        assert P[0]<P[1]and P[2]<P[3]
        for i in (0,4):assert max(abs(p[i][0]-P[0]),abs(p[i][0]-P[1]),abs(p[i][1]-P[2]),abs(p[i][1]-P[3]))<=r
        excluded=[]
        for mask in range(1,S+1):
            b,mass=bbox[mask];rm=minradius(b,P);slack=4*rm*rm-mass;assert slack>0
            excluded.append(dict(internal_subset_mask=mask,minimum_radius=str(rm),capture_cost_minus_mass=str(slack)))
        rr=minradius(bbox[16][0],P);rslack=4*rr*rr-m[4];assert rslack>0
        pa=(P[1]-P[0])*(P[3]-P[2]);area+=pa;certs.append(dict(patch=list(map(str,P)),area=str(pa),all_15_internal_subsets_excluded=excluded,residual_exclusion_slack=str(rslack)))
    assert patches[1][3]<patches[0][2]and area==F(5105041,25000000)>F(1,5)
    margin=F(129,100)**2-F(35,33)*sum(m);assert margin==F(8051,110000)>0
    assert F(99,100)>F(3,4)and F(297,100)>F(5,2)
    return dict(status='PASS_R69_INDEPENDENT_PURE_MIXED_WHOLE_PATCH',all_nonempty_subsets=rows,maximum_masks=maximum,L='1',S_rank=4,both_axes_fat=True,connected=True,nonterminal_margin=str(margin),MST_cost='297/100',common_witness_radius=str(r),common_witness_mass=str(witness_mass),common_witness_slack=str(witness_mass-4*r*r),patch_certificates=certs,union_patch_area=str(area),strict_excess_over_Lover5=str(area-F(1,5)),scope='Root exact construction independently certified via all15 internal subsets and constrained-center radii on both entire rectangles. Refutes pure P<=L/5 only; no net P-credit/Delta/Wprime3.5 or C2 refutation; no full halo area.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not sys.flags.optimize and Path('/Volumes/data').is_mount()
    out=a.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    t=time.perf_counter();d=run(json.loads(a.input.read_text()));d.update(entry_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,pid=os.getpid(),elapsed_seconds=time.perf_counter()-t);out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({k:v for k,v in d.items()if k not in ('all_nonempty_subsets','patch_certificates')},indent=2))
if __name__=='__main__':main()
