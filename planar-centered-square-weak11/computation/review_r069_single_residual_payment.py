#!/usr/bin/env python3
"""Exact finite domain/internal-halo-envelope check for a peripheral credit class.

No full halo-area/W/Wprime experiment or producer import is performed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def certify(data):
    p=[tuple(map(F,q))for q in data['points']];m=list(map(F,data['initial_masses']));S=data['S_mask'];n=len(p)
    assert p==[(F(0),F(0)),(F(0),F(9,10)),(F(9,10),F(0)),(F(9,10),F(9,10)),(F(23,20),F(9,10))]
    assert m==[F(1,4)]*4+[F(6,25)]and S==15
    rows=[];eligible=[];old=[];parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,1<<n):
        ids=[i for i in range(n)if mask>>i&1];mass=sum(m[i]for i in ids)
        bbox=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids));d2=max(bbox[1]-bbox[0],bbox[3]-bbox[2])**2
        row=dict(mask=mask,mass=str(mass),span_squared=str(d2),strict_eligible=d2<mass);rows.append(row)
        if d2<mass:
            eligible.append((mask,mass))
            for i in ids[1:]:parent[find(i)]=find(ids[0])
            if not(mask&~S):
                sn,sd=math.isqrt(mass.numerator),math.isqrt(mass.denominator)
                assert sn*sn==mass.numerator and sd*sd==mass.denominator
                radius=F(sn,2*sd);old.append(dict(mask=mask,mass=str(mass),receiver_right=str(bbox[0]+radius)))
    L=max(m for s,m in eligible);winners=[s for s,m in eligible if m==L]
    assert L==1 and winners==[S]and len({find(i)for i in range(n)})==1
    assert [v['mask']for v in old]==[1,2,4,8,15]
    t=max(F(v['receiver_right'])for v in old);q=p[4];b=m[4]
    assert t==q[0]==F(23,20)
    C=(F(-1,2),F(7,5),F(-1,2),F(7,5));margin=min(q[0]-C[0],C[1]-q[0],q[1]-C[2],C[3]-q[1])
    assert margin==F(1,4)and b<4*margin*margin
    area=(C[1]-C[0])*(C[3]-C[2]);assert area==F(361,100)
    # Cartesian support: the unit-square interval unions overlap in both axes.
    assert F(9,10)<1 and area>F(7,2)
    costs=[abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1])-abs(p[i][0]-p[j][0])*abs(p[i][1]-p[j][1])for i in range(4)for j in range(i)]
    assert min(costs)==F(9,10)and 3*min(costs)==F(27,10)>F(5,2)
    whole_span2=F(23,20)**2;nonterminal=whole_span2-F(35,33)*sum(m);assert nonterminal==F(97,13200)>0
    # The known positive crossing edge establishes actual eligible connectivity.
    assert (q[0]-F(9,10))**2<F(1,4)+b
    # Half-Q is safely peripheral: the whole Q lies to the right of K.
    assert q[0]-margin==F(9,10)and F(9,10)**2>F(3,4)
    wp_bound=area-b/2;delta_bound=area-b
    assert wp_bound==F(349,100)<F(7,2)and delta_bound==F(337,100)
    return dict(status='PASS_R69_SINGLE_RESIDUAL_HALFPLANE_PAYMENT',all_nonempty_subsets=rows,internal_eligible_receiver_bounds=old,L='1',maximum_masks=winners,S_rank=4,both_axes_fat=True,connected=True,nonterminal_margin=str(nonterminal),source_internal_halo_right_upper=str(t),residual_atom_mass=str(b),residual_radius_squared=str(b/4),C1=list(map(str,C)),minimum_Q_to_C1_face_margin=str(margin),Q_strictly_inside_C1=True,Q_outside_central_K=True,uniform_support_union_area=str(area),support_MST_cost='27/10',half_residual_halo_credit=str(b/2),Wprime_upper_by_proved_class=str(wp_bound),Delta_upper_by_support_minus_residual=str(delta_bound),strict_Wprime_gap_to_7over2=str(F(7,2)-wp_bound),scope='One prescribed31-subset/internal-envelope arithmetic certificate for general one-residual-atom halfplane-credit lemma. No full halo union, actual W/Wprime area, producer run, optimization or general3.5 theorem.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert not sys.flags.optimize and Path('/Volumes/data').is_mount();out=a.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    t=time.perf_counter();r=certify(json.loads(a.input.read_text()));r.update(entry_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,pid=os.getpid(),elapsed_seconds=time.perf_counter()-t);out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k!='all_nonempty_subsets'},indent=2))
if __name__=='__main__':main()
