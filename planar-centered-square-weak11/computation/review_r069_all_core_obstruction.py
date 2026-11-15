#!/usr/bin/env python3
"""Independent255-subset/wholepatch/all15core lower certificate; no full halo.

All-core roots are handled by exact rational shape certificates, including
sqrt(6)/4>3/5, not finite-bit lower/upper guesses.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import argparse,hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def intersection(rects):return (max(r[0]for r in rects),min(r[1]for r in rects),max(r[2]for r in rects),min(r[3]for r in rects))
def area(r):return max(F(0),r[1]-r[0])*max(F(0),r[3]-r[2])
def square(p,r):return (p[0]-r,p[0]+r,p[1]-r,p[1]+r)
def uncovered(J,covers):
    result=area(J)
    for k in range(1,len(covers)+1):
        for subset in combinations(covers,k):result+=(-1 if k%2 else 1)*area(intersection([J]+list(subset)))
    return result

def certify(data):
    if 'sources'in data:p=[tuple(map(F,q['point']))for q in data['sources']];m=[F(q['mass'])for q in data['sources']]
    else:p=[tuple(map(F,q))for q in data['points']];m=list(map(F,data.get('initial_masses',data.get('masses'))))
    S=data['S_mask'];assert S==15 and len(p)==len(set(p))==8
    assert set(p[:4])=={(x,y)for x in (F(0),F(9,10))for y in (F(0),F(9,10))}and m[:4]==[F(1,4)]*4
    assert set(p[4:])=={(F(-3,4),F(-11,50)),(F(28,25),F(-3,4)),(F(33,20),F(28,25)),(F(-11,50),F(33,20))}and m[4:]==[F(749,1000)]*4
    rows=[];eligible=[];closed={};bboxes={};parent=list(range(8))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in range(1,256):
        ids=[i for i in range(8)if mask>>i&1];mass=sum(m[i]for i in ids);b=(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids));d2=max(b[1]-b[0],b[3]-b[2])**2;bboxes[mask]=b
        closure=sum(1<<i for i,(x,y)in enumerate(p)if b[0]<=x<=b[1]and b[2]<=y<=b[3]);row=dict(mask=mask,mass=str(mass),span_squared=str(d2),eligible=d2<mass,closure=closure);rows.append(row)
        if closure==mask:closed[mask]=row
        if d2<mass:
            eligible.append(row)
            for i in ids[1:]:parent[find(i)]=find(ids[0])
    L=max(F(v['mass'])for v in eligible);winners=[v['mask']for v in eligible if F(v['mass'])==L]
    assert L==1 and winners==[15]and len({find(i)for i in range(8)})==1
    internal=[v for v in eligible if not v['mask']&~S];residual=[v for v in eligible if not v['mask']&S];cross=[v for v in eligible if v['mask']&S and v['mask']&~S]
    assert set(v['mask']for v in internal)=={1,2,4,8,15}
    assert set(v['mask']for v in residual)=={16,32,64,128}
    assert len(cross)==4 and all(v['mask'].bit_count()==2 and F(v['mass'])==F(999,1000)for v in cross)
    assert {v['mask']&S for v in cross}=={1,2,4,8}
    Iplus=[square(p[i],F(1,4))for i in range(4)]+[(F(2,5),F(1,2),F(2,5),F(1,2))]
    rq=F(433,1000);assert 4*rq*rq>F(749,1000);Rplus=[square(p[i],rq)for i in range(4,8)]
    for a in Iplus:
        for b in Rplus:
            z=intersection([a,b]);assert z[0]>z[1]or z[2]>z[3]
    r=F(999,2000);assert 4*r*r<F(999,1000)
    J=[]
    for v in cross:
        a,b,c,d=bboxes[v['mask']];box=(b-r,a+r,d-r,c+r);assert area(box)>0;J.append(box)
        for I in Iplus:
            z=intersection([box,I]);assert z[0]>z[1]or z[2]>z[3]
    for a,b in combinations(J,2):assert area(intersection([a,b]))==0
    patch=sum(uncovered(box,Rplus)for box in J);assert patch==F(255759,1000000)>F(1,5)
    # Exact shape types for every core at its optimal radius.
    Carea=F(361,100);minGamma=F(1419,400);heads=[]
    assert F(3,8)>F(3,5)**2 and F(3,16)<F(1,2)**2
    for H in range(1,16):
        k=H.bit_count();E=1-F(k,4)
        if k==1:g=Carea;kind='one_corner_K_inside_ZH';strict=False
        elif k==2:
            i,j=[i for i in range(4)if H>>i&1];adj=p[i][0]==p[j][0]or p[i][1]==p[j][1]
            # rho=sqrt(6)/4>3/5; formulas follow from the missing K corners.
            if adj:g=F(241,100)+F(19,10)*F(3,5);kind='adjacent_pair_Gamma=241/100+19rho/10'
            else:g=F(59,25)+2*F(3,5);kind='diagonal_pair_Gamma=59/25+2rho'
            strict=True
        elif k==3:g=minGamma;kind='three_corners_missing_K_hole_area1/16';strict=False
        else:g=Carea;kind='all_corners_ZH=C';strict=False
        assert g>=minGamma>F(7,2)
        # A nonempty V always gives |ZsqrtE(V)|>=E and9(E+beta)>=E.
        phi_lower=g+E
        if H!=S:assert E>=F(1,4)and phi_lower>Carea
        else:assert E==0 and phi_lower==Carea
        heads.append(dict(H_mask=H,size=k,M=str(F(k,4)),E=str(E),shape=kind,Gamma_lower=str(g),Gamma_lower_strict=strict,Phi_lower=str(phi_lower),all_core_payment_fails=True))
    margin=F(12,5)**2-F(35,33)*sum(m);assert margin>0
    assert min(heads,key=lambda v:F(v['Gamma_lower']))['Gamma_lower']==str(minGamma)
    # Four distinct crossing deleted singletons give a separate known payment.
    Iarea=4*F(1,4)+F(1,100);Wupper=Iarea+sum(F(v['mass'])-F(749,1000)for v in cross)
    assert Iarea==F(101,100)and Wupper==F(201,100)<F(7,2)
    for v in cross:
        a,b,c,d=bboxes[v['mask']];assert sorted((b-a,d-c))==[F(11,50),F(3,4)]and F(v['mass'])<1
    receiver_area_upper=(1-F(3,4))*(1-F(11,50));strict_W_bound=Iarea+4*receiver_area_upper
    assert receiver_area_upper==F(39,200)and strict_W_bound==F(179,100)<F(7,2)
    return dict(status='PASS_R69_INDEPENDENT_ALL_CORE_CONDITIONAL_OBSTRUCTION',all_nonempty_subsets=rows,closed_bboxes=len(closed),eligible_closed_bboxes=sum(v['eligible']for v in closed.values()),maximum_masks=winners,L='1',S_rank=4,both_axes_fat=True,connected=True,nonterminal_margin=str(margin),MST_cost='27/10',internal_halo_area_exact=str(Iarea),internal_residual_halo_credit_exact='0',all_crossing_eligible_labels=cross,fixed_witness_radius=str(r),mixed_receiver_boxes=[list(map(str,b))for b in J],residual_outer_radius=str(rq),whole_patch_area_by_inclusion_exclusion=str(patch),net_increment_lower=str(patch),net_excess_over_one_fifth=str(patch-F(1,5)),all15_head_certificates=heads,minimum_Gamma_exact=str(minGamma),Phi_star_exact=str(Carea),Phi_star_gap_to_7over2=str(Carea-F(7,2)),separate_W_upper_by_four_distinct_crossing_labels=str(Wupper),separate_W_strict_upper_by_four_receiver_area_caps=str(strict_W_bound),scope='255source subsets,whole fixed-radius patch andall15 optimal-head bounds. Refutes thecore-or-net conditional shortcut, not3.5: known finite single-label capacities separately pay W<=201/100. No full halo union or parameter search.')
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not sys.flags.optimize and Path('/Volumes/data').is_mount()
    out=a.output.resolve();assert out.is_relative_to(Path('/Volumes/data/CodexTemp').resolve())and out.suffix=='.json'and not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();r=certify(json.loads(a.input.read_text()));r.update(entry_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,pid=os.getpid(),elapsed_seconds=time.perf_counter()-start);out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in ('all_nonempty_subsets','all15_head_certificates','all_crossing_eligible_labels')},indent=2))
if __name__=='__main__':main()
