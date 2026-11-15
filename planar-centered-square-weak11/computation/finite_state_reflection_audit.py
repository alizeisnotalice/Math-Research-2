#!/usr/bin/env python3
"""Small exact counterexamples and finite-state thresholds; no template halo run."""
from fractions import Fraction as F
from math import ceil,isqrt
from pathlib import Path
import hashlib,json,os,platform,sys


def root_bracket(q,bits=60):
    D=1<<bits; t=isqrt(q.numerator*D*D//q.denominator)
    lo=F(t,D); hi=lo if lo*lo==q else F(t+1,D)
    assert lo*lo<=q<=hi*hi
    return lo,hi


def union_area(rects):
    xs=sorted(set(x for r in rects for x in r[:2])); ans=F(0)
    for l,r in zip(xs,xs[1:]):
        mid=(l+r)/2
        ys=sorted((q[2],q[3]) for q in rects if q[0]<mid<q[1])
        if not ys: continue
        a,b=ys[0]; length=F(0)
        for c,d in ys[1:]:
            if c>b: length+=b-a; a,b=c,d
            elif d>b: b=d
        ans+=(r-l)*(length+b-a)
    return ans


def row_area_bounds(N):
    inner=[]; outer=[]
    for p in range(1,N+1):
        if (p-1)**2>=4*p: continue
        lo,hi=root_bracket(F(p))
        for j in range(N-p+1):
            inner.append((F(j+p-1)-lo,F(j)+lo,-lo,lo))
            outer.append((F(j+p-1)-hi,F(j)+hi,-hi,hi))
    return N*union_area(inner),N*union_area(outer)


def main():
    assert os.path.ismount('/Volumes/data')
    src=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11')
    templates=[]; sig=[]
    for name in ['r035_final_0134','r041_final_0036']:
        p=src/'computation/inputs'/f'{name}.json'; data=json.loads(p.read_text())
        px=F(data['Px']);py=F(data['Py']);lam=F(data['lam'])
        pts=[tuple(map(F,q)) for q in data['pts']]; mass=sum(q[2] for q in pts)
        assert len(pts)==204 and lam==1
        assert all(0<=x<px and 0<=y<py and w>0 for x,y,w in pts)
        rho=mass/(px*py); old=rho*(1+px/20)*(1+py/20)
        reflected=rho*(1+px/15)*(1+py/30)
        assert old<1 and reflected<1 and 2*rho>1
        Lx,Ly=ceil(10/px),ceil(10/py)
        axis=max(ceil(20/px)+1,ceil(20/py)+1)
        depth=ceil(30/px); seamLy=ceil(15/py)
        threshold=max(axis,depth,2*seamLy)
        assert (Lx,Ly,axis,depth,seamLy,threshold)==(2,2,5,5,3,6)
        sig.append((px,py,mass))
        templates.append({'input':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
            'atoms':len(pts),'Px':str(px),'Py':str(py),'mass':str(mass),'density':str(rho),
            'old_R10_tail_bound':str(old),'reflected_R15_halfplane_bound':str(reflected),
            'full_double_periodic_density':str(2*rho),
            'old_cell_padding':[Lx,Ly],'old_area_polynomial_N_threshold':2*max(Lx,Ly),
            'axis_N_threshold':axis,'seam_horizontal_source_cells':depth,
            'seam_vertical_padding':seamLy,'reflected_and_increment_N_threshold':threshold})
    assert sig[0]==sig[1]
    # Alternating columns: the real period is 20, but truncation uses steps 10.
    parity=[]
    for N in range(1,11):
        weights=[1 if j%2==0 else 4 for j in range(N)]
        area=N*sum(weights)
        assert area==F(5,2)*N*N-F(3,2)*N*(N%2)
        w=weights[-1]
        parity.append({'N':N,'full_area':str(area),'axis_offset_from_last_source':('1/2' if w==1 else '1'),
                       'excess_coefficient_of_3_minus_2sqrt2':str(N*w)})
    parity_diff=F(parity[3]['full_area'])-3*F(parity[2]['full_area'])+3*F(parity[1]['full_area'])-F(parity[0]['full_area'])
    assert parity_diff==15
    assert F(5,200)*(1+F(20,4))*(1+F(10,4))<1 # R=2 cap.
    assert F(5,200)*(1+F(20,2))*(1+F(10,4))<1 # reflected R=2 cap.
    # Fixed full-cell phase, but the first three sizes precede all local patterns.
    fixed=[]
    assert F(4,20)*(1+F(1,20))*(1+F(20,20))<1 # R=10 cap.
    for N in range(1,9):
        lo,hi=row_area_bounds(N)
        fixed.append({'N':N,'area_lower':str(lo),'area_upper':str(hi)})
    lo4,hi4=map(F,(fixed[3]['area_lower'],fixed[3]['area_upper']))
    lo3,hi3=map(F,(fixed[2]['area_lower'],fixed[2]['area_upper']))
    lo2,hi2=map(F,(fixed[1]['area_lower'],fixed[1]['area_upper']))
    lo1,hi1=map(F,(fixed[0]['area_lower'],fixed[0]['area_upper']))
    dlo=lo4-3*hi3+3*lo2-hi1; dhi=hi4-3*lo3+3*hi2-lo1
    assert F(7,10)<dlo<=dhi<F(3,4)
    # Axis drift under fixed full-cell truncation: mass64, Px1, Py1000.
    axis_tail=F(64,1000)*(1+F(1,100))*(1+F(1000,100))
    assert axis_tail<1
    axis=[{'N':N,'c_formula':f'4*sqrt({N})' if N<=4 else str(N+4),
           'stabilized':N>=4} for N in range(1,9)]
    # Exterior and corner loss for disjoint unit-atom halos.
    crop=[]
    for N in range(1,6):
        full=F(N*N); cut=(N-F(1,2))**2
        assert full-cut==N-F(1,4)
        crop.append({'N':N,'full_area':str(full),'cropped_area':str(cut),'lost_area':str(full-cut)})
    output={'status':'PASS','scope':'Exact tails, safe thresholds, and small counterexamples only; R35/R41 axes and coefficient values were not computed.',
        'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
        'single_process':True,'arithmetic':'Fraction and directed dyadic square roots','sqrt_bits':60,
        'random_seed':None,'floating_tolerance':None,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'templates':templates,'phase_parity_cases':parity,'parity_third_difference_at_N1':str(parity_diff),
        'fixed_phase_small_N_cases':fixed,'fixed_phase_third_difference_interval':[str(dlo),str(dhi)],
        'axis_drift_tail_bound':str(axis_tail),'axis_drift_cases':axis,'exterior_crop_cases':crop}
    path=Path(__file__).parent/'results/finite-state-reflection/exact-check.json'
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','templates':2,'old_area_threshold':4,'axis_threshold':5,
                      'increment_threshold':6,'fixed_phase_fit_failure_sign':1 if dlo>0 else -1,
                      'output':str(path)},ensure_ascii=False))

if __name__=='__main__': main()
