#!/usr/bin/env python3
"""Exact small toys for R54 double-reflection closure; no template N scan."""
from fractions import Fraction as F
from math import ceil,floor,isqrt
from pathlib import Path
import hashlib,json,os,sys


def root(q,bits=60):
    D=1<<bits;t=isqrt(q.numerator*D*D//q.denominator)
    lo=F(t,D);hi=lo if lo*lo==q else F(t+1,D)
    assert lo*lo<=q<=hi*hi
    return lo,hi


def union_area(rects):
    xs=sorted({x for q in rects for x in q[:2]});ans=F(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2;ys=sorted((q[2],q[3]) for q in rects if q[0]<mid<q[1])
        if not ys:continue
        lo,hi=ys[0];dy=F(0)
        for c,d in ys[1:]:
            if c>hi:dy+=hi-lo;lo,hi=c,d
            elif d>hi:hi=d
        ans+=(b-a)*(dy+hi-lo)
    return ans


def halo_area(xs,win):
    # All toys have a certified cap2 and all horizontal gaps are 1,9,10.
    # Hence at most a singleton or one adjacent gap1 pair can be captured;
    # distinct vertical rows are distance10 and cannot interact.
    inner=[];outer=[];x0,x1,y0,y1=win
    for y in [F(0),F(10)]:
        for i,a in enumerate(xs):
            boxes=[(a,a,F(1))]
            if i+1<len(xs) and xs[i+1]-a<4:boxes.append((a,xs[i+1],F(2)))
            for a,b,m in boxes:
                if (b-a)**2>=m:continue
                rlo,rhi=root(m/4)
                for r,store in [(rlo,inner),(rhi,outer)]:
                    q=(max(x0,b-r),min(x1,a+r),max(y0,y-r),min(y1,y+r))
                    if q[0]<q[1] and q[2]<q[3]:store.append(q)
    return union_area(inner),union_area(outer)


def toy(alphas):
    P=F(10);Py=F(10);N=7;m=len(alphas)
    old=sorted(F(a)+i*P for i in range(N) for a in alphas)
    ell=old[0]-F(1,2);right=old[-1]+F(1,2);W=right-ell;period=2*W
    assert W>(N-1)*P and W>6*P and W>30
    base=sorted(old+[2*right-x for x in old])
    assert len(base)==2*N*m
    x0,x1=ell,ell+period
    xs=sorted({x+k*period for k in [-1,0,1] for x in base if x0-2<=x+k*period<=x1+2})
    full=halo_area(xs,(x0,x1,F(0),Py));old_area=halo_area(old,(ell,right,F(0),Py))
    slo,shi=root(F(2),80);klo=3-2*shi;khi=3-2*slo
    if len(alphas)==1:
        target=(2*N+2*klo,2*N+2*khi);old_target=(F(N),F(N));bx='0'
    else:
        target=(4*N+2*N*klo,4*N+2*N*khi)
        old_target=(2*N+(N-1)*klo,2*N+(N-1)*khi);bx='-(3-2*sqrt(2))'
    assert full[0]<=target[0]<=target[1]<=full[1]
    assert old_area[0]<=old_target[0]<=old_target[1]<=old_area[1]
    assert full[0]-2*old_area[1]<=2*klo<=2*khi<=full[1]-2*old_area[0]
    rho=F(m,100)
    assert rho*(1+P/2)*(1+Py/4)<1 # short-radius cap2.
    assert rho*F(3,2)*(1+Py/(6*P))<1 # all long radii.
    # Per-label closed-interval counts, bounded deterministic small test.
    tests=0
    for alpha in alphas:
        bb=[F(alpha)+i*P for i in range(N)]
        bb+= [2*right-(F(alpha)+i*P) for i in range(N)]
        starts=[ell-W,ell,F(0),right-F(1,2),right+1,ell+period-F(1,2)]
        lengths=[F(1),F(2),F(3),F(15),F(30),W-F(1,2),W,W+F(1,2),2*W,3*W+F(1,2)]
        for a in starts:
            for L in lengths:
                z={x+k*period for x in bb for k in range(floor((a-max(bb))/period)-1,ceil((a+L-min(bb))/period)+2)}
                count=sum(a<=x<=a+L for x in z)
                bound=L/P+2 if L<W else L/P+L/W+2
                assert count<=bound;tests+=1
    return {'alphas':list(map(str,alphas)),'N':N,'mass_per_base_cell':m,'W':str(W),
            'double_period':str(period),'double_mass':2*N*m,
            'full_area_interval':list(map(str,full)),'old_strip_area_interval':list(map(str,old_area)),
            'bx':bx,'gamma_left':'3/2-sqrt(2)','gamma_right':'3/2-sqrt(2)',
            'coefficient_expected':'3-2*sqrt(2)' if m==1 else '0',
            'count_tests':tests,'full_area_source_points':2*len(xs)}


def main():
    assert os.path.ismount('/Volumes/data')
    src=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json')
    d=json.loads(src.read_text());px,py=F(d['Px']),F(d['Py']);mass=sum(F(p[2]) for p in d['pts']);rho=mass/(px*py)
    short=rho*(1+px/15)*(1+py/30)
    long=rho*F(3,2)*(1+py/(6*px))
    assert short<1 and long<1 and 6*px>30
    toys=[toy([F(0)]),toy([F(0),F(9)])]
    # Naive canonical x-flip of left source with alpha=0,9 changes the phase:
    # true flipped rightmost source is0; naive uniform i<=0 adds alpha'=1.
    assert F(2,100)*(1+F(10,4))**2<1
    def max_right(points):
        points=sorted(points);lower=[];upper=[]
        for i,x in enumerate(points):
            for j in range(i,len(points)):
                if points[j]-x>=4:break
                mass=F(j-i+1)
                if (points[j]-x)**2>=mass:continue
                lo,hi=root(mass/4);lower.append(x+lo);upper.append(x+hi)
        return max(lower),max(upper)
    true_flip={-a-10*i for a in [F(0),F(9)] for i in range(4)}
    naive_flip={a+10*i for a in [F(0),F(1)] for i in range(-3,1)}
    assert max_right(true_flip)==(F(1,2),F(1,2))
    assert max_right(naive_flip)==(F(3,2),F(3,2))
    phase={'original_alphas':['0','9'],'true_flipped_cutoffs':{'0':'i<=0','1':'i<=-1'},
           'naive_flipped_cutoff':'both i<=0','true_right_halo_edge':'1/2','naive_right_halo_edge':'3/2',
           'actual_directed_root_bbox_check':'PASS'}
    ops=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/reflection')
    bulk=json.loads((ops/'r041-bulk.json').read_text())
    left=json.loads((ops/'r041-left.json').read_text());right=json.loads((ops/'r041-right.json').read_text())
    prior=Path(__file__).parent/'results/reflection-edge-state-review/r041-independent-coefficients.json'
    certificate=json.loads(prior.read_text())
    bx=(F(left['inner_area'])+F(right['inner_area'])-4*F(bulk['outer_area']),
        F(left['outer_area'])+F(right['outer_area'])-4*F(bulk['inner_area']))
    gr=tuple(map(F,certificate['gamma_rational_axis_interval']));eg=F(certificate['gamma_axis_error'])
    partial=(bx[0]+gr[0]-eg,bx[1]+gr[1]+eg)
    existing={'template':'R41','bx_interval':list(map(str,bx)),
              'bx_plus_true_gammaR_interval':list(map(str,partial)),
              'gammaL_break_even_interval':list(map(str,(-partial[1],-partial[0]))),
              'source_hashes':{str(prior):hashlib.sha256(prior.read_bytes()).hexdigest(),
                 **{str(ops/f'r041-{name}.json'):hashlib.sha256((ops/f'r041-{name}.json').read_bytes()).hexdigest()
                    for name in ['bulk','left','right']}}}
    out={'status':'PASS','scope':'Exact count constants and two N=7 toys only, no real-template N scan or actual left-seam area run.',
         'python':sys.version,'single_process':True,'arithmetic':'Fraction and directed dyadic roots',
         'sqrt_bits':60,'reference_sqrt_bits':80,'random_seed':None,'floating_tolerance':None,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
         'short_R15_tail_bound':str(short),'long_radius_tail_bound':str(long),
         'long_radius_tail_bound_float_display':float(long),'six_Px_minus_30':str(6*px-30),
         'toys':toys,'phase_counterexample':phase,'existing_R41_partial_coefficient':existing}
    path=Path(__file__).parent/'results/double-reflection-closure/exact-check.json'
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','toys':2,'N':7,'interval_count_tests':sum(t['count_tests'] for t in toys),
                      'long_tail_bound_display':float(long),'output':str(path)}))

if __name__=='__main__':main()
