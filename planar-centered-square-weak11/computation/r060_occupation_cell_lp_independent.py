#!/usr/bin/env python3
"""Read-only R60 occupation formula/geometry audit; no MATLAB or LP."""
from pathlib import Path
from fractions import Fraction as F
import datetime,hashlib,json,math,os,sys,time
import mpmath as mp
import numpy as np
BASE=Path(__file__).resolve().parent.parent
SOURCE=BASE/'computation/cauchy_occupation_cell_lp.m'
if not SOURCE.is_file():SOURCE=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-kernel/planar-centered-square-weak11/computation/cauchy_occupation_cell_lp.m')
PROOF=BASE/'proofs/r060-positive-cell-occupation-interface-independent.md'
if not PROOF.is_file():PROOF=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/proofs/r060-positive-cell-occupation-interface-independent.md')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def orbits(N,delta):
    out=[]
    for i in range(N+1):
        for j in range(i+1):
            x,y=delta*i,delta*j
            pts=set()
            for sx in [-1,1]:
                for sy in [-1,1]:pts.update([(sx*x,sy*y),(sx*y,sy*x)])
            out.append(((x,y),sorted(pts)))
    return out
def prob(d,delta,h):return mp.atan2(h*delta,h*h+d*d-delta*delta/4)/mp.pi
def pp(x,y,h):return h*h/(mp.pi**2*(h*h+x*x)*(h*h+y*y))
def mpF(v):return mp.mpf(v.numerator)/v.denominator

SCALE=1<<192
def updiv(a,b):return -((-a)//b)
def floor_fixed(x):return x.numerator*SCALE//x.denominator
def ceil_fixed(x):return updiv(x.numerator*SCALE,x.denominator)
def imul(a,b):return (a[0]*b[0]//SCALE,updiv(a[1]*b[1],SCALE))
def atan_small(x):
    assert 0<=x<=F(1,2)
    power=(floor_fixed(x),ceil_fixed(x));sq=imul(power,power);lo=hi=0
    for n in range(80):
        k=2*n+1;tl=power[0]//k;th=updiv(power[1],k)
        if n%2:lo-=th;hi-=tl
        else:lo+=tl;hi+=th
        power=imul(power,sq)
    rem=updiv(power[1],161)
    return lo-rem,hi+rem
A5=atan_small(F(1,5));A239=atan_small(F(1,239))
PILO=16*A5[0]-4*A239[1];PIHI=16*A5[1]-4*A239[0]
# Exact tangent verification of the Machin angle identity.
tan2=2*F(1,5)/(1-F(1,5)**2);tan4=2*tan2/(1-tan2**2)
assert (tan4-F(1,239))/(1+tan4*F(1,239))==1

def atan_pos(x):
    assert x>=0
    if x>1:
        lo,hi=atan_pos(1/x);return PILO//2-hi,updiv(PIHI,2)-lo
    if x>F(1,2):
        lo,hi=atan_small((1-x)/(1+x));return PILO//4-hi,updiv(PIHI,4)-lo
    return atan_small(x)
def prob_upper_fixed(dist,delta,h):
    Y=h*delta;X=h*h+dist*dist-delta*delta/4
    if X>0:lo,hi=atan_pos(Y/X)
    elif X==0:lo,hi=PILO//2,updiv(PIHI,2)
    else:
        al,ah=atan_pos(Y/(-X));lo,hi=PILO-ah,PIHI-al
    assert hi>0
    return updiv(hi*SCALE,PILO)
def completed_results():
    inp=BASE/'computation/inputs/r060-occupation-cell-lp-independent/completed-cases.json';data=json.loads(inp.read_text());rows=[]
    for c in data['cases']:
        hh,dd=map(F,(c['h'],c['delta']));lam=F(c['source_origin_mass']);aa=list(map(F,c['density_coefficients']));assert lam>=0 and all(a>=0 for a in aa);oo=orbits(c['N'],dd)
        assert len(oo)==len(aa)==len(c['orbit_sizes'])
        assert all(len(pts)==n for (_,pts),n in zip(oo,c['orbit_sizes']))
        mu=lam+dd*dd*sum(a*len(pts)for a,(_,pts)in zip(aa,oo));cost=mp.pi*mpF(hh)*mpF(mu)
        assert abs(cost-c['cost_from_exact_mass_formula_float'])<mp.mpf('1e-12')
        x,y=map(F,c['source_residual_sample_point']);target=c['source_residual_cell_index']-1
        residual=mpF(lam)*pp(mpF(x),mpF(y),mpF(hh))-mpF(aa[target])
        for a,(_,pts)in zip(aa,oo):
            if not a:continue
            residual+=mpF(a)*sum(prob(abs(mpF(x-u)),mpF(dd),mpF(hh))*prob(abs(mpF(y-v)),mpF(dd),mpF(hh))for u,v in pts)
        assert abs(residual-c['source_residual_sample_min'])<mp.mpf('1e-10')
        rows.append(dict(h=c['h'],delta=c['delta'],N=c['N'],mode=c['mode'],source=c['source'],mu_mass_exact=str(mu),halfmass_high_precision=mp.nstr(cost,45),halfmass_printed=c['cost_from_exact_mass_formula_float'],source_constraint_min_float=c['source_residual_constraint_min_float'],source_sample_min_saved=c['source_residual_sample_min'],source_witness_reconstructed=mp.nstr(residual,40),source_sample_point=c['source_residual_sample_point']))
    best=min((c for c in data['cases'] if c['mode']=='midpoint'),key=lambda c:c['cost_from_exact_mass_formula_float'])
    hh,dd=map(F,(best['h'],best['delta']));lam=F(best['source_origin_mass']);aa=list(map(F,best['density_coefficients']));oo=orbits(best['N'],dd);target=best['source_residual_cell_index']-1
    rep=oo[target][0];corner=tuple(map(F,best['source_residual_sample_point']));eps=F(1,10**8)
    x,y=tuple((1-eps)*a+eps*b for a,b in zip(corner,rep));margin=min(dd/2-abs(x-rep[0]),dd/2-abs(y-rep[1]));assert margin>0
    pi_lower=F(PILO,SCALE)
    ph_upper=lam*hh*hh/(pi_lower*pi_lower*(hh*hh+x*x)*(hh*hh+y*y));upper=ceil_fixed(ph_upper)
    for a,(_,pts)in zip(aa,oo):
        if not a:continue
        for u,v in pts:
            px=prob_upper_fixed(abs(x-u),dd,hh);py=prob_upper_fixed(abs(y-v),dd,hh)
            upper+=ceil_fixed(a*F(px*py,SCALE*SCALE))
    residual_upper=F(upper,SCALE)-aa[target];assert residual_upper<-1
    mu=lam+dd*dd*sum(a*len(pts)for a,(_,pts)in zip(aa,oo));lip=2*mu/(pi_lower*pi_lower*hh**3)
    radius=min(margin/2,-residual_upper/(2*lip));assert radius>0
    cert=dict(status='EXACT interval upper bound strictly negative at interior point and a positive-area neighborhood',h_exact=str(hh),delta_exact=str(dd),case_source=best['source'],target_orbit_index_1based=target+1,target_center=list(map(str,rep)),interior_point=list(map(str,(x,y))),distance_to_target_cell_boundary=str(margin),lambda_exact=str(lam),target_density_exact=str(aa[target]),source_residual_upper_exact=str(residual_upper),source_residual_upper_numeric=float(residual_upper),rho_neighborhood_radius_exact=str(radius),rho_neighborhood_radius_numeric=float(radius),pi_lower_exact=str(pi_lower),pi_upper_exact=str(F(PIHI,SCALE)),fixedpoint_bits=192,atan_alternating_terms=80,rounding='Each multiplication/division outward-rounded in exact integer arithmetic; final weighted terms rounded upward on common2^192 lattice.',neighborhood_argument='Within the targetcell f is constant; each partial derivative of lambda*p_h+P_h*f is at most2*mu_mass/(pi_lower^2*h^3). The stated rho-ball remains insidecell and keeps source residual negative.')
    return dict(input_path=str(inp),input_sha256=sha(inp),receipt=data['receipt'],case_count=len(rows),cases=rows,rigorous_midpoint_invalidity=cert)

def main():
    assert os.path.ismount('/Volumes/data');start=time.monotonic();mp.mp.dps=80
    assert sha(SOURCE)=='dab5992d41449cfeb047455d472d106ee05ce691e82f8f41f7e65bc241c9c1fd'
    geometry=[]
    for N in [1,2,4]:
        delta=F(1,8);oo=orbits(N,delta);flat=[p for _,points in oo for p in points]
        assert len(set(flat))==len(flat)==(2*N+1)**2
        assert len(oo)==(N+1)*(N+2)//2
        assert set(flat)=={(i*delta,j*delta)for i in range(-N,N+1)for j in range(-N,N+1)}
        geometry.append(dict(N=N,delta=str(delta),orbits=len(oo),physical_cells=len(flat),exact_tiling=True))
    delta=mp.mpf('0.125');h=mp.mpf('0.3');oo=orbits(2,F(1,8));sizes=[len(p) for _,p in oo]
    matrices={};res=[]
    for mode in ['midpoint','cellinf']:
        L=[]
        for (x,y),_ in oo:
            x,y=mpF(x),mpF(y);row=[]
            for _,pts in oo:
                val=mp.mpf(0)
                for u,v in pts:
                    dx=abs(x-mpF(u));dy=abs(y-mpF(v))
                    if mode=='cellinf':dx+=delta/2;dy+=delta/2
                    val+=prob(dx,delta,h)*prob(dy,delta,h)
                row.append(val)
            L.append(row)
        defect=max(abs(sizes[i]*L[i][j]-sizes[j]*L[j][i]) for i in range(len(oo))for j in range(len(oo)))
        assert defect<mp.mpf('1e-70')
        if mode=='midpoint':
            B=(mp.mpf(2)+mp.mpf('.5'))*delta
            for i,((x,y),_)in enumerate(oo):
                x,y=mpF(x),mpF(y)
                direct=(mp.atan((B-x)/h)-mp.atan((-B-x)/h))*(mp.atan((B-y)/h)-mp.atan((-B-y)/h))/mp.pi**2
                assert abs(sum(L[i])-direct)<mp.mpf('1e-70')
        matrices[mode]=L
        res.append(dict(mode=mode,detailed_balance_max_defect=mp.nstr(defect,8),row_sums=[mp.nstr(sum(row),30)for row in L]))
    assert all(matrices['cellinf'][i][j]<=matrices['midpoint'][i][j]for i in range(len(oo))for j in range(len(oo)))
    lower_tests=[]
    # Independently check targetcell corner/center values against each entrywise lower sourcecell term.
    for (x,y),_ in oo:
        x,y=mpF(x),mpF(y)
        for u,v in [oo[0][1][0],oo[-1][1][-1]]:
            u,v=mpF(u),mpF(v)
            low=prob(abs(x-u)+delta/2,delta,h)*prob(abs(y-v)+delta/2,delta,h)
            vals=[prob(abs(x+sx*delta/2-u),delta,h)*prob(abs(y+sy*delta/2-v),delta,h) for sx,sy in [(0,0),(-1,-1),(-1,1),(1,-1),(1,1)]]
            assert all(low<=val+mp.mpf('1e-70')for val in vals)
            lower_tests.append(dict(target=[mp.nstr(x),mp.nstr(y)],source=[mp.nstr(u),mp.nstr(v)],lower=mp.nstr(low,30),min_sample=mp.nstr(min(vals),30)))
    node,weight=np.polynomial.legendre.leggauss(24);node=(node+1)/2;weight=weight/2
    qtests=[]
    for delta,h,dx,dy in [('0.125','.3','0','0'),('.0625','.6','0','0'),('.0625','.6','.03125','.03125'),('.0625','.6','.03125001','.03125001'),('.125','.15','.25','.125')]:
        de,hh,xx,yy=map(mp.mpf,(delta,h,dx,dy))
        integrand=lambda s:prob(xx,de,s)*prob(yy,de,s)
        pts={mp.mpf(0),hh}
        for dist in [abs(xx-de/2),abs(xx+de/2),abs(yy-de/2),abs(yy+de/2)]:
            if 0<dist<hh:pts.add(dist)
        exact=2*mp.pi*mp.quad(integrand,sorted(pts))
        def pf(d,s):return np.arctan2(s*float(de),s*s+float(d*d-de*de/4))/math.pi
        ss=float(hh)*node
        GL=2*math.pi*float(hh)*np.dot(weight,pf(xx,ss)*pf(yy,ss))
        qtests.append(dict(delta=delta,h=h,dx=dx,dy=dy,reference=mp.nstr(exact,50),GL24_numeric=GL,relative_error=float((mp.mpf(GL)-exact)/exact)))
    completed=completed_results()
    out=dict(completed_case_audit=completed,status='PASS bounded occupation geometry/matrix/normalization checks',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,executable=sys.executable,mpmath=mp.__version__,numpy=np.__version__,precision_digits=80,workers=1,script_sha256=sha(Path(__file__)),sources=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)for p in [SOURCE,PROOF]],geometry=geometry,matrix_checks=res,entrywise_lower_samples=lower_tests,Qcell_GL24_checks=qtests,normalization='Each a_j is density per physical cell; halfmass pi*h*(lambda+delta^2*sum_j orbit_size_j*a_j). Matrix rows target, columns sourceorbit.',scope='No MATLAB/LP launched. Twelve saved results reconstructed. Exact directed certificate proves best-midpoint C4 invalidity only; quadrature and other result reconstructions remain nondirected. Source-proof density/C4 analyzed separately.',elapsed_seconds=time.monotonic()-start)
    p=BASE/'computation/results/r060-occupation-cell-lp-independent/checks.json';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],output=str(p),output_sha256=sha(p),script_sha256=out['script_sha256'],geometry=geometry,GL24_relative_errors=[r['relative_error']for r in qtests],sources=out['sources'],completed_case_count=completed['case_count'],rigorous_midpoint_upper=completed['rigorous_midpoint_invalidity']['source_residual_upper_numeric'],rigorous_midpoint_neighborhood_radius=completed['rigorous_midpoint_invalidity']['rho_neighborhood_radius_numeric'],elapsed_seconds=out['elapsed_seconds']),indent=2))
if __name__=='__main__':main()
