"""Floating local search for D,B3,B5,B7,B9,P11,P19 and one Fourier body.

The hull is formed from boundary samples, therefore this is only a candidate
placement generator for later exact area certification.
"""
from __future__ import annotations
import hashlib, json, math, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

CERT7 = Path("/Volumes/data/CodexProjects/T2/.local/workspaces/lebesgue-main/lebesgue-universal-cover/computation/010-seven-body-certificate.json")
OUT = Path(__file__).with_name("eight_body_numerical_result.json")
SEED = 20261117
DEADLINE_SECONDS = 155.0
ARC_SAMPLES = 48
DISK_SAMPLES = 384
FOURIER_SAMPLES = 512
A3 = 1.0 / 80.0
B5 = -19.0 / 1200.0

class TimeLimit(Exception): pass
def unit(t): return np.column_stack((np.cos(t),np.sin(t)))
def reuleaux(n,rho,tx,ty):
    R=1/(2*math.cos(math.pi/(2*n))); ans=[]
    for ell in range(n):
        phi=rho+2*math.pi*ell/n; c=np.array((tx+R*math.cos(phi),ty+R*math.sin(phi)))
        a=np.linspace(phi+math.pi-math.pi/(2*n),phi+math.pi+math.pi/(2*n),ARC_SAMPLES+1)
        ans.append(c+unit(a))
    return np.vstack(ans)
def polygon(n,rho,tx,ty):
    R=1/(2*math.cos(math.pi/(2*n))); return np.array((tx,ty))+R*unit(rho+2*math.pi*np.arange(n)/n)
def fourier(rho,tx,ty):
    t=2*math.pi*np.arange(FOURIER_SAMPLES)/FOURIER_SAMPLES
    h=.5+A3*np.cos(3*t)+B5*np.cos(5*t)
    hp=-3*A3*np.sin(3*t)-5*B5*np.sin(5*t)
    p=np.column_stack((h*np.cos(t)-hp*np.sin(t),h*np.sin(t)+hp*np.cos(t)))
    c,s=math.cos(rho),math.sin(rho)
    return p@np.array(((c,s),(-s,c)))+np.array((tx,ty))
def points(x):
    assert x.shape==(20,)
    out=[.5*unit(np.linspace(0,2*math.pi,DISK_SAMPLES,endpoint=False)),reuleaux(3,0,x[0],x[1])]
    k=2
    for n in (5,7,9):out.append(reuleaux(n,x[k],x[k+1],x[k+2]));k+=3
    out += [polygon(11,x[11],x[12],x[13]),polygon(19,x[14],x[15],x[16]),fourier(x[17],x[18],x[19])]
    return np.vstack(out)
def area(x):return float(ConvexHull(points(x)).volume)

def main():
    global ARC_SAMPLES,DISK_SAMPLES,FOURIER_SAMPLES
    source_bytes=CERT7.read_bytes();source=json.loads(source_bytes)
    base=np.array([float(v) for v in source['parameter_decimals_interpreted_exactly']])
    assert base.shape==(17,)
    started=time.monotonic();calls=0;best=np.r_[base,[0.,0.,0.]];best_area=area(best)
    bounds=[(-.1,.1),(-.1,.1)]
    for n in (5,7,9):bounds += [(0,2*math.pi/n),(-.1,.1),(-.1,.1)]
    bounds += [(0,2*math.pi/11),(-.1,.1),(-.1,.1),(0,2*math.pi/19),(-.1,.1),(-.1,.1),(0,2*math.pi),(-.1,.1),(-.1,.1)]
    assert len(bounds)==20
    def obj(x):
        nonlocal calls,best,best_area
        calls+=1
        if time.monotonic()-started>DEADLINE_SECONDS:raise TimeLimit
        v=area(x)
        if v<best_area:best_area,best=v,np.array(x,copy=True)
        return v
    phases=[]
    # The Fourier body is reflection-symmetric, so one orientation branch is
    # enough.  A small orientation multi-start is useful because sampled-hull
    # area is nonsmooth; unlike bounded Powell, Nelder--Mead retains its start.
    try:
        local=[]
        for k in range(12):
            start=np.array((2*math.pi*k/12,0.,0.))
            sol=minimize(lambda z:obj(np.r_[base,z]),start,method='Nelder-Mead',options={'maxiter':350,'xatol':2e-8,'fatol':2e-10})
            local.append({'start_angle':float(start[0]),'success':bool(sol.success),'message':str(sol.message),'fun':float(sol.fun),'nit':int(sol.nit),'nfev':int(sol.nfev)})
        phases.append({'phase':'fourier-only-multistart','starts':local})
        seed=best
    except TimeLimit:
        phases.append({'phase':'fourier-only','stopped':'deadline'});seed=best
    rng=np.random.default_rng(SEED);lo=np.array([z[0] for z in bounds]);hi=np.array([z[1] for z in bounds])
    for i,start in enumerate([seed,np.clip(seed+rng.normal(0,.0015,20),lo,hi)]):
        try:
            sol=minimize(obj,start,method='Powell',options={'maxiter':1050,'xtol':2e-8,'ftol':2e-10})
            phases.append({'phase':f'joint-{i}','success':bool(sol.success),'message':str(sol.message),'fun':float(sol.fun),'nit':int(sol.nit),'nfev':int(sol.nfev)})
        except TimeLimit:
            phases.append({'phase':f'joint-{i}','stopped':'deadline'});break
    coarse=area(best);ARC_SAMPLES,DISK_SAMPLES,FOURIER_SAMPLES=192,1536,2048;fine=area(best)
    result={'status':'completed','scope':'floating local search only; sampled-boundary hull and Powell do not certify m(F8), a lower bound, or global placement coverage','seed':SEED,'elapsed_seconds':time.monotonic()-started,'objective_calls':calls,'parameterization':'disk fixed at origin; B3 angle fixed at zero; x=[t3x,t3y,rho5,t5x,t5y,rho7,t7x,t7y,rho9,t9x,t9y,rho11,t11x,t11y,rho19,t19x,t19y,rhoF,tFx,tFy]','degrees_of_freedom':20,'input_seven_certificate':str(CERT7),'input_seven_certificate_sha256':hashlib.sha256(source_bytes).hexdigest(),'input_seven_parameters_exact_strings':source['parameter_decimals_interpreted_exactly'],'fourier_shape':'h(t)=1/2+cos(3t)/80-19*cos(5t)/1200','curvature_radius_lower_bound':.02,'best_parameters':best.tolist(),'coarse_sample_area':coarse,'fine_sample_area':fine,'coarse_sampling':{'arc_segments_per_Reuleaux_arc':48,'disk_vertices':384,'fourier_boundary_points':512},'fine_sampling':{'arc_segments_per_Reuleaux_arc':192,'disk_vertices':1536,'fourier_boundary_points':2048},'phases':phases,'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
