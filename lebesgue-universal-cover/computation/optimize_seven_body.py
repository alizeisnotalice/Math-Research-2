"""Floating local search for hull(D,B3,B5,B7,B9,P11,P19).

The objective is a sampled-boundary convex hull.  This is local numerical
exploration, not a certified value for the seven-body test family.
"""
from __future__ import annotations
import hashlib, json, math, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

OUT = Path(__file__).with_name("seven_body_numerical_result.json")
SEED = 20261117
DEADLINE_SECONDS = 155.0
ARC_SAMPLES = 48
DISK_SAMPLES = 384
# Exact decimal parameters used by the round-7 certificate, converted here to
# binary floats only for this numerical experiment.
X_SIX = np.array([
 -0.002116717736623155,-0.002172190513432106,
 0.03369871051209932,0.025092651107061735,-0.0038921881207567045,
 0.7087849854618435,0.0035628282693877652,-0.01262876252254122,
 0.34909300190743575,-0.0006339143657618253,-0.0006497687888177878,
 0.11685175287929492,-0.0009822210255093589,0.006365090528266086,
])
# Earlier P19 relaxation is only a numerical starting guess.
P19_START = np.array([0.2039279441803901,0.00015808203494954365,-0.000864531481407799])
class TimeLimit(Exception): pass

def unit(a): return np.column_stack((np.cos(a),np.sin(a)))
def reuleaux(n,rho,tx,ty):
    R=1/(2*math.cos(math.pi/(2*n))); groups=[]
    for ell in range(n):
        phi=rho+2*math.pi*ell/n; c=np.array((tx+R*math.cos(phi),ty+R*math.sin(phi)))
        ang=np.linspace(phi+math.pi-math.pi/(2*n),phi+math.pi+math.pi/(2*n),ARC_SAMPLES+1)
        groups.append(c+unit(ang))
    return np.vstack(groups)
def polygon(n,rho,tx,ty):
    R=1/(2*math.cos(math.pi/(2*n))); return np.array((tx,ty))+R*unit(rho+2*math.pi*np.arange(n)/n)
def points(x):
    assert x.shape==(17,)
    z=[.5*unit(np.linspace(0,2*math.pi,DISK_SAMPLES,endpoint=False)),reuleaux(3,0,x[0],x[1])]
    k=2
    for n in (5,7,9): z.append(reuleaux(n,x[k],x[k+1],x[k+2]));k+=3
    z += [polygon(11,x[11],x[12],x[13]),polygon(19,x[14],x[15],x[16])]
    return np.vstack(z)
def area(x): return float(ConvexHull(points(x)).volume)

def main():
    global ARC_SAMPLES,DISK_SAMPLES
    began=time.monotonic(); calls=0; best=np.r_[X_SIX,P19_START]; best_area=area(best)
    bounds=[(-.1,.1),(-.1,.1)]
    for n in (5,7,9): bounds += [(0,2*math.pi/n),(-.1,.1),(-.1,.1)]
    bounds += [(0,2*math.pi/11),(-.1,.1),(-.1,.1),(0,2*math.pi/19),(-.1,.1),(-.1,.1)]
    assert len(bounds)==17
    def obj(x):
        nonlocal calls,best,best_area
        calls+=1
        if time.monotonic()-began>DEADLINE_SECONDS: raise TimeLimit
        v=area(x)
        if v<best_area: best_area,best=v,np.array(x,copy=True)
        return v
    phases=[]
    try:
        sol=minimize(lambda q:obj(np.r_[X_SIX,q]),P19_START,method='Powell',bounds=bounds[14:],options={'maxiter':480,'xtol':2e-9,'ftol':2e-11})
        phases.append({'phase':'P19-only','success':bool(sol.success),'message':str(sol.message),'fun':float(sol.fun),'nit':int(sol.nit),'nfev':int(sol.nfev)})
        first=np.r_[X_SIX,sol.x]
    except TimeLimit:
        phases.append({'phase':'P19-only','stopped':'deadline'});first=best
    rng=np.random.default_rng(SEED);lo=np.array([b[0] for b in bounds]);hi=np.array([b[1] for b in bounds])
    for i,start in enumerate([first,np.clip(first+rng.normal(0,.0015,17),lo,hi)]):
        try:
            sol=minimize(obj,start,method='Powell',bounds=bounds,options={'maxiter':1000,'xtol':2e-8,'ftol':2e-10})
            phases.append({'phase':f'joint-{i}','success':bool(sol.success),'message':str(sol.message),'fun':float(sol.fun),'nit':int(sol.nit),'nfev':int(sol.nfev)})
        except TimeLimit:
            phases.append({'phase':f'joint-{i}','stopped':'deadline'});break
    coarse=area(best);ARC_SAMPLES,DISK_SAMPLES=192,1536;fine=area(best)
    result={'status':'completed','scope':'floating local search only; sampled boundary hull and Powell do not certify m(F7), a lower bound, or global placement coverage','seed':SEED,'elapsed_seconds':time.monotonic()-began,'objective_calls':calls,'parameterization':'disk fixed at origin; B3 angle fixed at zero; x=[t3x,t3y,rho5,t5x,t5y,rho7,t7x,t7y,rho9,t9x,t9y,rho11,t11x,t11y,rho19,t19x,t19y]','degrees_of_freedom':17,'input_six_parameters':X_SIX.tolist(),'best_parameters':best.tolist(),'coarse_sample_area':coarse,'fine_sample_area':fine,'coarse_sampling':{'arc_segments_per_Reuleaux_arc':48,'disk_vertices':384},'fine_sampling':{'arc_segments_per_Reuleaux_arc':192,'disk_vertices':1536},'phases':phases,'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__': main()
