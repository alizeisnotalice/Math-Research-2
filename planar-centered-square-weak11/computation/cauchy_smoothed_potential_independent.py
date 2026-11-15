#!/usr/bin/env python3
"""R55 independent scalar check; nondirected arithmetic is not an enclosure."""
from pathlib import Path
import datetime,hashlib,json,math,os,platform,sys,time
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp
import sympy as sp
import numpy as np
BASE=Path(__file__).resolve().parent.parent
SOURCE_CANDIDATES=[
    ('proofs/cauchy-tail-mass-obstruction-and-heavy-interface-independent.md',Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/proofs/cauchy-tail-mass-obstruction-and-heavy-interface-independent.md')),
    ('proofs/cauchy-potential-positive-part-independent.md',Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/proofs/cauchy-potential-positive-part-independent.md')),
    ('computation/cauchy_heavy_tail_fixed_screen.m',Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-kernel/planar-centered-square-weak11/computation/cauchy_heavy_tail_fixed_screen.m')),
]
# Integrated project files take precedence; fallbacks support the helper clone only.
SOURCES=[BASE/relative if (BASE/relative).is_file() else fallback for relative,fallback in SOURCE_CANDIDATES]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atanc(z):return math.atan(z)/z if z else 1.
def smooth(a,b,t):
    a,b=sorted((abs(a),abs(b)))
    if b==0:return 2/(math.pi*t)
    h=t*t+a*b
    return 2/math.pi*(math.atan2(b,t)+a*t/h*atanc(t*(b-a)/h))/(a+b)
def raw(a,b,t):
    if a==b:
        return 2/(math.pi*t) if a==0 else (math.atan(a/t)/a+t/(t*t+a*a))/math.pi
    return 2/math.pi*(b*math.atan(b/t)-a*math.atan(a/t))/(b*b-a*a)
def q_closed(a,b,t):
    a,b=sorted((abs(a),abs(b)))
    if b==0:return math.inf
    h=t*t+a*b
    return 2/math.pi*(math.atan2(t,b)-a*t/h*atanc(t*(b-a)/h))/(a+b)
def reference(a,b,t):
    a,b,t=map(lambda z:mp.mpf(str(z)),(a,b,t))
    if a==b:
        return 2/(mp.pi*t) if a==0 else (mp.atan(a/t)/a+t/(t*t+a*a))/mp.pi
    return 2/mp.pi*(b*mp.atan(b/t)-a*mp.atan(a/t))/(b*b-a*a)
def rel(x,y):return float(abs(mp.mpf(x)-y)/abs(y))
def orbit_checks():
    inp=BASE/'computation/inputs/cauchy-smoothed-potential-independent/fixed-positive-orbits.json'
    data=json.loads(inp.read_text(),parse_float=F)
    atoms=defaultdict(F)
    for o in data['orbits']:
        assert F(o['coefficient'])>=0
        for x,y in o['points']:atoms[F(x),F(y)]+=F(o['coefficient'])/len(o['points'])
    atoms={p:a for p,a in atoms.items() if a}
    for (x,y),a in atoms.items():
        for p in [(y,x),(-x,y),(x,-y),(-x,-y)]:assert atoms.get(p)==a
    A=sum(atoms.values());assert all(sum(a*p[k] for p,a in atoms.items())==0 for k in (0,1))
    L=max(max(abs(x),abs(y)) for x,y in atoms);S=max(abs(x)+abs(y) for x,y in atoms);assert S<64
    points=[(F(64),F(0)),(F(64),F(1,8)),(F(40),F(40)),(F(10**8),F(10**8)),(F(1,8),F(128))]
    checks=[]
    for x,y in points:
        r=x+y
        ds={p:abs(x-p[0])+abs(y-p[1])-r for p in atoms}
        d1=sum(a*ds[p] for p,a in atoms.items())
        paired=sum(2*a*max(p[0]-x,0) for p,a in atoms.items() if p[0]>0)+sum(2*a*max(p[1]-y,0) for p,a in atoms.items() if p[1]>0)
        assert d1==paired
        direct=A/r-sum(a/(r+ds[p]) for p,a in atoms.items())
        stable=paired/r**2-sum(a*ds[p]**2/(r**2*(r+ds[p])) for p,a in atoms.items())
        assert direct==stable
        checks.append(dict(point=[str(x),str(y)],first_moment=str(d1),Wcompact=str(stable)))
    positive=[dict(point=[str(x),str(y)],mass=str(a)) for (x,y),a in sorted(atoms.items())]
    return dict(path=str(inp),sha256=sha(inp),bytes=inp.stat().st_size,original_source=data['source'],positive_atom_count=len(atoms),total_mass=str(A),L=str(L),S=str(S),D4_exact=True,zero_position_mean_exact=True,nonnegative_exact=True,far_source_guard_exact=True,first_moment_and_compact_checks=checks,minimal_input=positive)
def main():
    assert os.path.ismount('/Volumes/data')
    start=time.monotonic();mp.mp.dps=100
    a,b,t,s,h=sp.symbols('a b t s h',positive=True)
    exact={}
    def check(name,expr):
        val=sp.factor(sp.together(expr));assert val==0,(name,val);exact[name]='exact zero'
    F=(b*sp.atan(b/t)-a*sp.atan(a/t))/(b*b-a*a)
    check('time_derivative',sp.diff(F,t)+t*t/((t*t+a*a)*(t*t+b*b)))
    check('diagonal_limit',sp.limit(F,b,a)-(sp.atan(a/t)/a+t/(t*t+a*a))/2)
    check('axis_limit',sp.limit(F,a,0)-sp.atan(b/t)/b)
    D=3*sp.sqrt(3)/(8*sp.pi)
    check('derivative_maximizer',sp.diff(a/(s*s+a*a)**2,a).subs(a,s/sp.sqrt(3)))
    check('derivative_maximum_value',(a/(s*s+a*a)**2).subs(a,s/sp.sqrt(3))-9/(16*sp.sqrt(3)*s**3))
    check('compact_cancellation',h/(s*(s+h))-(h/s**2-h*h/(s*s*(s+h))))
    tests=[]
    values=[(0,0,1),(0,1,1),(1,1,1),(1,1+1e-14,1),(1e-10,2e-10,1),(1e4,1e4*(1+1e-12),.1),(1e8,1e8+1,.01),(0,1e8,.01),(1,2,1e6),(1e-8,2e-8,1e-10)]
    for aa,bb,tt in values:
        ref=reference(aa,bb,tt);stable=smooth(aa,bb,tt)
        tests.append(dict(a=aa,b=bb,t=tt,reference=mp.nstr(ref,55),stable_float=stable,stable_relative_error=rel(stable,ref),raw_relative_error=rel(raw(aa,bb,tt),ref)))
    quads=[]
    for aa,bb,tt in [(0,0,1),(0,2,.3),(.5,.5,.1),(.5,.50001,.1),(1,2,1),(8,8,1),(8,80,1)]:
        aaa,bbb,ttt=map(mp.mpf,map(str,(aa,bb,tt)))
        fn=lambda ss:ss*ss/((ss*ss+aaa*aaa)*(ss*ss+bbb*bbb))
        split=max(aaa,bbb,ttt)
        integ=2/mp.pi*mp.quad(fn,[ttt,split,mp.inf]) if split>ttt else 2/mp.pi*mp.quad(fn,[ttt,mp.inf])
        rr=reference(aa,bb,tt)
        assert abs(integ-rr)<mp.mpf('1e-85')*abs(rr)
        quads.append(dict(a=aa,b=bb,t=tt,relative_discrepancy=mp.nstr(abs(integ-rr)/abs(rr),8)))
    # Gauss nodes are independently generated, not copied from MATLAB.
    nodes,weights=np.polynomial.legendre.leggauss(8);nodes=(nodes+1)/2;weights=weights/2
    qtests=[]
    for aa,bb in [(8,8),(8,80),(8.0001,10000),(10,10),(100,100),(1e8,1e8)]:
        ref=2/mp.pi*mp.quad(lambda ss:ss*ss/((ss*ss+mp.mpf(str(aa))**2)*(ss*ss+mp.mpf(str(bb))**2)),[0,1])
        gauss=2/math.pi*sum(ww*xx*xx/((xx*xx+aa*aa)*(xx*xx+bb*bb)) for xx,ww in zip(nodes,weights))
        qtests.append(dict(a=aa,b=bb,t=1,reference=mp.nstr(ref,45),gauss8_relative_error=rel(gauss,ref),closed_relative_error=rel(q_closed(aa,bb,1),ref)))
    orbits=orbit_checks()
    out=dict(fixed_orbit_audit=orbits,status='PASS exact scalar checks and nondirected numerical cross-checks',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,executable=sys.executable,platform=platform.platform(),sympy=sp.__version__,mpmath=mp.__version__,numpy=np.__version__,precision_digits=100,workers=1,seed=None,script_sha256=sha(Path(__file__)),sources=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in SOURCES],symbolic_checks=exact,stable_V_checks=tests,independent_integral_checks=quads,short_Q_checks=qtests,boundary_threshold=dict(C_exact='3*sqrt(3)/(8*pi)',C_numeric=str(sp.N(D,35)),t_threshold_numeric=str(sp.N(sp.sqrt(D),35)),condition='sum ak=A>0, ak>=0, common t; t^2>C is sufficient for strict ray decrease on D'),scope='No MATLAB launch, no LP, no 2D integration, no rigorous floating enclosure, no general C2 upper bound.',elapsed_seconds=time.monotonic()-start)
    path=BASE/'computation/results/cauchy-smoothed-potential-independent/checks.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],symbolic_checks=len(exact),max_stable_V_relative_error=max(z['stable_relative_error'] for z in tests),max_raw_V_relative_error=max(z['raw_relative_error'] for z in tests),max_Q_gauss8_relative_error=max(z['gauss8_relative_error'] for z in qtests),boundary_threshold=out['boundary_threshold'],sources=out['sources'],fixed_orbit_audit={k:v for k,v in orbits.items() if k!='minimal_input' and k!='first_moment_and_compact_checks'},output=str(path),output_sha256=sha(path),elapsed_seconds=out['elapsed_seconds']),indent=2))
if __name__=='__main__':main()
