#!/usr/bin/env python3
"""R60 bounded independent rectangle/Jensen verification; no directed certificate."""
from pathlib import Path
import datetime,hashlib,json,math,os,sys,time
import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_jacobi
import sympy as sp
BASE=Path(__file__).resolve().parent.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def xlog(x):return x*mp.log(x) if x else mp.mpf(0)
def H(x,y):
    a,b=abs(x),abs(y)
    return mp.sign(x)*mp.sign(y)*(xlog(a+b)-xlog(a)-xlog(b))
def Urect(l,h,p,q):return H(h,q)-H(l,q)-H(h,p)+H(l,p)
def angle(s,l,h):
    if not s:
        return mp.pi if l<0<h else mp.pi/2 if l==0 or h==0 else mp.mpf(0)
    return mp.atan2(s*(h-l),s*s+l*h)
def D(u,l,h,t):return (h-l)/t if not u else mp.atan2(t*u*(h-l),t*t+l*h*u*u)/u
def Hfloat(v,l,h,p,q,t):
    u=np.sqrt(v)
    def dd(a,b):
        res=np.empty_like(u);nz=u!=0
        res[nz]=np.arctan2(t*u[nz]*(b-a),t*t+a*b*u[nz]*u[nz])/u[nz]
        res[~nz]=(b-a)/t
        return res
    return dd(l,h)*dd(p,q)
def gauss(n):
    x,w=roots_jacobi(n,0,-.5)
    return (x+1)/2,w/(2*math.sqrt(2))
def radau(n):
    if n==1:return np.array([0.]),np.array([1.])
    x,w=roots_jacobi(n-1,0,.5);v=(x+1)/2;ww=w/(4*math.sqrt(2)*v)
    return np.r_[0.,v],np.r_[1-sum(ww),ww]
def uniform_radau(n):
    if n==1:return np.array([0.]),np.array([1.])
    x,w=roots_jacobi(n-1,0,1);v=(x+1)/2;ww=w/(4*v)
    return np.r_[0.,v],np.r_[1-sum(ww),ww]
def split_sandwich(n,l,h,p,q,t):
    scale=max(abs(l),abs(h),abs(p),abs(q));eps=min(1.,(t/scale)**2) if scale else 1.
    vg,wg=gauss(n);vr,wr=radau(n)
    lower=math.sqrt(eps)*np.dot(wg,Hfloat(eps*vg,l,h,p,q,t))
    upper=math.sqrt(eps)*np.dot(wr,Hfloat(eps*vr,l,h,p,q,t))
    ng,wg=np.polynomial.legendre.leggauss(n);ng=(ng+1)/2;wg=wg/2;nr,wr=uniform_radau(n)
    a=eps;pieces=1
    while a<1:
        b=min(4*a,1.);v=a+(b-a)*ng;w=a+(b-a)*nr
        lower+=(b-a)*np.dot(wg,Hfloat(v,l,h,p,q,t)/(2*np.sqrt(v)))
        upper+=(b-a)*np.dot(wr,Hfloat(w,l,h,p,q,t)/(2*np.sqrt(w)))
        a=b;pieces+=1
    return 2*t/math.pi*lower,2*t/math.pi*upper,pieces
def breaks(t,ends):
    xs={mp.mpf(0),t}
    for x in ends:
        if 0<abs(x)<t:xs.add(abs(x))
    return sorted(xs)
def main():
    assert os.path.ismount('/Volumes/data')
    start=time.monotonic();mp.mp.dps=90
    a,b=sp.symbols('a b',positive=True)
    primitive=(a+b)*sp.log(a+b)-a*sp.log(a)-b*sp.log(b)
    assert sp.simplify(sp.diff(primitive,a,b)-1/(a+b))==0
    cases=[('origin_corner','0','1','0','1','.3'),('crossing','-.25','.5','-.5','2','.1'),('far','2','3','4','5','.2'),('thin_source_edge','.999','1.001','-.001','.001','.03'),('axis_strip','-8','-7.5','0','.01','1')]
    rows=[]
    for name,*vals in cases:
        l,h,p,q,t=map(mp.mpf,vals)
        points={mp.mpf(0),mp.mpf(1)}
        for e in (l,h,p,q):
            if abs(e)>t:points.add(t/abs(e))
        for x,y in [(l,h),(p,q)]:
            if x*y<0 and t/mp.sqrt(-x*y)<1:points.add(t/mp.sqrt(-x*y))
        V=2*t/mp.pi*mp.quad(lambda u:D(u,l,h,t)*D(u,p,q,t),sorted(points))
        Q=2/mp.pi*mp.quad(lambda s:angle(s,l,h)*angle(s,p,q),breaks(t,(l,h,p,q)))
        U=Urect(l,h,p,q)
        discrepancy=abs(U-Q-V)
        assert discrepancy<mp.mpf('1e-75')*max(1,abs(U)),(name,discrepancy)
        n=256;g=[D(mp.mpf(k)/n,l,h,t)*D(mp.mpf(k)/n,p,q,t) for k in range(n+1)]
        assert all(g[k]>=g[k+1] for k in range(n))
        lower=2*t/mp.pi*sum(g[1:])/n;upper=2*t/mp.pi*sum(g[:-1])/n
        assert lower<=V<=upper
        gap_formula=2*t/(mp.pi*n)*(g[0]-g[-1]);assert abs(upper-lower-gap_formula)<mp.mpf('1e-75')
        quadrature=[]
        for n in [4,8,16,32]:
            vg,wg=gauss(n);vr,wr=radau(n)
            args=tuple(map(float,(l,h,p,q,t)))
            G=2*float(t)/math.pi*np.dot(wg,Hfloat(vg,*args));R=2*float(t)/math.pi*np.dot(wr,Hfloat(vr,*args))
            momentG=max(abs(np.dot(wg,vg**k)-1/(2*k+1)) for k in range(2*n))
            momentR=max(abs(np.dot(wr,vr**k)-1/(2*k+1)) for k in range(2*n-1))
            assert min(wg)>0 and min(wr)>0 and max(momentG,momentR)<1e-12
            # Floating quadratures are evidence only; theorem determines sides.
            quadrature.append(dict(n=n,Gauss_V_lower_numeric=float(G),left_Radau_V_upper_numeric=float(R),width_over_reference=float((mp.mpf(R)-mp.mpf(G))/V),Gauss_relative_error=float((V-mp.mpf(G))/V),Radau_relative_error=float((mp.mpf(R)-V)/V),moment_residual_G=momentG,moment_residual_R=momentR))
        floatQ=[]
        for n in [16,32]:
            x,w=np.polynomial.legendre.leggauss(n);ss=float(t)*(x+1)/2
            def aa(ll,hh):return np.arctan2(ss*(float(hh)-float(ll)),ss*ss+float(ll*hh))
            approx=float(t)/math.pi*np.dot(w,aa(l,h)*aa(p,q))
            floatQ.append(dict(n=n,Q_numeric=float(approx),relative_error=float((mp.mpf(approx)-Q)/Q)))
        L8,U8,pieces=split_sandwich(8,*map(float,(l,h,p,q,t)))
        split=dict(n=8,piece_count=pieces,lower_numeric=float(L8),upper_numeric=float(U8),width_over_reference=float((mp.mpf(U8)-mp.mpf(L8))/V))
        rows.append(dict(name=name,split_sandwich=split,direct_Q_GL=floatQ,intervals=[vals[:2],vals[2:4]],t=vals[4],U_rectangle=mp.nstr(U,60),V_rectangle=mp.nstr(V,60),Q_rectangle=mp.nstr(Q,60),U_minus_V_minus_Q=mp.nstr(discrepancy,8),monotone_Riemann_n=256,monotone_Riemann_lower=mp.nstr(lower,50),monotone_Riemann_upper=mp.nstr(upper,50),quadrature=quadrature))
    out=dict(status='PASS bounded exact-symbolic and nondirected rectangle checks',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,executable=sys.executable,sympy=sp.__version__,mpmath=mp.__version__,numpy=np.__version__,scipy=scipy.__version__,precision_digits=90,workers=1,random_seed=None,script_sha256=sha(Path(__file__)),symbolic_check='mixed derivative of U rectangle primitive exactly 1/(a+b)',cases=rows,scope='Five deterministic rectangles, quadrature/moment scalar checks only; no MATLAB, LP, candidate, directed enclosure or new C2 bound.',elapsed_seconds=time.monotonic()-start)
    p=BASE/'computation/results/r060-cauchy-rectangle-jensen/checks.json';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],case_count=len(rows),script_sha256=out['script_sha256'],output_sha256=sha(p),output=str(p),elapsed_seconds=out['elapsed_seconds'],n16_widths={r['name']:r['quadrature'][2]['width_over_reference']for r in rows},split8_widths={r['name']:r['split_sandwich']['width_over_reference']for r in rows},Q_GL16_errors={r['name']:r['direct_Q_GL'][0]['relative_error']for r in rows}),indent=2))
if __name__=='__main__':main()
