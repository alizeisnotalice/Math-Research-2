"""Uniform two-cap offset lower bound by Minkowski interpolation of rational trigonometric data."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from flint import arb,ctx
def ea(v):
    f=F(v);return arb(f.numerator)/f.denominator

def verify(path):
    dat=json.loads(path.read_bytes());focus=dat['focus'];assert focus in (3,5)
    ks=dat['frequencies'];assert ks==[k for k in range(5,192,2) if k%3]
    coef=[ea(z) for z in dat['coefficients']]
    m=len(ks);a=coef[0];cc,ss=coef[1:1+m],coef[1+m:]
    assert (len(cc),len(ss))==(m,m) and a<0
    D2=9*abs(a)+sum(k*k*(abs(c)+abs(v)) for k,c,v in zip(ks,cc,ss,strict=True))
    N=4096
    def curve(t):
        return a*(3*t).sin()+sum(c*(k*t).cos()+v*(k*t).sin() for k,c,v in zip(ks,cc,ss,strict=True))
    # The curvature perturbation is anti-periodic at pi. Linear interpolation
    # contributes at most D2 * mesh^2/8 between consecutive exact grid values.
    mesh=2*arb.pi()/N
    peaks=[abs(curve(i*mesh)).upper() for i in range(N//2+1)]
    cmax=max(peaks)+D2*mesh**2/8
    raw=min(arb(1),ea('.4999')/cmax)
    mult=F(int((raw.lower()*10**10).floor().unique_fmpz()),10**10);fac=ea(mult)
    curv=ea('.5')-fac*cmax
    assert curv>0
    def support(t):
        v=-a*(3*t).sin()/8
        ct,st=t.cos(),t.sin()
        for k,c,y in zip(ks,cc,ss,strict=True):
            sign=1 if k%6==1 else -1
            v+=(c*((k*t).cos()-ct)+y*((k*t).sin()-sign*st))/(1-k*k)
        return fac*v
    cfg=[(3,'31.50',[1,1,-1]),(5,'31.54',[1,-1,-1]),(3,'28.50',[1,1,-1]),(5,'28.46',[-1,1,-1])]
    bounds={3:[],5:[]}
    for cap,r,signs in cfg:
        for j,sg in enumerate(signs):
            ang=(ea(r)+120*j)*arb.pi()/180
            bounds[cap].append((sg*support(ang)).lower())
    actual={k:min(vals) for k,vals in bounds.items()}
    assert all(x>0 for x in actual.values())
    return dict(focus=focus,scale=str(mult),curvature_positive_lower=str(curv.lower()),curvature_global_upper=str(cmax.upper()),
      signed_caps={k:str(val) for k,val in actual.items()},source_sha256=sha256(path.read_bytes()).hexdigest()),actual[focus]

def main(p3,p5):
    ctx.prec=192
    X,a=verify(p3);Y,b=verify(p5)
    M=min(a,b)
    budget=F(int((M.lower()*10**12).floor().unique_fmpz()),10**12)
    rt=arb(3).sqrt();deg=arb.pi()/180
    def area(e3,e5):
        def cap(s,e):
            c=(ea(s)*deg).cos()
            return rt*(c/rt-ea('.5')-e)**2/(4*c*c-3)
        return rt/2-cap('1.50',e3)-cap('1.54',e5)
    base=area(arb(0),arb(0))
    low=min(area(ea(budget),arb(0)).lower(),area(arb(0),ea(budget)).lower())
    assert low>base
    return dict(complete=True,precision=ctx.prec,cosine_sine_degrees=191,curvature_mesh=4096,
      data=[X,Y],necessary_sum_strict_lower=str(budget),min_first_stage_area_lower=str(low),
      min_added_area_lower=str((low-base.upper()).lower()),verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser()
    for k in ('p3','p5','out'):p.add_argument(k,type=Path)
    a=p.parse_args();r=main(a.p3,a.p5)
    a.out.write_text(json.dumps(r,indent=1)+'\n')
    print(r['necessary_sum_strict_lower'],r['min_first_stage_area_lower'],r['min_added_area_lower'])
