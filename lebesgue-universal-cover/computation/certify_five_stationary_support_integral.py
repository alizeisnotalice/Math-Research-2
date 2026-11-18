"""Check the two-angle stationary point with the support-function area integral."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,arb_series,ctx

def support_area(st,a,b):
    intervals,*_=st.model(a,b)
    start=0*a;twice=0*a
    for _,end,_,_,rho,c in intervals:
        cx,cy=c
        twice+=rho**2*(end-start)
        twice+=2*rho*(cx*(end.sin()-start.sin())-cy*(end.cos()-start.cos()))
        twice+=(cx**2-cy**2)*((2*end).sin()-(2*start).sin())/2
        twice-=cx*cy*((2*end).cos()-(2*start).cos())
        start=end
    return twice/2

def deriv(st,a,b):
    a0,b0=arb_series([a]),arb_series([b]);aa,bb=arb_series([a,1]),arb_series([b,1])
    f=support_area(st,aa,b0);g=support_area(st,a0,bb);d=support_area(st,aa,bb)
    mix=d[2]-f[2]-g[2]
    return (f[1],g[1]),((2*f[2],mix),(mix,2*g[2]))

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("contact",type=Path)
    a=p.parse_args();sys.path.insert(0,str(a.contact.resolve()))
    import certify_five_body_contact_stationary as st
    ctx.prec=224;ctx.cap=3
    X=[st.interval(*row)for row in st.BOX]
    mid=[st.q((F(l)+F(h))/2)for l,h in st.BOX];rad=[st.q((F(h)-F(l))/2)for l,h in st.BOX]
    topology=st.verify_topology(*X)
    g,H=deriv(st,*X);gm,Hm=deriv(st,*mid)
    det=lambda A:A[0][0]*A[1][1]-A[0][1]*A[1][0]
    assert H[0][0]>0 and det(H)>0
    d=det(Hm);C=[[Hm[1][1]/d,-Hm[0][1]/d],[-Hm[1][0]/d,Hm[0][0]/d]]
    C=[[arb(round(float(u.mid())*2**36))/2**36 for u in row]for row in C];assert abs(det(C))>0
    E=[[int(i==j)-sum(C[i][k]*H[k][j]for k in (0,1))for j in (0,1)]for i in (0,1)]
    contraction=max(sum(abs(v).abs_upper()for v in row)for row in E);assert contraction<1
    for i in (0,1):
        err=abs(sum(C[i][k]*gm[k]for k in (0,1)))+sum(abs(E[i][j])*rad[j]for j in (0,1))
        assert err<rad[i]
    value=support_area(st,*mid)
    err=sum(abs(gm[i])*rad[i]for i in (0,1))+sum(abs(H[i][j])*rad[i]*rad[j]/2 for i in (0,1)for j in (0,1))
    M=arb(value,err.abs_upper());ref=json.loads((a.contact/'five-body-contact-stationary.json').read_bytes())
    assert arb(ref['area_ball']).overlaps(M) and st.q(".83729500622102524")<M<st.q(".83729500622102525")
    print(json.dumps(dict(support_identity="A = integral(h^2 - derivative(h)^2)/2",pieces=len(st.WORD),
        stationary_angle_box=st.BOX,normal_partition=topology,contraction_upper=str(contraction),
        Hessian_determinant=str(det(H)),area_ball=str(M),Green_area_ball=ref['area_ball'],
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        chart_sha256=sha256((a.contact/"certify_five_body_contact_stationary.py").read_bytes()).hexdigest()),indent=2))
if __name__=="__main__":main()
