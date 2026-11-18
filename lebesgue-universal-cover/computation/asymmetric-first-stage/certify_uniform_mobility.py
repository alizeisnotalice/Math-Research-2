"""A sharp continuous necessary offset sum for unequal positive Gibbs tilts in either order."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from flint import arb,ctx
def ea(v):
    f=F(v);return arb(f.numerator)/f.denominator

def verify(ks,data,increasing):
    a=ea(data['third_curvature']); assert a<0
    cur=list(map(ea,data['sine_curvature']));assert len(cur)==len(ks)
    coef=[v/(1-k**2) for k,v in zip(ks,cur,strict=True)]
    ss=[1 if k%6==1 else -1 for k in ks]
    N=4096;mesh=2*arb.pi()/N
    def curve(t):
        return a*(3*t).sin()+sum(v*(k*t).sin() for k,v in zip(ks,cur,strict=True))
    B2=9*abs(a)+sum(k*k*abs(v) for k,v in zip(ks,cur,strict=True))
    curve_max=max(abs(curve(i*mesh)).upper() for i in range(N//2+1))+B2*mesh**2/8
    def base(t):
        st=t.sin()
        return -a*(3*t).sin()/8+sum(v*((k*t).sin()-sign*st) for k,v,sign in zip(ks,coef,ss,strict=True))
    def der(t):
        ct=t.cos()
        return -3*a*(3*t).cos()/8+sum(v*(k*(k*t).cos()-sign*ct) for k,v,sign in zip(ks,coef,ss,strict=True))
    g=lambda t:(5*t).sin()+t.sin()
    dg=lambda t:5*(5*t).cos()+t.cos()
    lo,hi=F('1.50'),F('1.54');deg=arb.pi()/180
    parts=100;roots=[];spans=[]
    pivot=1 if increasing else 0
    for i in range(parts):
        s=ea(lo+(hi-lo)*i/parts).union(ea(lo+(hi-lo)*(i+1)/parts))
        t=[(30+s+120*j)*deg for j in range(3)]
        cp=t[pivot];roots.append(-base(cp)/g(cp));spans.append(t)
    beta=roots[0]
    for val in roots[1:]:beta=beta.union(val)
    max_before=curve_max+24*abs(beta)
    R=F(int(((ea('.4999')/max_before).lower()*10**9).floor().unique_fmpz()),10**9)
    R=min(F(1),R);fac=ea(R)
    curvature=ea('.5')-fac*max_before
    assert curvature>0
    slopes=[];stables=[]
    for t in spans:
        cp=t[pivot]
        slopes.append((fac*(-1 if increasing else 1)*(der(cp)+beta*dg(cp))).lower())
        signs={0:1,2:-1} if increasing else {1:-1,2:1}
        for j,sg in signs.items():
            val=fac*sg*(base(t[j])+beta*g(t[j]))
            stables.append(val.lower())
    slope=min(slopes)
    stable=min(stables)
    coefficient=min(slope,(stable/(ea(hi-lo)*deg)).lower()).lower()
    assert coefficient>ea('.117')
    # The vertical reflection of the three even directions permutes j by (1,0,2).
    signA=([1,1,-1] if increasing else [1,-1,1])
    signB=([1,-1,-1] if increasing else [-1,-1,1])
    for mirror in (False,True):
        A,B=signA[:],signB[:]
        if mirror:A,B=([A[1],A[0],A[2]],[B[1],B[0],B[2]])
        for j in range(3):assert A[j]*B[(j+(-1 if mirror else 1))%3]<0
    return dict(order='s3<s5' if increasing else 's3>s5',critical_even_index=pivot,scale=str(R),
        beta_interval=str(beta),curvature_positive_lower=str(curvature.lower()),
        critical_slope_lower=str(slope),stable_absolute_lower=str(stable),coefficient_lower=str(coefficient))

def run(path):
    ctx.prec=192;q=json.loads(path.read_bytes())
    ks=q['frequencies'];assert ks==[k for k in range(5,192,2) if k%3]
    rows=[verify(ks,q['cases'][n],n=='increasing') for n in ('increasing','decreasing')]
    return dict(complete=True,interval=['1.50','1.54'],precision=ctx.prec,maximum_degree=191,
        coefficient_strictly_exceeds='.117',cases=rows,source_sha256=sha256(path.read_bytes()).hexdigest(),
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser();p.add_argument('data',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();r=run(a.data);a.out.write_text(json.dumps(r,indent=1)+'\n')
    for row in r['cases']:print(row['order'],row['coefficient_lower'])
