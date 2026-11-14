"""Independent expanded-tensor Arb evaluation of the one-point ray cap.
No hl_core imports; finite differences are applied to each B-spline separately.
The singular closed form is reimplemented from the displayed analytic formula.
This is not a C1/C2/C4 continuum replay.
"""
import hashlib, json, platform, sys
from fractions import Fraction as F
from pathlib import Path
import flint
from flint import arb, acb, fmpq, ctx
ctx.prec=192
root=Path(sys.argv[1]).resolve()
out=Path(sys.argv[2]).resolve()

def A(q):
    q=F(q)
    return arb(fmpq(q.numerator,q.denominator))

def pow_(x,q):
    if x==0: return arb(0)
    return (A(q)*A(x).log()).exp()

def B(q):
    q=abs(F(q))
    if q>=2: return F(0)
    if q<=1: return (4-6*q*q+3*q*q*q)/6
    return (2-q)**3/6

weights=(1,-4,6,-4,1)
def LB(q,p):
    kappa=-2/(A(p)*A(1-p)*A(2-p)*A(3-p))
    return kappa*sum(w*pow_(abs(q-k),3-p) for k,w in zip(range(-2,3),weights))

def spline_A(d,u,v):
    p,h,N=F(d['alpha']),F(d['h']),int(d['N'])
    coeff={tuple(map(int,k.split(','))):F(val) for k,val in d['a'].items()}
    def coef(i,j): return coeff.get((max(abs(i),abs(j)),min(abs(i),abs(j))),F(0))
    # Direct tensor expansion over all positions, unlike the original D_m evaluator.
    b_u={i:B(u/h-i) for i in range(-N,N+1)}
    b_v={j:B(v/h-j) for j in range(-N,N+1)}
    x=sum(A(sum(coef(i,j)*b_v[j] for j in b_v if b_v[j]))*LB(u/h-i,p)
          for i in range(-N,N+1))
    y=sum(A(sum(coef(i,j)*b_u[i] for i in b_u if b_u[i]))*LB(v/h-j,p)
          for j in range(-N,N+1))
    return pow_(h,-p)*(x+y)

def singular_A(d,u,v):
    p,R=F(d['alpha']),F(d['R'])
    assert u>0 and v>0 and u+v<R
    alpha=A(p)
    C=-8*(1-alpha).gamma()*(2*alpha).gamma()*(arb.pi()*alpha/2).cos()**2/(alpha+1).gamma()
    def I(dist):
        z=A(-(R-dist)/dist)
        hg=acb(z).hypgeom_2f1(acb(alpha),acb(2*alpha),acb(2*alpha+1),bc=True).real
        return pow_(R,-p)*pow_(dist,-p)/alpha-pow_(dist,-2*p)*hg/(2*alpha)
    return C*pow_(u+v,-2*p)+4*I(R-u-v)+2*I(R+u-v)+2*I(R+v-u)

def ball(q):
    return {'display':str(q),'lo':str(q.lower().fmpq()),'hi':str(q.upper().fmpq())}

base=json.loads((root/'scaled.json').read_text())
cur=json.loads((root/'candidate.json').read_text())
assert all(base[k]==cur[k] for k in ['alpha','R','h','N','c'])
diff={k:F(cur['a'].get(k,'0'))-F(base['a'].get(k,'0')) for k in base['a']|cur['a']}
assert {k:v for k,v in diff.items() if v}=={'0,0':-F(1,200)}
u,v=F(35,64),F(5,256)
p,h,R,c=map(F,[cur['alpha'],cur['h'],cur['R'],cur['c']])
sing=singular_A(cur,u,v)
AK=A(c)*sing+spline_A(cur,u,v)
AP=pow_(h,-p)*(LB(u/h,p)*A(B(v/h))+A(B(u/h))*LB(v/h,p))
AB=A(c)*sing+spline_A(base,u,v)
assert (AK-(AB-A(F(1,200))*AP)).contains(0)
assert AK>0 and AP>A(F(9,100)) and AK<A(F(23,25000))
cap=AK/(2048*AP)
assert cap<A(F(1,200000))
coeff={tuple(map(int,k.split(','))):F(val) for k,val in cur['a'].items()}
N=int(cur['N'])
def coef(i,j): return coeff.get((max(abs(i),abs(j)),min(abs(i),abs(j))),F(0))
mass=A(c)*A(p)*pow_(R,2-p)/A(2-p)+A(h*h/2)*A(sum(coef(i,j) for i in range(-N,N+1) for j in range(-N,N+1)))
assert mass<A(F(3590865312933,10**12))
# Entirely rational coarse endpoint deduction, independent of printed decimals.
coarse=F(23,25000)/(2048*F(9,100))
assert coarse<F(1,200000)
# General-order leading coefficients; coefficients omit c and the shown u powers.
a=F(6,5)
minus=-4*p/(a*(a-1)); plus=4/(p-1)
assert a-p==F(1,40) and minus==-F(235,12) and plus==F(160,7)
assert -4*p/(p*(p-1))+plus==0
result={'scope':'independent direct tensor expansion of a single necessary C4 row; analytic asymptotic coefficients only',
 'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':ctx.prec,
 'inputs':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['scaled.json','candidate.json']},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'q':[str(u),str(v)],'AK_current':ball(AK),'AP':ball(AP),'AK_base':ball(AB),
 'additional_halfmass_cap':ball(cap),'current_halfmass':ball(mass),'coarse_exact_cap':str(coarse),
 'coarse_strict_less':'1/200000',
 'asymptotic_exact':{'p':str(p),'a':str(a),'a_minus_p':str(a-p),'negative_coefficient':str(minus),'positive_coefficient':str(plus),'at_a_equal_p_leading_sum':'0'},
 'full_C1_C2_C4_replayed':False}
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
