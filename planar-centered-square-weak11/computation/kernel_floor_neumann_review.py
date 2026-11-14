"""Independent fixed-kernel-family Neumann certificate verifier.
Rebuilds the entire 1182x1182 transpose feature matrix from certificate points,
using direct orbit expansion and explicit singular closed forms. Does not import
hl_core, kernel_family_floor, NumPy, SciPy, or call any matrix solve/inverse.
R is read as exact binary64 dyadic rationals; y0 is rational endpoint midpoint.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
from flint import arb, acb, arb_mat, fmpq, ctx
import hashlib,json,platform,struct,sys,time,flint
ctx.prec=192
start=time.time()
cp, kp, rp, op=map(Path,sys.argv[1:5])
cert=json.loads(cp.read_text()); candidate=json.loads(kp.read_text())
classes=[tuple(x) for x in cert['variables']['classes']]
assert classes==sorted(tuple(map(int,s.split(','))) for s in candidate['a'])
assert len(classes)==1181 and len(set(classes))==1181
for s in ['alpha','R','h']: assert candidate[s]==cert['variables'][s]
assert hashlib.sha256(kp.read_bytes()).hexdigest()==cert['candidate_sha256']
rows=cert['active_constraints']; n=len(classes)+1; assert len(rows)==n==1182
p,R,h=map(F,[candidate['alpha'],candidate['R'],candidate['h']])
# Self-contained necessary-row preflight; no unarchived full menu is needed.
assert len({(row['kind'],tuple(row.get('point',[]))) for row in rows})==n
for row in rows:
    assert row['kind'] in ('K','A')
    u,v=map(F,row['point'])
    assert u>=v>=0 and u+v>0
    if row['kind']=='K':
        assert row['rhs']==int(u+v<=1)
    else:
        assert v>0 and u+v!=R and row['rhs']==0

def A(q):
    q=F(q); return arb(fmpq(q.numerator,q.denominator))

def power(q, exponent):
    q=F(q); assert q>0
    return (A(exponent)*A(q).log()).exp()

@lru_cache(None)
def B(q):
    q=abs(q)
    if q>=2:return F(0)
    return (4-6*q*q+3*q*q*q)/6 if q<=1 else (2-q)**3/6
alpha=A(p)
kappa=-2/(alpha*A(1-p)*A(2-p)*A(3-p))
scale=power(h,-p)*kappa
@lru_cache(None)
def exponent_power(q):
    return power(q,3-p) if q else arb(0)
@lru_cache(None)
def LB(q):
    q=abs(q)
    return scale*sum(w*exponent_power(abs(q-k)) for k,w in zip(range(-2,3),(1,-4,6,-4,1)))

C=-8*(1-alpha).gamma()*(2*alpha).gamma()*(arb.pi()*alpha/2).cos()**2/(alpha+1).gamma()
Rpower=power(R,-p)
def hyp(a,b,c,z,**flags):return acb(z).hypgeom_2f1(acb(a),acb(b),acb(c),**flags).real
@lru_cache(None)
def tail_I(d):
    return Rpower*power(d,-p)/alpha-power(d,-2*p)*hyp(alpha,2*alpha,2*alpha+1,A(1-R/d),bc=True)/(2*alpha)

def outside_J(x,y):
    if y>=R:return arb(0)
    assert x+y>R and x>0 and y>0
    q=1+p; r=x+y; e=1-p
    def first_primitive(z):
        return power(z,e)/A(e)*hyp(A(e),A(q),A(e+1),A(z),ac=True)
    first=power(r,1-p-q)*(first_primitive(R/r)-first_primitive(y/r))
    e=p+q-1; c=x-y
    def second_primitive(z):
        return power(z,e)/A(e)*hyp(A(q),A(e),A(e+1),A(-c*z),bc=True)
    second=second_primitive(1/y)-second_primitive(1/R)
    L=R-y
    constant=-Rpower*(power(x-L,-p)-power(x+L,-p))/alpha
    return first+second+constant

def singular_generator(u,v):
    assert u>0 and v>0 and u+v!=R
    if u+v<R:
        return C*power(u+v,-2*p)+4*tail_I(R-u-v)+2*tail_I(R+u-v)+2*tail_I(R+v-u)
    return 2*outside_J(u,v)+2*outside_J(v,u)

def orbit(i,j):return sorted({(s*x,t*y) for x,y in [(i,j),(j,i)] for s in [-1,1] for t in [-1,1]})
orbits=[orbit(i,j) for i,j in classes]
maxindex=max(i for i,j in classes)
indices=range(-maxindex,maxindex+1)
matrix=arb_mat(n,n)
for col,row in enumerate(rows):
    u,v=map(F,row['point']); assert u>=v>=0 and u+v>0
    X,Y=u/h,v/h
    bx={i:B(X-i) for i in indices}; by={i:B(Y-i) for i in indices}
    if row['kind']=='K':
        assert row['rhs']==int(u+v<=1)
        matrix[0,col]=power(u+v,-p)-Rpower if u+v<R else arb(0)
        for j,orb in enumerate(orbits,1):
            val=sum(bx[a]*by[b] for a,b in orb)
            if val:matrix[j,col]=A(val)
    else:
        assert row['kind']=='A' and v>0 and u+v!=R and row['rhs']==0
        matrix[0,col]=singular_generator(u,v)
        for j,orb in enumerate(orbits,1):
            val=arb(0)
            for a,b in orb:
                if by[b]:val+=LB(X-a)*A(by[b])
                if bx[a]:val+=A(bx[a])*LB(Y-b)
            if not val.is_zero():matrix[j,col]=val
    assert matrix[0,col].is_finite()
    if (col+1)%300==0:print(json.dumps({'stage':'independent_matrix','done':col+1,'n':n}),flush=True)
mass=arb_mat(n,1)
mass[0,0]=alpha*power(R,2-p)/A(2-p)
for i,orb in enumerate(orbits,1):mass[i,0]=A(h*h*len(orb)/2)
y0=[(F(row['dual']['lo'])+F(row['dual']['hi']))/2 for row in rows]
assert all(q>0 for q in y0)
Y0=arb_mat(n,1,[A(q) for q in y0])
raw=rp.read_bytes();assert len(raw)==n*n*8
preconditioner=arb_mat(n,n,[A(F(t[0])) for t in struct.iter_unpack('<d',raw)])
build_time=time.time()-start
print(json.dumps({'stage':'independent_neumann_product','n':n,'build_seconds':build_time}),flush=True)
E=preconditioner*matrix
for i in range(n):E[i,i]-=1
q=max(sum(abs(E[i,j]) for j in range(n)).upper() for i in range(n))
assert q<1
residual=mass-matrix*Y0
correction=preconditioner*residual
eta=max(abs(correction[i,0]).upper() for i in range(n))
delta=eta/(1-q)
assert all(A(t)>delta for t in y0)
obj0=sum(F(row['rhs'])*t for row,t in zip(rows,y0))
lower=A(obj0)-sum(abs(row['rhs']) for row in rows)*delta
floor=F(1794671,500000)
assert lower>A(floor)
def ball(q):return {'display':str(q),'lo':str(q.lower().fmpq()),'hi':str(q.upper().fmpq())}
result={'scope':'fixed 1182-variable family floor; never C2 bound',
 'proof':'independent full matrix reconstruction and Neumann existence bound; no solve/inverse',
 'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':ctx.prec,
 'variable_count':n,'active_count':len(rows),'all_active_rows_necessary':True,
 'active_scope_preflight':'K/A labels, ordered nonnegative nonorigin coordinates, exact rhs, regular nonaxial A points; no full-menu dependency',
 'all_free_spline_equalities_certified':True,
 'q_upper':ball(q),'eta_upper':ball(eta),'uniform_solution_error_bound':ball(delta),
 'minimum_y0':str(min(y0)),'objective_lower':ball(lower),'certified_floor':str(floor),
 'build_seconds':build_time,'total_seconds':time.time()-start,
 'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [cp,kp,rp]},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
op.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
