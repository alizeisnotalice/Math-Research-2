"""Exact standard-library saturation replay for degrees 23 and 24."""
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json,sys,time

sys.dont_write_bytecode=True
from verify_degrees_21_22 import T,ZERO,setup,determinant,positive_low_rows,proof_exists_positive
from verify_degree_19 import Element,norm,cross,element
from verify_degrees_17_18 import P,H,Interval,real_bounds,evaluate,common_values

LO,HI=map(Q,real_bounds()['t_isolation'])
INTERVAL=Interval(LO,HI)


def generate(j):
    for abc in product('ABX',repeat=3):
        for eleven in 'ABXZ':
            early=''.join(abc)+eleven
            for twelve in ('AB' if j%2==0 else (None,)):
                for bits in product((0,1),repeat=j-19):
                    yield (early,twelve,bits)


def sequence(j,shape):
    early,twelve,bits=shape
    prefix=list(early)
    if j%2==0:prefix.append(twelve)
    prefix.extend({'X':'Z','Z':'F','A':'A','B':'B'}[i] for i in early[::-1])
    assert len(prefix)==j-15
    return prefix


def model(t,j,shape,values):
    A,B=setup(t)
    remainder=[]
    variables=iter(values)
    sequ=sequence(j,shape)
    cs=[1]*4+list(shape[2])
    assert len(sequ)==len(cs)
    for k,(side,c) in enumerate(zip(sequ,cs,strict=True),start=8):
        s=c-cross(A,B,k)
        if side=='A':a,b=s,0
        elif side=='B':a,b=0,s
        elif side in ('X','F'):
            a=next(variables)
            b=s-a
        else:
            assert side=='Z'
            a=b=0
            remainder.append(s)
        A.append(a);B.append(b)
    assert not list(variables)
    return remainder,A,B


def low_positive(early):
    ap,bp=positive_low_rows()
    for i,position in enumerate(early,start=8):
        if position in ('A','X'):ap.add(i)
        if position in ('B','X'):bp.add(i)
    return ap,bp


def pure_bound(j,shape):
    _,A,B=model(INTERVAL,j,shape,())
    for k in range(8,len(A)):
        for side,coeffs in (('A',A),('B',B)):
            a=coeffs[k]
            if isinstance(a,Interval) and a.hi<0:
                return side,k,str(a.hi)
    return None


def variables(j,shape):
    return tuple(k for k,side in enumerate(sequence(j,shape),start=8) if side in ('X','F'))


def system(j,shape):
    coords=variables(j,shape)
    dimension=len(coords)
    low=(0,)*dimension
    base=model(T,j,shape,low)[0]
    assert len(base)==dimension
    nonlinear=(j==24 and shape[0][0]=='X')
    # The final row 16 includes -A8*B8, contributing +u^2 (u=A8).
    lhs=[]
    for k in range(dimension):
        vals=low[:k]+(1,)+low[k+1:]
        lhs.append([a-b for a,b in zip(model(T,j,shape,vals)[0],base,strict=True)])
    if nonlinear:
        assert coords[0]==8
        lhs[0][-1]=lhs[0][-1]-1
    return nonlinear,coords,list(zip(*lhs,strict=True)),[-b for b in base]


def linear_test(coords,M,rhs):
    den=determinant(M)
    if den==ZERO:return None
    dn=norm(den)
    for i,k in enumerate(coords):
        num=determinant([r[:i]+(b,)+r[i+1:] for r,b in zip(M,rhs,strict=True)])
        nn=norm(num)
        if Q(nn,dn).denominator>1:return (k,nn,dn)
    return None


def mod_value(e):
    return evaluate(e.coefficients,-6)%P


def local_linear_rejection(M,rhs):
    for row,b in zip(M,rhs,strict=True):
        nonzero=[a for a in row if a!=ZERO]
        if not nonzero and b!=ZERO:
            return ('zero_linear_row',b.coefficients)
        if all(not mod_value(a) for a in row) and mod_value(b):
            return ('modular_linear_row',mod_value(b))
        if len(nonzero)==1:
            nn,nd=norm(b),norm(nonzero[0])
            if Q(nn,nd).denominator>1:
                return ('local_nonintegral_norm',nn,nd)
    return None


def quadratic(M,rhs):
    N=[r[1:] for r in M[:-1]]
    d=element(determinant(N))
    if d==ZERO:
        return None
    const=[]
    slope=[]
    for i in range(len(N)):
        const.append(determinant([r[:i]+(z,)+r[i+1:] for r,z in zip(N,rhs[:-1],strict=True)]))
        slope.append(determinant([r[:i]+(-z[0],)+r[i+1:] for r,z in zip(N,M[:-1],strict=True)]))
    b=d*M[-1][0]+sum(a*w for a,w in zip(M[-1][1:],slope,strict=True))
    c=-d*rhs[-1]+sum(a*w for a,w in zip(M[-1][1:],const,strict=True))
    return d,b,c,const,slope


def narrow_t():
    lo,hi=LO,HI
    for _ in range(65):
        mid=(lo+hi)/2
        if evaluate(H,mid)>0:
            lo=mid
        else:
            hi=mid
    assert evaluate(H,lo)>0>evaluate(H,hi)
    return Interval(lo,hi)


NARROW=narrow_t()


def quadratic_boxes(data,shape,depth=17):
    d,b,c,const,slope=data
    interval_of=lambda a:evaluate(a.coefficients,NARROW)
    d,b,c=map(interval_of,(d,b,c))
    variables=[(interval_of(a),interval_of(v)) for a,v in zip(const,slope,strict=True)]
    assert d.lo>0 or d.hi<0
    inverse=Interval(1/d.hi,1/d.lo)
    q_max=common_values(NARROW)[4].hi
    nodes=[(Interval(Q(0),q_max),0)]
    finished=[]
    rejected=Counter()
    count=0
    while nodes:
        u,k=nodes.pop()
        count+=1
        h=(d*u+b)*u+c
        if h.lo>0 or h.hi<0:
            rejected['quadratic_sign']+=1
            continue
        assigned=[u]+[(a+w*u)*inverse for a,w in variables]
        _,A,B=model(NARROW,24,shape,assigned)
        if any(isinstance(z,Interval) and z.hi<0 for z in A[8:]+B[8:]):
            rejected['coefficient_negative']+=1
            continue
        if k==depth:
            finished.append((u,A,B))
            continue
        mid=(u.lo+u.hi)/2
        nodes.append((Interval(mid,u.hi),k+1))
        nodes.append((Interval(u.lo,mid),k+1))
    return count,dict(rejected),finished


def late_prefix(j,A,B,digits,free):
    A=list(A);B=list(B)
    free=iter(free)
    for k,d in enumerate(digits,start=j-7):
        s=d-cross(A,B,k)
        if k==j-1:
            assert not list(free)
            return s
        if k in (j-2,j-5):
            a,b=0,s
        elif k in (j-3,j-4):
            a,b=s,0
        else:
            a=next(free)
            b=s-a
        A.append(a);B.append(b)
    raise ValueError('missing final convolution')


def positive_late(A,B,j):
    ap,bp=positive_low_rows()
    ap|={k for k,a in enumerate(A) if isinstance(a,Interval) and a.lo>0}
    bp|={k for k,b in enumerate(B) if isinstance(b,Interval) and b.lo>0}
    return tuple(k for k in range(j-7,j) if proof_exists_positive(ap,bp,k))


def late_pure(j,shape):
    _,A,B=model(-6,j,shape,[])
    _,at,bt=model(INTERVAL,j,shape,[])
    forced=positive_late(at,bt,j)
    zeros=(0,)*7
    c=late_prefix(j,A,B,zeros,(0,0))%P
    weights=(120,112,161,11,78,6,1)
    survivors=[]
    for bits in product((0,1),repeat=7):
        r=(c+sum(a*b for a,b in zip(weights,bits,strict=True)))%P
        for free in ((0,0),(1,0),(0,1)):
            assert late_prefix(j,A,B,bits,free)%P==r
        if r==0:
            survivors.append(bits)
            assert not all(bits[k-j+7] for k in forced)
    return {'constant':c,'forced':forced,'weights':weights,'residue_survivors':survivors}


def polynomial(values):
    a=[x%P for x in values]
    while a and not a[-1]:
        a.pop()
    return ModPoly(tuple(a))


@dataclass(frozen=True)
class ModPoly:
    coefficients: tuple[int,...]

    def __add__(self,other):
        other=as_mod_poly(other)
        a,b=self.coefficients,other.coefficients
        return polynomial((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b))))

    __radd__=__add__

    def __neg__(self):
        return polynomial(-a for a in self.coefficients)

    def __sub__(self,other):
        return self+-as_mod_poly(other)

    def __rsub__(self,other):
        return as_mod_poly(other)+-self

    def __mul__(self,other):
        other=as_mod_poly(other)
        a,b=self.coefficients,other.coefficients
        if not a or not b:
            return polynomial(())
        result=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                result[i+j]+=x*y
        return polynomial(result)

    __rmul__=__mul__


def as_mod_poly(a):
    return a if isinstance(a,ModPoly) else polynomial((a,))


def gcd(a,b):
    a=list(a.coefficients);b=list(b.coefficients)
    while b:
        norm=pow(b[-1],-1,P)
        while len(a)>=len(b):
            i=len(a)-len(b)
            c=a[-1]*norm
            for j,x in enumerate(b):
                a[i+j]=(a[i+j]-c*x)%P
            while a and not a[-1]:
                a.pop()
        a,b=b,a
    return polynomial(x*pow(a[-1],-1,P) for x in a)


def quadratic_late(data,shape,boxes):
    assert boxes
    t=-6
    d,b,c,const,slope=data
    d,b,c=map(mod_value,(d,b,c))
    assert d
    quin=polynomial((c,b,d))
    aa=polynomial((0,1))
    values=[aa]+[(mod_value(v)+mod_value(w)*aa)*pow(d,-1,P) for v,w in zip(const,slope,strict=True)]
    _,A,B=model(t,24,shape,values)
    A=list(map(as_mod_poly,A));B=list(map(as_mod_poly,B))
    cert=[]
    for u,ar,br in boxes:
        forced=positive_late(ar,br,24)
        assert forced==tuple(range(17,24))
        bits=(1,)*7
        r=late_prefix(24,A,B,bits,(0,0))
        assert r==late_prefix(24,A,B,bits,(1,0))==late_prefix(24,A,B,bits,(0,1))
        assert gcd(quin,r).coefficients==(1,)
        cert.append({'u_interval':[str(u.lo),str(u.hi)],'forced_rows':forced,
                     'quadratic_mod_163_ascending':quin.coefficients,
                     'late_residual_mod_163_ascending':r.coefficients,'gcd':1})
    return cert


def check(j):
    census=Counter()
    late=[]
    certified=[]
    quadratic_nodes=0
    for shape in generate(j):
        census['all']+=1
        a,b=low_positive(shape[0])
        if any(c==0 and proof_exists_positive(a,b,k) for k,c in enumerate(shape[2],start=12)):
            census['positive_forces_1']+=1
            continue
        if not variables(j,shape):
            bound=pure_bound(j,shape)
            if bound:
                census['negative_pure']+=1
                certified.append((shape,'negative_pure',bound))
            else:
                census['pure_modular']+=1
                late.append((shape,late_pure(j,shape)))
            continue
        nonlinear,coords,M,rhs=system(j,shape)
        if not nonlinear:
            cert=linear_test(coords,M,rhs)
            if cert:
                census['nonintegral_norm']+=1
                certified.append((shape,'nonintegral_norm',cert))
                continue
            cert=local_linear_rejection(M,rhs)
            if cert:
                census[cert[0]]+=1
                certified.append((shape,*cert))
                continue
            assert 'X' not in shape[0][:3]
            _,A,B=model(INTERVAL,j,shape,(0,)*len(coords))
            defect=cross(A,B,11)
            assert isinstance(defect,Interval) and defect.lo>1
            census['row11_overload']+=1
            certified.append((shape,'row11_overload',str(defect.lo)))
            continue
        data=quadratic(M,rhs)
        if data is None:
            cert=local_linear_rejection(M[:-1],rhs[:-1])
            assert cert and cert[0]=='local_nonintegral_norm'
            census['quadratic_singular_linear_norm']+=1
            certified.append((shape,'quadratic_singular_linear_norm',cert[1:]))
            continue
        n,leaves,boxes=quadratic_boxes(data,shape)
        quadratic_nodes+=n
        if boxes:
            evidence=quadratic_late(data,shape,boxes)
            census['quadratic_modular']+=1
        else:
            evidence=[]
            census['quadratic_intervals']+=1
        certified.append((shape,'quadratic',{'nodes':n,'rejected_leaves':leaves,'late_boxes':evidence}))
    return {'degree':j,'counts':dict(census),'pure_modular_128_codes_each':late,
            'quadratic_box_nodes':quadratic_nodes,'early_certificates':certified}


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--out',required=True,type=Path)
    options=p.parse_args()
    began=time.monotonic()
    assert (P,evaluate(H,-6)%P)==(163,0) and all(P%i for i in range(2,13))
    rec=[check(j) for j in (23,24)]
    for r in rec:
        print(r['degree'],r['counts'],r['quadratic_box_nodes'])
    names=(Path(__file__).name,'verify_degrees_21_22.py','verify_degree_19.py','verify_degrees_17_18.py')
    base=Path(__file__).parent
    result={'status':'verified','degrees':rec,'live':0,'python':sys.version,'seconds':time.monotonic()-began,
            'narrow_t':[str(NARROW.lo),str(NARROW.hi)],
            'source_sha256':{name:sha256((base/name).read_bytes()).hexdigest() for name in names}}
    with options.out.open('x') as f:f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':main()
