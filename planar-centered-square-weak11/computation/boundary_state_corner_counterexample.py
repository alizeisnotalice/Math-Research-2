#!/usr/bin/env python3
"""Exact finite common-radius square counterexample to x/y-marginal closure.
No floating computations enter proof checks. Deterministic; no dependencies.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import isqrt
import hashlib, json, platform, sys
from pathlib import Path

A = [(F(-3,4),F(-3,4),F(1)),(F(-1,4),F(-1,4),F(1))]
B = [(F(-3,4),F(-1,4),F(1)),(F(-1,4),F(-3,4),F(1))]
C = [(F(0),F(0),F(1,4))]
R = F(1)
BLOCK = (F(-1),F(0),F(-1),F(0))
U = (F(-3,16),F(-1,16),F(17,64),F(19,64))

def marginal(atoms,k):
    out={}
    for p in atoms: out[p[k]]=out.get(p[k],F(0))+p[2]
    return sorted(out.items())

def rectangles(atoms,bits,side,clip=None):
    """All closed source boxes, all enclosed mass, strict receiving boxes."""
    D=1<<bits
    xs=sorted({p[0] for p in atoms});ys=sorted({p[1] for p in atoms})
    out=set(); count=0
    for a,b in combinations_with_replacement(xs,2):
      for c,d in combinations_with_replacement(ys,2):
        mass=sum((m for x,y,m in atoms if a<=x<=b and c<=y<=d),F(0))
        if not mass: continue
        count+=1
        q=D*D*mass/4
        k=isqrt(q.numerator//q.denominator)
        assert F(k*k)<=q<F((k+1)*(k+1))
        t=min(R,F(k if side=='inner' else k+1,D))
        x0,x1,y0,y1=b-t,a+t,d-t,c+t
        if clip:
          x0=max(x0,clip[0]);x1=min(x1,clip[1]);y0=max(y0,clip[2]);y1=min(y1,clip[3])
        if x0<x1 and y0<y1:out.add((x0,x1,y0,y1))
    return sorted(out),count

def area(rects):
    xs=sorted({v for r in rects for v in r[:2]})
    ans=F(0)
    for x0,x1 in zip(xs,xs[1:]):
        intervals=sorted((y0,y1) for a,b,y0,y1 in rects if a<x1 and x0<b)
        total=F(0)
        if intervals:
            lo,hi=intervals[0]
            for a,b in intervals[1:]:
                if a>hi:total+=hi-lo;lo,hi=a,b
                else:hi=max(hi,b)
            total+=hi-lo
        ans+=(x1-x0)*total
    return ans

def interval(atoms,bits=32,clip=None):
    ri,n=rectangles(atoms,bits,'inner',clip);ro,n2=rectangles(atoms,bits,'outer',clip)
    return area(ri),area(ro),n

def mixed_interval(atoms,bits=32):
    # |H(atoms+C) \ (H(atoms) union H(C))|
    ji,_=rectangles(atoms+C,bits,'inner');jo,_=rectangles(atoms+C,bits,'outer')
    ai,_=rectangles(atoms,bits,'inner');ao,_=rectangles(atoms,bits,'outer')
    ci,_=rectangles(C,bits,'inner');co,_=rectangles(C,bits,'outer')
    return area(ji)-area(ao+co),area(jo)-area(ai+ci)

def event_values(atoms,x,y):
    return [(r,sum((m for a,b,m in atoms if max(abs(x-a),abs(y-b))<=r),F(0)))
      for r in sorted({max(abs(x-a),abs(y-b)) for a,b,m in atoms})]

def radical_interval(constant, c2=F(0), c5=F(0), bits=40):
    lo=hi=F(constant)
    for rad,c in [(2,c2),(5,c5)]:
        D=1<<bits;k=isqrt(rad*D*D)
        a,b=F(k,D),F(k+1,D)
        lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi

def check():
    assert sum(p[2] for p in A)==sum(p[2] for p in B)==2
    assert all(marginal(A,k)==marginal(B,k) for k in (0,1))
    assert sorted((x,-1-y,m) for x,y,m in A)==sorted(B)
    # Exact open rectangle with one fixed common-radius A+C witness.
    x0,x1,y0,y1=U;r=F(35,64)
    assert max(abs(x0+F(1,4)),abs(x1+F(1,4)),y1+F(1,4))<=r
    assert max(abs(x0),abs(x1),y1)<=r
    assert F(5,4)>4*r*r
    # All-radius exclusion in A alone, C alone, and B+C, from ordered distance events.
    assert y0+F(1,4)>F(1,2)
    assert y0+F(3,4)>1
    assert F(2)<4*(y0+F(3,4))**2
    assert F(1,4)<4*y0*y0
    first_B_min=x0+F(3,4)
    assert first_B_min>y1+F(1,4)
    assert F(5,4)<4*first_B_min**2
    assert F(9,4)<4*(y0+F(3,4))**2
    assert (x1-x0)*(y1-y0)==F(1,256)
    assert sum(p[2] for p in A+C)==F(9,4)<4*R*R
    res={}
    for name,atoms in [('A',A),('B',B),('A+C',A+C),('B+C',B+C)]:
      lo,hi,n=interval(atoms)
      res[name]={'area_lower':str(lo),'area_upper':str(hi),'area_decimal':float((lo+hi)/2),'nonzero_boxes':n}
    for name,atoms in [('A',A),('B',B)]:
      lo,hi,n=interval(atoms,clip=BLOCK)
      res[name]['internal_area_lower']=str(lo);res[name]['internal_area_upper']=str(hi)
      ml,mh=mixed_interval(atoms)
      res[name]['pure_mixed_lower']=str(ml);res[name]['pure_mixed_upper']=str(mh)
      res[name]['pure_mixed_decimal']=float((ml+mh)/2)
    alo,ahi,_=interval(A+C);blo,bhi,_=interval(B+C)
    assert ahi<blo or bhi<alo
    res['joined_area_difference_interval']=[str(alo-bhi),str(ahi-blo)]
    ma=mixed_interval(A);mb=mixed_interval(B)
    res['pure_mixed_difference_interval']=[str(ma[0]-mb[1]),str(ma[1]-mb[0])]
    res['pure_mixed_different']=ma[1]<mb[0] or mb[1]<ma[0]
    analytic={
      'A':radical_interval(F(13,4),c2=F(-1)),
      'A internal':radical_interval(F(19,8),c2=F(-1)),
      'A+C':radical_interval(F(67,16),c2=F(-1),c5=F(-3,8)),
      'B+C':radical_interval(F(19,4),c2=F(-1,2),c5=F(-7,8)),
      'A pure mixed':radical_interval(F(15,16),c5=F(-3,8)),
      'B pure mixed':radical_interval(F(2),c5=F(-7,8)),
    }
    checks=[(analytic['A'],interval(A)[:2]),(analytic['A internal'],interval(A,clip=BLOCK)[:2]),
      (analytic['A+C'],(alo,ahi)),(analytic['B+C'],(blo,bhi)),
      (analytic['A pure mixed'],ma),(analytic['B pure mixed'],mb)]
    assert all(a[0]<=b[1] and b[0]<=a[1] for a,b in checks)
    assert F(89,40)**2<5 and F(7,5)**2<2 and F(9,4)**2>5
    res['analytic_formula_intervals']={k:[str(a),str(b)] for k,(a,b) in analytic.items()}
    res['analytic_geometric_intervals_overlap']='PASS'
    p=Path(__file__).resolve()
    return {'python':sys.version,'platform':platform.platform(),'script_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
      'lambda':'1','R':'1','dyadic_bits':32,'seed':None,'tolerance':None,
      'A':[[str(v) for v in p] for p in A],'B':[[str(v) for v in p] for p in B],'C':[[str(v) for v in p] for p in C],
      'witness_rectangle':[str(v) for v in U],'witness_radius':'35/64','witness_average':'256/245',
      'results':res,'exact_checks':'PASS'}
if __name__=='__main__':
    result=check();s=json.dumps(result,indent=2)
    if len(sys.argv)>1:Path(sys.argv[1]).write_text(s+'\n')
    print(s)
