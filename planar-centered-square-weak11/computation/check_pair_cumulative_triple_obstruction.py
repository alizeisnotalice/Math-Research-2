#!/usr/bin/env python3
"""Exact finite all-subset verification of a second-order cumulative strip defect."""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
from hashlib import sha256
import json, os, sys

PTS=[(F(-1),F(0),F(4)),(F(1),F(3,2),F(4)),(F(2),F(0),F(4))]
R=F(2);XCLIP=(F(0),F(1,2));BITS=44

def rectangles(indices,side):
    D=1<<BITS;out=[];examined=0
    for n in range(1,len(indices)+1):
      for chosen in combinations(indices,n):
        pts=[PTS[i] for i in chosen];mass=sum((p[2] for p in pts),F(0));q=mass*D*D/4
        k=isqrt(q.numerator//q.denominator)
        assert F(k*k)<=q<F((k+1)*(k+1));examined+=1
        t=min(R,F(k+(side=='outer'),D))
        a=min(p[0] for p in pts);b=max(p[0] for p in pts)
        c=min(p[1] for p in pts);d=max(p[1] for p in pts)
        x0=max(XCLIP[0],b-t);x1=min(XCLIP[1],a+t)
        y0=d-t;y1=c+t
        if x0<x1 and y0<y1:out.append((x0,x1,y0,y1))
    return out,examined

def area(rects):
    xs=sorted({v for r in rects for v in r[:2]});total=F(0)
    for left,right in zip(xs,xs[1:]):
      intervals=sorted((a,b) for x0,x1,a,b in rects if x0<right and left<x1)
      ylength=F(0)
      if intervals:
        lo,hi=intervals[0]
        for a,b in intervals[1:]:
          if a>hi:ylength+=hi-lo;lo,hi=a,b
          else:hi=max(hi,b)
        ylength+=hi-lo
      total+=(right-left)*ylength
    return total

def bounds(indices):
    ri,ni=rectangles(indices,'inner');ro,no=rectangles(indices,'outer')
    assert ni==no
    return (area(ri),area(ro)),ri,ro,ni

def radical_bounds(c0,c2=F(0),c3=F(0),c6=F(0),bits=48):
    lo=hi=F(c0);D=1<<bits
    for n,c in [(2,c2),(3,c3),(6,c6)]:
      k=isqrt(n*D*D);a=F(k,D);b=F(k+1,D)
      assert a*a<=n<b*b
      lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi

def main():
    assert os.path.ismount('/Volumes/data')
    data={};rects={};count=0
    for n in (1,2,3):
      for inds in combinations(range(3),n):
        interval,ri,ro,examined=bounds(inds);count+=examined
        data[inds]=interval;rects[inds]=(ri,ro)
    singles=[(0,),(1,),(2,)];pairs=[(0,1),(0,2),(1,2)]
    Slo=sum((data[t][0] for t in singles),F(0));Shi=sum((data[t][1] for t in singles),F(0))
    Plo=sum((data[t][0] for t in pairs),F(0));Phi=sum((data[t][1] for t in pairs),F(0))
    full=data[(0,1,2)]
    residual=(full[0]-Phi+Slo,full[1]-Plo+Shi)
    assert residual[0]>F(1,10)
    ui=[r for t in pairs for r in rects[t][0]]
    uo=[r for t in pairs for r in rects[t][1]]
    purepair_union=(area(ui),area(uo))
    beyond_pair_union=(full[0]-purepair_union[1],full[1]-purepair_union[0])
    assert beyond_pair_union[0]>F(1,10)
    # Independent radical formulas for full strip and the Mobius defect.
    analytic_full=radical_bounds(F(7,2),F(2),F(-3,2),F(-1))
    analytic_res=radical_bounds(F(-1,2),F(4),F(-3,2),F(-1))
    assert analytic_full[0]<=full[1] and full[0]<=analytic_full[1]
    assert analytic_res[0]<=residual[1] and residual[0]<=analytic_res[1]
    # Strict rational witness, not point sampling.
    assert F(23,16)**2>2 and F(3,2)**2>2
    assert F(31,32)>0 and F(47,32)>1
    r=F(13,8)
    assert max(F(3,2),F(49,32),F(25,16))<r
    assert F(12)/(4*r*r)==F(192,169)>1
    assert F(1,16)*F(1,16)==F(1,256)
    # Short analytic sign certificate for the exact radical residual.
    assert F(140,99)**2<2 and F(26,15)**2>3 and F(49,20)**2>6
    explicit=4*F(140,99)-F(1,2)-F(3,2)*F(26,15)-F(49,20)
    assert explicit==F(211,1980)>F(1,10)
    # Periodic/full-scale tail: [2,13/4) one period, then orbit count.
    assert F(12)/(4*F(2)**2)==F(3,4)<1
    assert F(3,4)/F(13,4)+F(3)/(F(13,4)**2)==F(87,169)<1
    source=Path(__file__).resolve()
    result={'status':'PASS','python':sys.version,'script_sha256':sha256(source.read_bytes()).hexdigest(),
      'arithmetic':'Fraction and directed rational root enclosures; no floating comparisons',
      'seed':None,'tolerance':None,'root_bits':BITS,'primitive_subsets_evaluated':count,
      'input':[[str(v) for v in p] for p in PTS],'h':'1/2','T':'8','R':'2',
      'areas':{''.join(str(i+1) for i in t):[str(a),str(b)] for t,(a,b) in data.items()},
      'single_sum':[str(Slo),str(Shi)],'pair_sum':[str(Plo),str(Phi)],
      'second_order_prediction':[str(Plo-Shi),str(Phi-Slo)],
      'purepair_union':[str(a) for a in purepair_union],
      'mobius_residual':[str(a) for a in residual],
      'beyond_purepair_union':[str(a) for a in beyond_pair_union],
      'explicit_residual_lower_bound':'211/1980','residual_greater_than':'1/10',
      'strict_witness_area':'1/256','strict_witness_radius':'13/8','strict_witness_average':'192/169',
      'scope':'specified three sources and endpoint-charged second-order cumulative template only; not all pair potentials'}
    out=source.parent/'results/pair-cumulative-triple-obstruction/exact-check.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','primitive_subsets_evaluated','second_order_prediction',
      'mobius_residual','beyond_purepair_union','explicit_residual_lower_bound')},indent=2))
    print(out)
if __name__=='__main__':main()
