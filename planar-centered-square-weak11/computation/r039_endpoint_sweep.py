#!/usr/bin/env python3
"""R39 conditional endpoint-order and dyadic y-union slope enclosure.

Inputs must already be the COMPLETE active-box set on a source-stable slab;
this module does not certify that completeness, rejection margins, or tails.
Expressions are (rational offset at s0, rational velocity, sign, radicand).
Any unresolved comparison fails closed rather than guessing an order.
"""
from fractions import Fraction as F
from math import isqrt
from functools import cmp_to_key
from r039_interval_primitives import sqrt_interval

def expr(a,v=0,sgn=0,q=0):
    a,v,q=F(a),F(v),F(q)
    if q<0 or sgn not in (-1,0,1): raise ValueError('invalid radical')
    if not q or not sgn: return (a,v,0,F(0))
    n,d=isqrt(q.numerator),isqrt(q.denominator)
    if n*n==q.numerator and d*d==q.denominator:
        return (a+sgn*F(n,d),v,0,F(0))
    return (a,v,sgn,q)

def interval(e,bits):
    a,_,sign,q=e
    if not sign: return a,a
    lo,hi=sqrt_interval(q,bits)
    return (a+lo,a+hi) if sign==1 else (a-hi,a-lo)

def compare(a,b,bits=80):
    if a==b: return 0
    al,ah=interval(a,bits); bl,bh=interval(b,bits)
    if ah<bl: return -1
    if bh<al: return 1
    # Equality at s0 with differing velocities is a wall, not a persistent tie.
    if a[0]==b[0] and a[2:]==b[2:]:
        raise ValueError('nonpersistent endpoint collision at center')
    raise ValueError('unresolved endpoint comparison; raise precision or change center')

def union_length(intervals):
    z=sorted((a,b) for a,b in intervals if a<b)
    if not z: return F(0)
    total=F(0); a,b=z[0]
    for c,d in z[1:]:
        if c>b: total+=b-a; a,b=c,d
        else: b=max(b,d)
    return total+b-a

class UnionTree:
    """Exact dynamic union of fixed rational intervals."""
    def __init__(self, intervals):
        self.ys=sorted({z for a,b in intervals for z in (a,b)})
        self.index={z:i for i,z in enumerate(self.ys)}
        self.n=max(0,len(self.ys)-1)
        self.cover=[0]*(4*self.n+4); self.length=[F(0)]*(4*self.n+4)
    def update(self,a,b,delta):
        if a>=b: return
        aa,bb=self.index[a],self.index[b]
        def change(node,lo,hi):
            if aa<=lo and hi<=bb: self.cover[node]+=delta
            else:
                mid=(lo+hi)//2
                if aa<mid: change(2*node,lo,mid)
                if bb>mid: change(2*node+1,mid,hi)
            if self.cover[node]<0: raise ValueError('inconsistent interval events')
            if self.cover[node]: self.length[node]=self.ys[hi]-self.ys[lo]
            elif hi-lo==1: self.length[node]=F(0)
            else: self.length[node]=self.length[2*node]+self.length[2*node+1]
        change(1,0,self.n)
    def value(self): return self.length[1] if self.n else F(0)

def slope_certificate(boxes,xleft,xright,yleft,yright,bits=80):
    """Return conditional slope enclosure and sufficient endpoint stability radius.

    Each box is (xlo expression, xhi expression, ylo expression, yhi expression).
    y velocities must be zero. Domain bounds are rational. Empty clipped x/y boxes
    must already have been safely excluded for the candidate slab.
    """
    xl,xr,yl,yr=map(F,(xleft,xright,yleft,yright))
    if xl>=xr or yl>=yr: raise ValueError('invalid domain')
    left,right=expr(xl),expr(xr)
    allx={left,right}; raw_refs={left:[('domain','left')],right:[('domain','right')]}
    for i,(l,r,y0,y1) in enumerate(boxes):
        if y0[1] or y1[1]: raise ValueError('y endpoints must be fixed')
        allx.update((l,r))
        raw_refs.setdefault(l,[]).append((i,'left'))
        raw_refs.setdefault(r,[]).append((i,'right'))
    order=sorted(allx,key=cmp_to_key(lambda a,b: compare(a,b,bits)))
    rank={e:i for i,e in enumerate(order)}
    delta=None; minpair=None
    for a,b in zip(order,order[1:]):
        al,ah=interval(a,bits); bl,bh=interval(b,bits)
        gap=bl-ah
        if gap<=0: raise ValueError('unresolved ordered endpoint gap')
        dv=abs(b[1]-a[1])
        if dv:
            t=gap/(2*dv)
            if delta is None or t<delta: delta,minpair=t,(a,b)
    events={i:[] for i in range(len(order))}
    ybounds=[]
    for i,(l,r,y0,y1) in enumerate(boxes):
        lo=max(rank[l],rank[left]); hi=min(rank[r],rank[right])
        if lo>=hi: raise ValueError('empty x clip: needs rejection certificate')
        a,b=interval(y0,bits); c,d=interval(y1,bits)
        inner=(max(yl,b),min(yr,c)); outer=(max(yl,a),min(yr,d))
        if inner[0]>=inner[1]: raise ValueError('y clip not certified positive at chosen precision')
        ybounds.append((inner,outer)); events[lo].append((1,i)); events[hi].append((-1,i))
    inner=UnionTree([z[0] for z in ybounds]); outer=UnionTree([z[1] for z in ybounds])
    beta_lo=beta_hi=F(0); diagnostics=[]
    for j in range(len(order)):
        before=(inner.value(),outer.value())
        for sign,i in events[j]:
            inner.update(*ybounds[i][0],sign); outer.update(*ybounds[i][1],sign)
        after=(inner.value(),outer.value())
        velocity=order[j][1]
        if velocity:
            diff=(before[0]-after[1],before[1]-after[0])
            contribution=(velocity*diff[0],velocity*diff[1]) if velocity>0 else (velocity*diff[1],velocity*diff[0])
            if contribution!=(0,0):
                diagnostics.append({'endpoint':order[j],'raw_box_endpoint_refs':raw_refs[order[j]],'enter_leave_box_indices':events[j],
                                    'y_before':before,'y_after':after,'contribution':contribution})
        if j==len(order)-1: break
        dv=order[j+1][1]-order[j][1]
        if dv:
            l,h=after
            if dv>0: beta_lo+=dv*l; beta_hi+=dv*h
            else: beta_lo+=dv*h; beta_hi+=dv*l
    assert beta_lo<=beta_hi
    return {'beta_lower':beta_lo,'beta_upper':beta_hi,'endpoint_radius':delta,
            'endpoint_count':len(order),'box_count':len(boxes),'limiting_pair':minpair,'bits':bits,
            'moving_endpoint_diagnostics':diagnostics}

def selftest():
    e=expr
    # One rectangle grows rightwards from length 1 to length 2, height 2.
    r=slope_certificate([(e(1),e(2,1),e(1),e(3))],0,5,0,5)
    assert (r['beta_lower'],r['beta_upper'])==(2,2)
    # Overlap must be counted as a y-union, not a sum.
    r=slope_certificate([(e(1),e(2,1),e(0),e(2)),(e(1),e(2,1),e(1),e(3))],0,5,0,5)
    assert (r['beta_lower'],r['beta_upper'])==(3,3)
    try: slope_certificate([(e(1),e(1,1),e(0),e(1))],0,5,0,5)
    except ValueError: pass
    else: raise AssertionError('center wall was not rejected')
    # Two atoms with weights 9,16, level 1 and R=3; hand-supplied full boxes.
    b=[(e('-3/2'),e('3/2'),e('-3/2'),e('3/2')),
       (e(0,1),e(4,1),e(-2),e(2)),
       (e('-1/2',1),e('5/2'),e('-5/2'),e('5/2'))]
    r=slope_certificate(b,-10,10,-10,10)
    assert (r['beta_lower'],r['beta_upper'],r['endpoint_radius'])==(2,2,F(1,2))
    b[1]=(e('9/4',1),e('25/4',1),e(-2),e(2))
    b[2]=(e('7/4',1),e('5/2'),e('-5/2'),e('5/2'))
    r=slope_certificate(b,-10,10,-10,10)
    assert (r['beta_lower'],r['beta_upper'],r['endpoint_radius'])==(-1,-1,F(1,8))
    import random
    rng=random.Random(3901)
    iv=[(F(a),F(b)) for a in range(12) for b in range(a+1,13)]
    tree=UnionTree(iv); active=set()
    for _ in range(2000):
        k=rng.randrange(len(iv)); a,b=iv[k]
        if k in active: active.remove(k); tree.update(a,b,-1)
        else: active.add(k); tree.update(a,b,1)
        assert tree.value()==union_length(iv[j] for j in active)
    print('PASS 5 toy checks and 2000 dynamic-union updates (seed 3901); conditional sweep only')
if __name__=='__main__': selftest()
