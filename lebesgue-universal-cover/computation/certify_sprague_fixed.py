#!/usr/bin/env python3
"""Proof-oriented interval checks for Gibbs's first three cuts at 1.5494 degrees.

This is deliberately not an area calculation.  It certifies the finite
inequalities used in audit-sprague/proof.md: face membership, the selected
circle roots, the endpoint-Voronoi tests, and the ``outside one disk'' fan
tests.  Every comparison is between rational interval endpoints.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from sigma_parameter import SIGMA_OUT, radians
import hashlib, json, sys
from rational_trig import PI, trig, scale

Q = 10**24

class I:
    def __init__(self, a, b=None):
        if isinstance(a, I): self.lo, self.hi = a.lo, a.hi; return
        self.lo, self.hi = F(a), F(a if b is None else b)
        assert self.lo <= self.hi
    def rnd(self): return I(F((self.lo*Q).__floor__(), Q), F((self.hi*Q).__ceil__(), Q))
    def __add__(self, b):
        b=I(b); return I(self.lo+b.lo, self.hi+b.hi).rnd()
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,b): return self+-I(b)
    def __rsub__(self,b): return I(b)+-self
    def __mul__(self,b):
        b=I(b); v=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(v),max(v)).rnd()
    __rmul__=__mul__
    def __truediv__(self,b):
        b=I(b); assert b.lo*b.hi>0
        return self*I(1/b.hi,1/b.lo).rnd()
    def sqrt(self):
        assert self.lo>=0
        lo=isqrt((self.lo*Q*Q).__floor__()); hi=isqrt((self.hi*Q*Q).__floor__())+1
        assert F(lo*lo,Q*Q)<=self.lo and F(hi*hi,Q*Q)>=self.hi
        return I(F(lo,Q),F(hi,Q))
    def out(self): return [str(self.lo),str(self.hi)]

def V(x,y): return I(x),I(y)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def mul(t,a): return tuple(t*x for x in a)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),I(0))
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def norm2(a): return dot(a,a)
def perp(a): return -a[1],a[0]
def unit(q, coefficient=0):
    a=radians(q, coefficient); return I(*trig(a,'cos')), I(*trig(a,'sin'))
# Cramer's-rule formula for a.x=b.x=1/2; no floating operations occur.
def meet(a,b):
    d=cross(a,b)
    return (I(F(1,2))*(b[1]-a[1])/d, I(F(1,2))*(a[0]-b[0])/d)
def project(c,n): return sub(c,mul(dot(n,c)-I(F(1,2)),n))
def cl(c,n,sign):
    d=dot(n,c)-I(F(1,2))
    return add(project(c,n),mul(I(sign)*(I(1)-d*d).sqrt(),perp(n)))
def cc(a,b,sign):
    q=sub(b,a); ds=norm2(q)
    return add(mul(I(F(1,2)),add(a,b)),mul(I(sign)*(I(1)-ds/I(4)).sqrt()/ds.sqrt(),perp(q)))
def param(x,a,b): return dot(sub(x,a),sub(b,a))/norm2(sub(b,a))

n=[unit(F(j,3)) for j in range(6)]
m=[unit(F(1,6)+F(j,3), 1) for j in range(6)]
A1=meet(n[1],n[2]); B1=meet(n[0],n[1]); D1=meet(n[4],n[5]); E1=meet(n[3],n[4]); F1=meet(n[2],n[3])
B3=meet(n[1],m[0]); B2=meet(n[0],m[0]); C3=meet(n[0],m[5]); C2=meet(n[5],m[5])
D3=meet(n[5],m[4]); E3=meet(n[4],m[3]); E2=meet(n[3],m[3]); F3=meet(n[3],m[2]); F2=meet(n[2],m[2])
K=project(F3,m[5]); G=cl(F3,n[5],1)
Ip=project(B3,m[3]); H=project(C3,n[3]); J=cc(B3,C3,-1); Jother=cc(B3,C3,1)
P=project(G,n[2]); Qp=project(E3,n[1]); X=cc(G,E3,-1); Xother=cc(G,E3,1)

# Chosen roots and all four tangent feet are on the source's stated segments.
for name,x,a,b in [('K',K,C3,C2),('G',G,D3,C2),('I',Ip,E3,E2),('H',H,E2,F3),
                   ('P',P,F1,A1),('Q',Qp,A1,B1)]:
    t=param(x,a,b); assert t.lo>=0 and t.hi<=1, (name,t.out())
assert Jother[0].lo>I(F(1,2)).hi       # the other J root is outside P
assert Xother[1].hi<D1[1].lo            # the other X root is outside P
P_slacks={}
for name,x in [('J',J),('X',X)]:
    vals=[]
    for tag,u in [(f'n{j}',n[j]) for j in range(6)]+[('m3',m[3]),('m5',m[5])]:
        slack=I(F(1,2))-dot(u,x); assert slack.lo>0, (name,tag,slack.out()); vals.append((tag,slack.out()))
    P_slacks[name]=vals
for c,x in [(F3,K),(F3,G),(B3,Ip),(B3,J),(C3,J),(C3,H),(G,P),(G,X),(E3,Qp),(E3,X)]:
    e=norm2(sub(x,c))-I(1); assert e.lo<=0<=e.hi, (c,x,e.out())

# Contact lemma for the two *opposite* internal cross-sections.  P has the
# C/E cuts m5.x<=1/2 and m3.x<=1/2.  Constant width gives h(m2)>=1/2 and
# h(m0)>=1/2, since m2=-m5 and m0=-m3.  A forced original-hexagon contact
# on B1C3 / E1F1 lies strictly on the other side, so convexity crosses the
# two displayed sections.  H intersected with those sections is F3F2/B3B2.
for tag,normal,a,b in [('CB_below_m2',m[2],B1,C3),('EF_below_m0',m[0],E1,F1)]:
    assert (dot(normal,a)-I(F(1,2))).hi<0 and (dot(normal,b)-I(F(1,2))).hi<0, tag

# An endpoint is nearest on [a,b] throughout the indicated Voronoi halfspace.
def vor(z,a,b): return dot(sub(z,a),sub(b,a))
def strict_nonpositive(x, tag):
    assert x.hi<0, (tag,x.out()); return x
def nonnegative(x, tag):
    assert x.lo>0, (tag,x.out()); return x
def arc_halfplane(name,c,u0,u1,a,b,tangent_index,cells=256):
    """Certify L(c+minor_arc(u0,u1))<=0 without point sampling.

    The radial projection of each rational t-cell's chord covers the complete
    minor arc.  A cell containing a known tangent equality is handled by the
    exact derivative of L(c+w/|w|), directed into the negative half-plane.
    """
    v=sub(b,a); f0=vor(add(c,u0),a,b); f1=vor(add(c,u1),a,b)
    fs=[f0,f1]; tangent=fs[tangent_index]; other=fs[1-tangent_index]
    assert tangent.lo<=0<=tangent.hi, (name,'specified tangent endpoint',tangent.out())
    assert other.hi<0, (name,'non-tangent endpoint',other.out())
    dw=sub(u1,u0); deriv=[]
    for j in range(cells):
        t=I(F(j,cells),F(j+1,cells)); w=add(u0,mul(t,sub(u1,u0))); r2=norm2(w)
        # derivative numerator; r^3 is positive and does not alter its sign
        num=dot(v,dw)*r2-dot(v,w)*dot(w,dw)
        if j==0 and tangent_index==0:
            assert num.hi<0, (name,'inward derivative at u0',num.out()); deriv.append(num.out()); continue
        if j==cells-1 and tangent_index==1:
            assert num.lo>0, (name,'inward derivative at u1',num.out()); deriv.append(num.out()); continue
        u=mul(I(1)/r2.sqrt(),w); f=vor(add(c,u),a,b)
        assert f.hi<0, (name,j,f.out())
    return {'endpoint_values':[f0.out(),f1.out()],'tangent_derivative_numerators':deriv}

# For each circular arc its endpoint vectors have positive dot product, hence
# it is the minor arc.  The proof file uses an exhaustive radial-chord
# interval cover and a one-sided exact derivative at a tangent endpoint.
# Equalities named below are exact from the projection definitions (opposite
# normals), not rounded numerical tests.
arcs=[
 ('CS',F3,K,G,C2, [(F3,F2)],0),
 ('ES-B',B3,Ip,J,E2, [(B3,B2)],0),
 ('ES-C',C3,J,H,E2, [(C3,B1)],1),
 ('AS-E',E3,Qp,X,A1, [(E3,D1)],0),
 ('AS-G',G,X,P,A1, [(G,D1)],1),
]
arc_rows=[]
fan_boxes={}
def arc_fan_box(c,u0,u1,apex,cells=128):
    """Outward interval enclosure of the whole minor-arc fan.

    Radial projection of (1-t)u0+t*u1 parametrizes the minor arc because
    u0.u1>0.  Cells cover all real t in [0,1], so this is an interval cover,
    not a collection of point samples.
    """
    lo=[min(apex[k].lo, u0[k].lo+c[k].lo, u1[k].lo+c[k].lo) for k in range(2)]
    hi=[max(apex[k].hi, u0[k].hi+c[k].hi, u1[k].hi+c[k].hi) for k in range(2)]
    for j in range(cells):
        t=I(F(j,cells),F(j+1,cells)); w=add(u0,mul(t,sub(u1,u0)))
        u=mul(I(1)/norm2(w).sqrt(),w); z=add(c,u)
        # A fan point is a convex combination of apex and z.
        for k in range(2):
            lo[k]=min(lo[k],z[k].lo,apex[k].lo); hi[k]=max(hi[k],z[k].hi,apex[k].hi)
    return I(lo[0],hi[0]),I(lo[1],hi[1])
for name,c,u0,u1,apex,segs,tangent_index in arcs:
    v0,v1=sub(u0,c),sub(u1,c)
    nonnegative(dot(v0,v1),name+' minor-arc dot')
    # Apex lies strictly outside this disk.  Tangency accounts for a zero
    # endpoint; the other endpoint gives the strict interval slack recorded.
    d=sub(apex,c); nonnegative(norm2(d)-I(1),name+' apex exterior')
    radial0=dot(d,v0)-I(1); radial1=dot(d,v1)-I(1)
    radial=[radial0,radial1]
    assert radial[tangent_index].lo<=0<=radial[tangent_index].hi, (name,'specified disk tangent',radial[tangent_index].out())
    assert radial[1-tangent_index].lo>0, (name,'strict non-tangent disk margin',radial[1-tangent_index].out())
    vals=[]
    for a,b in segs:
        wA=vor(apex,a,b); w0=vor(add(c,u0),a,b); w1=vor(add(c,u1),a,b)
        assert wA.hi<=0, (name,'apex vor',wA.out())
        # A zero endpoint is one of I,H,P,Q and follows exactly from its
        # orthogonal-projection definition; all other entries are strict.
        assert [w0, w1][1-tangent_index].hi<0, (name, "non-tangent Voronoi endpoint", w0.out(), w1.out())
        vals.append({'apex':wA.out(),'arc':arc_halfplane(name,c,v0,v1,a,b,tangent_index)})
    arc_rows.append({'arc':name,'minor_dot_lower':str(dot(v0,v1).lo),
                     'apex_disk_margin_lower':str((norm2(d)-I(1)).lo),
                     'tangent_endpoint_index':tangent_index,
                     'radial_endpoint_margins':[radial0.out(),radial1.out()],
                     'voronoi_endpoint_values':vals})
    fan_boxes[name]=arc_fan_box(c,v0,v1,apex)

# ES/AS each split along its apex-to-crossing diagonal.  The strict oriented
# ray tests put that diagonal between the two outer rays through an angle <pi,
# hence the two visible fan interiors do not overlap except on the diagonal.
def cross_row(apex,a,b): return cross(sub(a,apex),sub(b,apex))
fan_order={}
for tag,apex,left,mid,right,sgn in [('ES',E2,Ip,J,H,1),('AS',A1,Qp,X,P,-1)]:
    c1=cross_row(apex,left,mid); c2=cross_row(apex,mid,right); c3=cross_row(apex,left,right)
    if sgn>0: assert c1.lo>0 and c2.lo>0 and c3.lo>0, (tag,c1.out(),c2.out(),c3.out())
    else: assert c1.hi<0 and c2.hi<0 and c3.hi<0, (tag,c1.out(),c2.out(),c3.out())
    fan_order[tag]={'left_mid':c1.out(),'mid_right':c2.out(),'left_right':c3.out(),'orientation':'positive' if sgn>0 else 'negative'}
# The three closed cuts are separated at this fixed angle.
cs=fan_boxes['CS']; es=(I(min(fan_boxes['ES-B'][0].lo,fan_boxes['ES-C'][0].lo),max(fan_boxes['ES-B'][0].hi,fan_boxes['ES-C'][0].hi)),I(min(fan_boxes['ES-B'][1].lo,fan_boxes['ES-C'][1].lo),max(fan_boxes['ES-B'][1].hi,fan_boxes['ES-C'][1].hi)))
ass=(I(min(fan_boxes['AS-E'][0].lo,fan_boxes['AS-G'][0].lo),max(fan_boxes['AS-E'][0].hi,fan_boxes['AS-G'][0].hi)),I(min(fan_boxes['AS-E'][1].lo,fan_boxes['AS-G'][1].lo),max(fan_boxes['AS-E'][1].hi,fan_boxes['AS-G'][1].hi)))
assert cs[0].lo>I(F(1,4)).hi and es[0].hi<I(-F(1,4)).lo
assert ass[1].lo>I(F(1,4)).hi and cs[1].hi<I(0).lo and es[1].hi<I(0).lo

# The extra DC-face strengthening required before the AS contact claim.
vG,vC=sub(G,F3),sub(C2,F3); dC=sub(C2,F3)
nonnegative(dot(vG,vC)-I(1), 'GC2 radial chord exterior margin')
nonnegative(norm2(dC)-I(1), 'C2 exterior F3 disk')
strict_nonpositive(vor(G,F3,F2), 'G nearest F3 on F3F2')
strict_nonpositive(vor(C2,F3,F2), 'C2 nearest F3 on F3F2')

out={'status':'passed','scope':'fixed-angle coordinate/branch/distance certificate, including strict signs for the support/convexity F3F2 and B3B2 contact lemma.  Together with proof.md it proves Sprague Propositions 3--5/Lemma 1 at this angle.  It does not prove Gibbs Proposition 2 or later reflection cuts.',
 'python':sys.version,"sigma_degrees_interval": SIGMA_OUT,'rounding_denominator':str(Q),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'trig_helper_sha256':hashlib.sha256(Path(__file__).with_name('rational_trig.py').read_bytes()).hexdigest(),
 'coordinates':{k:[x.out() for x in p] for k,p in [('F3',F3),('F2',F2),('B3',B3),('B2',B2),('C3',C3),('C2',C2),('E3',E3),('E2',E2),('G',G),('K',K),('I',Ip),('H',H),('J',J),('P',P),('Q',Qp),('X',X)]},
 'selected_root_P_strict_slacks':P_slacks,
 'arcs':arc_rows,
 'fan_diagonal_order':fan_order,
 'contact_section_margins':{'CB_below_m2':[ (dot(m[2],z)-I(F(1,2))).out() for z in (B1,C3)],'EF_below_m0':[ (dot(m[0],z)-I(F(1,2))).out() for z in (E1,F1)]},
 'fan_interval_boxes':{k:[x.out() for x in b] for k,b in [('CS',cs),('ES',es),('AS',ass)]},
 'dc_open_part_data':{'C2_minus_F3_squared_minus_1':(norm2(dC)-I(1)).out(),
                      'G_minus_F3_dot_C2_minus_F3_minus_1':(dot(vG,vC)-I(1)).out(),
                      'Vor_F3_values_at_G_C2':[vor(G,F3,F2).out(),vor(C2,F3,F2).out()]}}
print(json.dumps(out,indent=2))
