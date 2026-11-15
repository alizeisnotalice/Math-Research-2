#!/usr/bin/env python3
"""Exact small algebra and exclusion diagnostics; no production halo kernel."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

start = time.process_time()
def st(x): return str(F(x))
def dot(a, b): return sum((x*y for x,y in zip(a,b)),F(0))
def affine_sum(weights, vectors):
    return [sum((a*v[j] for a,v in zip(weights,vectors)),F(0)) for j in range(len(vectors[0]))]
def f(x,y):
    return F(1,3)+dot([F(2),F(-1)],x)+dot([F(-2),F(3)],y)+x[0]*(2*y[0]-y[1])+x[1]*(y[0]+3*y[1])
# P is a one-dimensional segment in R^2; Q is a two-dimensional triangle.
P=[[F(-1),F(1)],[F(2),F(-2)]]
Q=[[F(0),F(0)],[F(1),F(0)],[F(0),F(2)]]
alpha=[F(2,5),F(3,5)];beta=[F(1,6),F(1,3),F(1,2)]
x,y=affine_sum(alpha,P),affine_sum(beta,Q)
vertex_values=[[f(u,v) for v in Q] for u in P]
convex_sum=sum((alpha[i]*beta[j]*vertex_values[i][j] for i in range(2) for j in range(3)),F(0))
assert f(x,y)==convex_sum<=max(z for row in vertex_values for z in row)
product={'P_vertices':[[st(z) for z in r] for r in P],
         'Q_vertices':[[st(z) for z in r] for r in Q],
         'interior_value':st(f(x,y)),'vertex_convex_sum':st(convex_sum),
         'maximum_vertex_pair':st(max(z for row in vertex_values for z in row)),
         'coupled_nonproduct_counterexample':'f(s,-s)=1-s^2; maximum1 at s0, endpoint values0'}

margins=[]
for lam,m,d0,u,H in [(F(1),F(1),F(3,2),F(-1),F(1)),
                     (F(3,2),F(1,10),F(1,2),F(1),F(1,8)),
                     (F(2),F(1,4),F(2),F(1),F(3,2))]:
    q=lam*d0*d0-m;assert q>0 and u!=0 and H>0
    L=lam*abs(u)*(2*abs(d0)+abs(u)*H)
    delta=min(H,q/(2*L))
    candidates=[-delta,delta]
    stationary=-d0/u
    if -delta<=stationary<=delta:candidates.append(stationary)
    minimum=min(lam*(d0+u*h)**2-m for h in candidates)
    assert minimum>=q/2>0
    margins.append({k:st(v) for k,v in [('lambda',lam),('m',m),('d0',d0),('u',u),('H',H),('q',q),('L',L),('delta',delta),('exact_min_residual',minimum),('certified_q_half',q/2)]})
# q=0 and nonzero span speed: one side becomes nonempty in the unclipped x test.
h=F(1,10)
assert (1-h)**2-1<0<(1+h)**2-1

# Source/window prefilter omitted at s0 can become relevant without span events.
# One fixed label in the slab I: z(s)=-11/10+s, t=1/2, R=1, F=[0,10].
R=F(1);z0=F(-11,10);t=F(1,2)
assert z0+R<0 and z0+t<0
s_later=F(4,5);z_later=z0+s_later
assert z_later+t>0 and F(0)<F(2)*t
window={'source':'z(s)=-11/10+s','R':'1','radius':'1/2','F':['0','10'],
        'span_for_every_s':'0','R_window_entry_event':'1/10',
        'actual_receiver_entry_event':'3/5','later_s':st(s_later),
        'later_clipped_receiver_length':st(z_later+t)}

# Sorting only currently clipped endpoints misses the raw endpoint's F crossing.
eps=F(1,10)
def clipped_length(s): return max(F(0),min(F(10),2-eps+s)-max(F(0),-eps+s))
base=clipped_length(F(0));later=F(1,5)
assert base==F(19,10) and clipped_length(later)==2
assert base+later!=clipped_length(later)
clip={'raw_left':'-1/10+s','raw_right':'19/10+s','F':['0','10'],
      'raw_left_crossing_event':'1/10','area_at_s0':st(base),
      'left_branch_slope_at_s0':'1','later_s':st(later),
      'actual_later_length':st(clipped_length(later)),
      'false_unchanged_branch_prediction':st(base+later)}

def union_length(intervals):
    intervals=sorted((a,b) for a,b in intervals if a<b)
    if not intervals:return F(0)
    total=F(0);a,b=intervals[0]
    for c,d in intervals[1:]:
        if b<c:total+=b-a;a,b=c,d
        else:b=max(b,d)
    return total+b-a

def ceil(q):return -((-q).__floor__())
def dyadic_union_bounds(intervals,bits):
    D=1<<bits
    inner=[(F(ceil(a*D),D),F((b*D).__floor__(),D)) for a,b in intervals]
    outer=[(F((a*D).__floor__(),D),F(ceil(b*D),D)) for a,b in intervals]
    lo,hi=union_length(inner),union_length(outer)
    true=union_length(intervals);assert lo<=true<=hi
    return lo,true,hi

sets=[[(F(0),F(5,3)),(F(3,2),F(3))],
      [(F(0),F(3,5)),(F(1),F(7,4))],
      [(F(1,9),F(2,9))]]
coeff=[1,-1,0]
signed=[]
for bits in [2,4,8]:
    bounds=[dyadic_union_bounds(intervals,bits) for intervals in sets]
    true=sum((c*tr for c,(_,tr,_) in zip(coeff,bounds)),F(0))
    lower=sum((c*(lo if c>=0 else hi) for c,(lo,_,hi) in zip(coeff,bounds)),F(0))
    upper=sum((c*(hi if c>=0 else lo) for c,(lo,_,hi) in zip(coeff,bounds)),F(0))
    assert lower<=true<=upper and lower>0
    signed.append({'bits':bits,'coefficients':coeff,
                   'Y_bounds':[[st(z) for z in r] for r in bounds],
                   'beta_lower':st(lower),'beta_exact':st(true),'beta_upper':st(upper)})
# Explicit wrong-sign accumulation can exclude the exact beta.
assert F(2)-F(3)==-1
wrong_lower=F(2)-F(11,4);wrong_upper=F(9,4)-F(3)
assert wrong_lower==wrong_upper==F(-3,4) and not wrong_lower<=F(-1)<=wrong_upper

result={'status':'passed','scope':'rational small algebra and omission diagnostics only; no real halo or production kernel',
        'product_vertex_identity':product,'q_margin_cases':margins,
        'q_zero_nonpersistent_wall':{'negative_residual':st((1-h)**2-1),'positive_residual':st((1+h)**2-1)},
        'missing_window_prefilter_event':window,'missing_raw_clip_event':clip,
        'signed_dyadic_slope_bounds':signed,
        'wrong_negative_coefficient_example':{'exact_beta':'-1','wrong_interval':['-3/4','-3/4'],'rejects_exact_value':True},
        'evidence_boundary':'General claims use the accompanying analytic proof; finite toy checks diagnose signs, omitted events and rounding directions.',
        'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-start,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'actual_command':[sys.executable,str(Path(__file__).resolve()),str(Path(sys.argv[1]).resolve())]}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','scope','cpu_seconds','script_sha256']},indent=2))
