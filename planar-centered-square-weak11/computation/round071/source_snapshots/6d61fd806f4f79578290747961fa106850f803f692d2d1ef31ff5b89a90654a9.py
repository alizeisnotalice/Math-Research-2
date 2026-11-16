#!/usr/bin/env python3
"""Exact actual-support audit and small strict deletion neighborhoods for R71.

If positive residual distance prefixes have mass M_d<4d^2, then any epsilon
<= min(d/2,(4d^2-M_d)/(16d)) keeps every subset mass strictly below the square
area at nearby receivers: max-distance decreases by at most epsilon. Combining
epsilon<=full_margin/4 with r=1/4 gives a strict full open witness of radius
r+epsilon throughout that open epsilon-square. No full halo area is computed.
"""
import hashlib,json,os,sys,time
from fractions import Fraction as F
from pathlib import Path

t=time.monotonic();src=Path("computation/results/round071-joint/one-column-result.json");d=json.loads(src.read_text())
O=[(tuple(map(F,p)),F(m)) for p,m in zip(d["primal"]["positions"],d["primal"]["weights"]) if F(m)>0]
S=[((F(6*i,25),F(6*j,25)),F(1,25)) for i in range(5) for j in range(5)];A=S+O
xs=sorted(set(p[0] for p,m in A));ys=sorted(set(p[1] for p,m in A));closures={}
for i,xl in enumerate(xs):
 for xh in xs[i:]:
  for j,yl in enumerate(ys):
   for yh in ys[j:]:
    hit=tuple(k for k,(p,m) in enumerate(A) if xl<=p[0]<=xh and yl<=p[1]<=yh)
    if not hit or hit in closures:continue
    ps=[A[k][0] for k in hit];mass=sum(A[k][1] for k in hit);span=max(max(p[0] for p in ps)-min(p[0] for p in ps),max(p[1] for p in ps)-min(p[1] for p in ps))
    assert mass<=max(F(1),span*span)
    closures[hit]=(mass,span)
eligible=[hit for hit,(m,s) in closures.items() if m>s*s];maxmass=max(closures[hit][0] for hit in eligible);assert maxmass==1
maxblocks=[hit for hit in eligible if closures[hit][0]==1]
parent=list(range(len(A)))
def root(x):
 while parent[x]!=x:x=parent[x]
 return x
for hit in eligible:
 for k in hit[1:]:parent[root(k)]=root(hit[0])
receivers=[];boxes=[];r=F(1,4)
for pp in d["centres"]:
 p=tuple(map(F,pp));ds=sorted(set(max(abs(q[0]-p[0]),abs(q[1]-p[1])) for q,m in O));prefix=[]
 for dist in ds:
  mass=sum(m for q,m in O if max(abs(q[0]-p[0]),abs(q[1]-p[1]))<=dist);slack=4*dist*dist-mass;assert slack>=0
  prefix.append({"radius":str(dist),"mass":str(mass),"slack":str(slack)})
 full_mass=sum(m for q,m in A if max(abs(q[0]-p[0]),abs(q[1]-p[1]))<=r);margin=full_mass-4*r*r;assert margin>0
 strict_open_radius=F(251,1000);assert full_mass>4*strict_open_radius**2
 active=[a for a in prefix if F(a["slack"])==0]
 eps=F(0)
 if not active:
  eps=min([F(1,4),margin/4]+[min(F(a["radius"])/2,F(a["slack"])/(16*F(a["radius"]))) for a in prefix])
  assert eps>0
  for a in prefix:assert F(a["mass"])<4*(F(a["radius"])-eps)**2
  assert full_mass>4*(r+eps)**2
  box=[p[0]-eps,p[0]+eps,p[1]-eps,p[1]+eps];boxes.append(box)
 receivers.append({"centre":[str(a) for a in p],"full_closed_radius":"1/4","full_captured_mass":str(full_mass),"full_margin":str(margin),"open_full_radius":"251/1000","open_full_margin_at_least":str(full_mass-4*strict_open_radius**2),"positive_residual_prefixes":prefix,"active_positive_residual_prefixes":active,"strict_deletion_open_box_epsilon":str(eps),"strict_deletion_open_box_bounds":[str(a) for a in box] if eps else None})
for i,b in enumerate(boxes):
 for c in boxes[i+1:]:assert b[1]<=c[0] or c[1]<=b[0] or b[3]<=c[2] or c[3]<=b[2]
area=sum((b[1]-b[0])*(b[3]-b[2]) for b in boxes);total=sum(m for p,m in A);span=max(xs[-1]-xs[0],ys[-1]-ys[0])
out={"python":sys.version,"pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"primal_source_sha256":hashlib.sha256(src.read_bytes()).hexdigest(),"S_count":25,"outside_positive_count":len(O),"actual_positive_count":len(A),"outside_total_mass":str(sum(m for p,m in O)),"actual_total_mass":str(total),"actual_bbox":[str(xs[0]),str(xs[-1]),str(ys[0]),str(ys[-1])],"actual_span":str(span),"actual_bbox_terminal_gap_d2_minus_35m_over33":str(span*span-F(35,33)*total),"complete_actual_bbox_closures":len(closures),"eligible_closures":len(eligible),"maximum_eligible_mass":str(maxmass),"maximum_blocks":[list(a) for a in maxblocks],"eligible_component_count":len(set(root(i) for i in range(len(A)))),"outside_positive_atoms":[{"position":[str(a) for a in p],"mass":str(m)} for p,m in O],"receivers":receivers,"strict_deletion_open_box_count":len(boxes),"strict_deletion_area_lower":str(area),"strict_deletion_area_lower_float":float(area),"seconds":time.monotonic()-t,"scope":"Exact support geometry, six point witnesses, five disjoint open deletion boxes only; no full halo area or general C2 bound; sixth centre is on residual boundary."}
Path("computation/results/round071-joint/primal-geometry-certificate.json").write_text(json.dumps(out,indent=2)+"\n")
Path("computation/inputs/r071/positive33-actual-input.json").write_text(json.dumps({"threshold":"1","background":"Lebesgue R2","geometry":"centered axis-parallel common-radius squares","S_indices":list(range(25)),"positions":[[str(a) for a in p] for p,m in A],"masses":[str(m) for p,m in A]},indent=2)+"\n")
print(json.dumps({k:out[k] for k in ["pid","actual_positive_count","complete_actual_bbox_closures","eligible_closures","maximum_eligible_mass","eligible_component_count","actual_bbox_terminal_gap_d2_minus_35m_over33","strict_deletion_open_box_count","strict_deletion_area_lower","strict_deletion_area_lower_float","seconds"]},indent=2))
