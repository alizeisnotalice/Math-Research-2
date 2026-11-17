#!/usr/bin/env python3
"""Exact finite checks for the R32 transvection avoidance reduction."""
from __future__ import annotations
import hashlib,itertools,json,platform,sys,time
from pathlib import Path
OUT=Path(__file__).with_name("transvection_avoidance_result.json")
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for s in range(0,len(a),2*h):
   for j in range(h):
    x,y=a[s+j],a[s+j+h];a[s+j],a[s+j+h]=x+y,x-y
  h*=2
 return a
def dmask(points):
 out=0
 for x in points:
  for y in points:out|=1<<(x^y)
 return out
def rank(vs):
 b={}
 for v in vs:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return len(b)
def anf_degree(mask,n):
 a=[(mask>>x)&1 for x in range(1<<n)]
 for i in range(n):
  for x in range(1<<n):
   if x>>i&1:a[x]^=a[x^(1<<i)]
 return max((x.bit_count() for x,v in enumerate(a) if v),default=0)
def two_subspaces_6():
 out=set()
 for u in range(1,64):
  for v in range(u+1,64):
   if v!=u:
    out.add((1<<0)|(1<<u)|(1<<v)|(1<<(u^v)))
 return sorted(out)
def geometry_check():
 kernels=two_subspaces_6();assert len(kernels)==651
 total=yes=degree2=degree2_yes=0
 for B in itertools.combinations(range(16),8):
  total+=1;DB=dmask(B);line=any(not(DB>>v&1) for v in range(1,16))
  A=[t|(1<<1)|(y<<2) for t in (0,1) for y in B]
  DA=dmask(A);plane=any(((K&~1)&DA)==0 for K in kernels)
  assert line==plane
  yes+=line
  Bmask=sum(1<<x for x in B)
  if anf_degree(Bmask,4)<=2:
   degree2+=1;degree2_yes+=line
 return {"eight_point_B_total":total,"two_dim_subspaces_F2_6":len(kernels),"projectable_B":yes,"nonprojectable_B":total-yes,"equivalence_checks_passed":total,"degree_le_2_B_count":degree2,"degree_le_2_all_projectable":degree2_yes==degree2}
def actual_eight_variable():
 def pi(w):
  indicator=int(((w>>2)&3)==0)
  return (w&~1)|((w&1)^1^indicator)
 assert len({pi(w) for w in range(16)})==16
 def phase(z,w):return ((z&pi(w)).bit_count()&1)
 signs=[];signsT=[];diff=[]
 for z in range(16):
  for w in range(16):
   signs.append(1 if phase(z,w)==0 else -1)
   zT=z^(w&1)
   signsT.append(1 if phase(zT,w)==0 else -1)
   if phase(z,w)^phase(zT,w):diff.append(z|(w<<4))
 A=[z|(w<<4) for z in range(16) for w in range(16) if (w&1) and ((w>>2)&3)==0]
 assert sorted(diff)==sorted(A) and len(A)==32
 WA,WT=fwht(signs),fwht(signsT)
 # Reduced B uses coordinates (z1,z2,z3,w1,w2,w3).
 B=[z3|(w3<<3) for z3 in range(8) for w3 in range(8) if ((w3>>1)&3)==0]
 DB=dmask(B);DA=dmask(A)
 KB=[1<<(3+1),1<<(3+2)]
 KA=[1<<4,1<<6,1<<7]
 KBspan={0,KB[0],KB[1],KB[0]^KB[1]}
 KAspan={0}
 for v in KA:KAspan|={x^v for x in list(KAspan)}
 return {"pi_is_permutation":True,"f_bent":all(abs(x)==16 for x in WA),"fT_bent":all(abs(x)==16 for x in WT),"difference_weight":len(A),"A_formula_verified":True,"reduced_B_weight":len(B),"reduced_B_anf_degree":anf_degree(sum(1<<x for x in B),6),"KB_rank":rank(KB),"KA_rank":rank(KA),"KB_avoids_DB":all(not(DB>>x&1) for x in KBspan if x),"KA_avoids_DA":all(not(DA>>x&1) for x in KAspan if x)}
def main():
 t=time.perf_counter();out={"scope":"Exact k=3 geometry and one explicit k=4 bent transvection only.","python":sys.version,"platform":platform.platform(),"k3_geometry":geometry_check(),"k4_actual_bent_transvection":actual_eight_variable(),"elapsed_seconds":time.perf_counter()-t,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");print(json.dumps(out,ensure_ascii=False,indent=2))
 assert out["k3_geometry"]["degree_le_2_B_count"]==870 and out["k3_geometry"]["degree_le_2_all_projectable"]
 assert all(out["k4_actual_bent_transvection"].values())
if __name__=="__main__":main()
