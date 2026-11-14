"""R16 bounded exact audit: strict last-exit examples and targeted geometry.
No random search, LP, imported capacity enumerator, or changes to R5/R9.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json,platform,sys,time
start=time.process_time(); root=Path(sys.argv[1]); output=Path(sys.argv[2])
path=root/'computation/results/mass-capacity-attempt/exact-targeted.json'; saved=json.loads(path.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def distance(p,q):return max(abs(p[k]-q[k]) for k in (0,1))
def masks_at_side(points,side):
 r=side/2; coordinates=[]
 for k in (0,1):
  vals=sorted({p[k]+s*r for p in points for s in (-1,1)})
  coordinates.append(vals+[(a+b)/2 for a,b in zip(vals,vals[1:])]+[vals[0]-1,vals[-1]+1])
 return {tuple(i for i,p in enumerate(points) if distance(p,(x,y))<=r) for x,y in product(*coordinates)}

def last_exit(points,masses,center):
 distances=sorted({distance(p,center) for p in points});closing=[]
 for idx,d in enumerate(distances):
  S=tuple(i for i,p in enumerate(points) if distance(p,center)<=d)
  mass=sum(masses[i] for i in S);q=mass/4
  following=distances[idx+1] if idx+1<len(distances) else None
  if d*d<q and (following is None or q<following*following):closing.append((q,S))
 assert closing
 return max(closing)

last_cases=[]
for name,points,masses,expect_sq,expect_label in [
 ('old_root_exactly_at_positive_entering_jump',[(F(0),F(0)),(F(1),F(0))],[F(4),F(1)],F(5,4),(0,1)),
 ('negative_gap_then_positive_again',[(F(0),F(0)),(F(2),F(0))],[F(1),F(20)],F(21,4),(0,1)),
 ('entering_jump_only_reaches_zero',[(F(0),F(0)),(F(1),F(0))],[F(1),F(3)],F(1,4),(0,)),
]:
 sq,S=last_exit(points,masses,points[0]);assert (sq,S)==(expect_sq,expect_label)
 assert all(distance(points[i],points[0])**2<sq for i in S)
 assert all(distance(points[i],points[0])**2>sq for i in range(len(points)) if i not in S)
 last_cases.append({'case':name,'points':[[str(x),str(y)] for x,y in points],'masses':list(map(str,masses)),'last_radius_squared':str(sq),'constant_label_near_last_exit':list(S)})

# Moving leftmost-anchor intervals: the same exact event comparison with ell^2.
xs=list(map(F,['0','1/5','4','401/100','20']));weights=list(map(F,['1/100','1/100','1/10000','1','1/2']))
remaining=list(range(5));parts=[]
while remaining:
 anchor=xs[remaining[0]];events=[xs[i]-anchor for i in remaining];closings=[]
 for k,d in enumerate(events):
  S=remaining[:k+1];mass=sum(weights[i] for i in S);following=events[k+1] if k+1<len(events) else None
  if d*d<mass and (following is None or mass<following*following):closings.append((mass,S[:]))
 assert closings
 side_sq,S=max(closings,key=lambda v:v[0]);assert S==remaining[:len(S)]
 assert all(xs[i]<anchor for old in parts for i in old['S'])
 parts.append({'S':S,'anchor':str(anchor),'critical_side_squared':str(side_sq)})
 remaining=remaining[len(S):]
assert [z['S'] for z in parts]==[[0],[1],[2,3],[4]]

p=[tuple(map(F,z)) for z in saved['points']];m=list(map(F,saved['masses']));M=sum(m)
assert M==F(saved['M'])==F(41,200)
diam=max(max(x[k] for x in p)-min(x[k] for x in p) for k in (0,1))
assert diam*diam==F(saved['diameter_squared'])==F(121,100)>M
caps={tuple(map(int,s.split(','))):F(c) for s,c in saved['all_capacities'].items()}
assert len(caps)==23 and all(0<c<=sum(m[i] for i in S) for S,c in caps.items())
partition=saved['exact_saturated_partition'];assert sorted(i for S in partition for i in S)==list(range(7))
assert all(caps[tuple(S)]==sum(m[i] for i in S) for S in partition)
# Exact geometric critical-radius gaps supporting the direct proof.
cluster=tuple(range(6));cluster_center=(F(0),F(0));cluster_q=sum(m[i] for i in cluster)/4
assert max(distance(p[i],cluster_center)**2 for i in cluster)<cluster_q<distance(p[6],cluster_center)**2
far_q=m[6]/4;assert far_q<min(distance(p[6],p[i])**2 for i in cluster)
assert caps[(0,)]==caps[(1,)]==F(1,25)<F(1,10)
fixed_sides={}
for side in [F(199,1000),F(1,5),F(201,1000)]:
 masks=masks_at_side(p,side)
 has=[(i,) in masks for i in (0,1)]
 assert has==([True,True] if side<F(1,5) else [False,False])
 fixed_sides[str(side)]=has
# r_n=1/10-1/(100n), n>=1, stays above 9/100>3/50;
# included internal point has distance3/50, all excluded points have >=1/10.
for i,center in [(0,(F(-2,25),F(0))),(1,(F(2,25),F(0)))]:
 assert distance(p[i],center)==F(3,50)<F(9,100)
 assert min(distance(p[j],center) for j in range(7) if j!=i)==F(1,10)
 assert 4*F(1,10)**2==F(1,25)<m[i]

four=saved['overlapping_saturated_union_failure'];p4=[tuple(map(F,z)) for z in four['points']];m4=list(map(F,four['masses']))
c4={tuple(map(int,s.split(','))):F(c) for s,c in four['all_capacities'].items()}
assert len(c4)==9
for S,center in [(tuple(four['saturated_S']),(F(47,100),F(0))),(tuple(four['saturated_T']),(F(137,100),F(0)))]:
 q=sum(m4[i] for i in S)/4
 assert max(distance(p4[i],center)**2 for i in S)<q<min(distance(p4[i],center)**2 for i in range(4) if i not in S)
 assert c4[S]==sum(m4[i] for i in S)==F(3,2)
span=max(max(z[k] for z in p4)-min(z[k] for z in p4) for k in (0,1))
assert span**2==F(four['union_diameter_squared'])==F(2116,625)>sum(m4)==F(2)
assert tuple(range(4)) not in c4
assert all(c4[(i,)]==m4[i] for i in range(4))
# Coincident labels always have identical incidence, including closed boundaries.
pdup=[(F(0),F(0)),(F(0),F(0)),(F(1),F(0))];pmerged=[(F(0),F(0)),(F(1),F(0))]
collision_checks=[]
for side in [F(1,10),F(1),F(3,2)]:
 original=masks_at_side(pdup,side);merged=masks_at_side(pmerged,side)
 assert all((0 in S)==(1 in S) for S in original)
 mapped={tuple(i for i,group in enumerate([(0,1),(2,)]) if any(j in S for j in group)) for S in original}
 assert mapped==merged
 collision_checks.append(str(side))
result={'scope':'analytic-proof audit plus bounded exact checks only; no search, LP, or full capacity-table re-enumeration',
 'python':platform.python_version(),'cpu_seconds':time.process_time()-start,'script_sha256':sha(Path(__file__)),
 'sources':{str(z):sha(z) for z in [path,root/'proofs/mass-capacity-equality-attempt.md',root/'rounds/016-mass-capacity-equality.md']},
 'last_exit_cases':last_cases,'horizontal_prefix_example':{'positions':list(map(str,xs)),'masses':list(map(str,weights)),'partition':parts},
 'seven_source':{'all_23_saved_capacities_positive_and_at_most_label_mass':True,'critical_singletons':str(F(1,25)),'near_boundary_fixed_side_exact_incidence':fixed_sides,'cluster_capacity':str(4*cluster_q),'far_singleton':str(m[6]),'partition':partition,'total':str(M)},
 'four_source':{'nine_saved_labels':'arithmetic checked; key saturated labels have direct critical-radius gap','S_and_T_capacity':'3/2','union_minimum_area':'2116/625','union_mass':'2','all_singletons_saturated':True},
 'collision_fixed_side_checks':collision_checks,'general_P_equals_M_proved':False}
output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
