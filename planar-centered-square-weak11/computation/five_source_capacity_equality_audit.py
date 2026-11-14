"""Independent bounded R21 exact checks; no search, LP, or coefficient imports.
The all-coordinate proof is reviewed analytically; finite cases exercise shared
coordinate walls, unequal diagonal offsets, strict endpoints and saved fixtures.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,platform,sys,time
start=time.process_time();root=Path(sys.argv[1]);out=Path(sys.argv[2])
proof=root/'proofs/five-source-capacity-equality.md';source=root/'computation/results/five-source-capacity/exact-checks.json';saved=json.loads(source.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def incidence(points,box):
 a,b,c,d=box;return tuple(i for i,(x,y) in enumerate(points) if a<=x<=b and c<=y<=d)
def full_quadrants(p,others):return all(any(s*(q[0]-p[0])>=0 and t*(q[1]-p[1])>=0 for q in others) for s,t in product([-1,1],repeat=2))
def masks_at_side(points,ell):
 radius=ell/2;coords=[]
 for k in (0,1):
  vals=sorted({p[k]+sgn*radius for p in points for sgn in [-1,1]})
  coords.append(vals+[(a+b)/2 for a,b in zip(vals,vals[1:])]+[vals[0]-1,vals[-1]+1])
 return {tuple(i for i,(a,b) in enumerate(points) if abs(x-a)<=radius and abs(y-b)<=radius) for x,y in product(*coords)}
# Unequal a,b ensure no nonisometric normalization was used; all t0/v boundaries.
a,b,s0,r0=F(2,3),F(3,5),F(1,7),F(2,9)
diagonal_checks=[]
for t0,v in product([F(0),F(1,8)],repeat=2):
 points=[(F(0),F(0)),(a,b),(-s0,F(0)),(a,b+r0),(a+t0,-v)]
 assert full_quadrants(points[0],points[1:]) and full_quadrants(points[1],[points[0],*points[2:]])
 Dx,Dy=a+s0+t0,b+r0+v
 for theta in [F(1,100),F(1,2),F(99,100)]:
  ell=theta*Dx
  low,high=max(-s0,-ell),min(F(0),a+t0-ell);assert low<high
  left=(low+high)/2;top=min(b,ell)/2
  bp=(left,left+ell,top-ell,top)
  assert incidence(points,bp)==(0,) and bp[1]-bp[0]==bp[3]-bp[2]
  ell=theta*Dy
  low,high=max(-v,b-ell),min(b,b+r0-ell);assert low<high
  bottom=(low+high)/2
  lowx,highx=max(F(0),a-ell),a;assert lowx<highx
  left=(lowx+highx)/2
  bq=(left,left+ell,bottom,bottom+ell)
  assert incidence(points,bq)==(1,) and bq[1]-bq[0]==bq[3]-bq[2]
 diagonal_checks.append({'a':str(a),'b':str(b),'s0':str(s0),'r0':str(r0),'t0':str(t0),'v':str(v),'Dx_squared':str(Dx*Dx),'Dy_squared':str(Dy*Dy),'isolating_square_checks':6})
# Shared-coordinate construction, including guards exactly on x=a. The q mass
# has its first root exactly at g+ entry: last-exit must continue past the jump.
aligned_checks=[]
for xplus,xminus in product([a,a+F(1,7)],[a,a+F(2,9)]):
 sp=F(1,5);yp,ym=F(3,8),-F(4,9)
 points=[(F(0),F(0)),(a,F(0)),(-sp,F(0)),(xplus,yp),(xminus,ym)]
 masses=[F(1,20),F(9,16),F(1,100),F(1,100),F(1,100)]
 assert full_quadrants(points[0],points[1:]) and full_quadrants(points[1],[points[0],*points[2:]])
 pair_mass=masses[0]+masses[2];assert pair_mass>sp*sp
 entries={1:F(0),3:max(xplus-a,2*yp),4:max(xminus-a,-2*ym)}
 closing=[]
 for entry in sorted(set(entries.values())):
  S=tuple(i for i,e in entries.items() if e<=entry);mass=sum(masses[i] for i in S)
  later=[e for e in entries.values() if e>entry];next_entry=min(later) if later else None
  if entry*entry<mass and (next_entry is None or mass<next_entry*next_entry):closing.append((mass,S))
 final_mass,T=max(closing);assert T==(1,3) and final_mass==F(229,400)
 assert masses[1]==entries[3]**2 and max(entries[i]**2 for i in T)<final_mass<entries[4]**2
 partition=[(0,2),T,(4,)];assert sorted(i for S in partition for i in S)==list(range(5))
 D=max(max(p[k] for p in points)-min(p[k] for p in points) for k in (0,1));assert sum(masses)<D*D
 aligned_checks.append({'xplus':str(xplus),'xminus':str(xminus),'pair_critical_side_squared':str(pair_mass),'q_right_halfplane_last_side_squared':str(final_mass),'partition':partition,'entry_root_coincidence_tested':True,'scope':'checks the conditional geometric subconstruction; not a claim these chosen p,q are both anomalous'})
# Saved examples: direct strict critical-radius gaps certify the partition;
# singleton upper bounds are reviewed analytically and checked at closed walls.
fixtures={}
for name,record in saved['exact_examples'].items():
 points=[tuple(map(F,p)) for p in record['points']];masses=list(map(F,record['masses']));caps={tuple(map(int,S.split(','))):F(c) for S,c in record['all_capacities'].items()}
 assert all(0<c<=sum(masses[i] for i in S) for S,c in caps.items())
 parts=record['saturated_partition'];assert sorted(i for S in parts for i in S)==list(range(5))
 assert all(caps[tuple(S)]==sum(masses[i] for i in S) for S in parts)
 boundary={}
 threshold=F(2) if name.startswith('axis') else F(3)
 for ell in [threshold-F(1,100),threshold,threshold+F(1,100)]:
  value=(0,) in masks_at_side(points,ell);assert value==(ell<threshold);boundary['p_'+str(ell)]=value
 threshold=F(4) if name.startswith('axis') else F(3)
 for ell in [threshold-F(1,100),threshold,threshold+F(1,100)]:
  value=(1,) in masks_at_side(points,ell);assert value==(ell<threshold);boundary['q_'+str(ell)]=value
 if name.startswith('axis'):
  assert caps[(0,)]==4<5 and caps[(1,)]==16<17
  assert F(record['M'])==F(223,10)<F(record['diameter_squared'])==25
  assert caps[(0,2)]==F(51,10)>1 and caps[(1,3,4)]==F(86,5)>9
 else:
  assert caps[(0,)]==caps[(1,)]==9<10
  assert caps[tuple(range(5))]==F(203,10)>9
 fixtures[name]={'M':record['M'],'diameter_squared':record['diameter_squared'],'partition':parts,'critical_singletons':[str(caps[(0,)]),str(caps[(1,)])],'fixed_side_closed_incidence':boundary,'all_saved_capacities_arithmetic_checked':len(caps),'full_continuous_capacity_table_reenumerated':False}
# Exact odd balanced-matrix obstruction: three explicit square families.
tri=[(F(0),F(0)),(F(2),F(0)),(F(1),F(2))];odd=[]
for ell in [F(21,10),F(12,5),F(17,7)]:
 assert 2<ell and ell*ell<6
 boxes=[((2-ell)/2,(2+ell)/2,-ell,F(0)),(1-ell,F(1),F(0),ell),(F(1),1+ell,F(0),ell)]
 labels=[incidence(tri,B) for B in boxes];assert labels==[(0,1),(0,2),(1,2)]
 odd.append({'side':str(ell),'incidences':labels})
matrix=[[int(i in S) for S in [(0,1),(0,2),(1,2)]] for i in range(3)]
assert all(sum(row)==2 for row in matrix) and all(sum(matrix[i][j] for i in range(3))==2 for j in range(3))
assert saved['unbalanced_hypergraph_example']['capacities']==['6','6','6']
result={'scope':'independent bounded exact examples supporting a separately reviewed analytic N<=5 theorem; no search, LP, or full capacity-table recomputation',
 'python':platform.python_version(),'cpu_seconds':time.process_time()-start,'script_sha256':sha(Path(__file__)),
 'sources':{str(p):sha(p) for p in [proof,source,root/'computation/five_source_capacity_check.py',root/'rounds/021.md']},
 'diagonal_unequal_offset_and_axis_boundary_checks':diagonal_checks,'aligned_halfplane_checks':aligned_checks,
 'saved_fixtures':fixtures,'odd_cycle_explicit_squares':odd,'odd_cycle_incidence_matrix':matrix,
 'all_N_capacity_equality_proved':False,'actual_strict_B_equals_M_claimed':False}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
