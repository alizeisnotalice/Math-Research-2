"""R19 bounded symbolic/rational checks for mass-cell structure.
No candidate optimization, current 204-source geometry, or capacity-equality work.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,platform,sys,time
import sympy as s
start=time.process_time(); out=Path(sys.argv[1])
M,N,lam=s.symbols('M N lam',positive=True);a,b,cv=s.symbols('a b cv',real=True)
t=s.sqrt(M)/(2*s.sqrt(lam))
assert s.simplify(s.diff(t,M)*cv-cv/(8*lam*t))==0
expanded=s.expand((a+t)*(b+t))
expected=a*b+(a+b)*s.sqrt(M)/(2*s.sqrt(lam))+M/(4*lam)
assert s.simplify(expanded-expected)==0
u=s.sqrt(N)/(2*s.sqrt(lam));cross=s.expand((a+t)*(b-u))
cross_expected=a*b+b*s.sqrt(M)/(2*s.sqrt(lam))-a*s.sqrt(N)/(2*s.sqrt(lam))-s.sqrt(M*N)/(4*lam)
assert s.simplify(cross-cross_expected)==0
h=s.symbols('h',real=True);w=s.symbols('w',positive=True)
activation_right=s.limit((1+h-s.sqrt(1+h))/h,h,0,dir='+');assert activation_right==s.Rational(1,2)
activation_continuity=s.limit(w-s.sqrt(w),w,1,dir='+');assert activation_continuity==0
# Endpoint/union regressions; these are mechanisms, not complete periodic halos.
tangent=s.diff(s.sqrt(1+h),h).subs(h,0);assert tangent==s.Rational(1,2)
transverse_plus=s.limit((s.sqrt(1+h)-1)/h,h,0,dir='+')
transverse_minus=s.limit((s.sqrt(1-h)-1)/h,h,0,dir='-')
assert (transverse_plus,transverse_minus)==(s.Rational(1,2),-s.Rational(1,2))
cap_lam=F(3,2);cap_R=F(1,2);cap_mass=4*cap_lam*cap_R**2
assert cap_mass==F(3,2) and cap_mass/cap_lam==4*cap_R**2
cap_left=1/cap_lam;cap_right=F(0)
# General normalized-objective derivative, via an actual symbolic quotient.
z=s.symbols('z',real=True);A0,Ad,T0,Td=s.symbols('A0 Ad T0 Td',real=True)
J=lam*(A0+z*Ad)/(T0+z*Td)
Jprime=s.simplify(s.diff(J,z).subs(z,0));assert s.simplify(Jprime-lam*(T0*Ad-A0*Td)/T0**2)==0
# Concrete counterexample to area-only improvement on nonconserved families.
# One point per 2x2 cell, lambda=1, w in [1/4,1]. For side L>=2,
# count/L^2 <= w*(1/2+1/L)^2 <= w<=1; hence only singleton halo contributes.
# Its area is w, while lambda*A/w is identically1. R=2 gives uniform G<=9/16.
normalization=[]
for mass in [F(1,4),F(1)]:
 root=F(1,2) if mass==F(1,4) else F(1)
 area=4*(root/2)**2;ratio=area/mass
 assert area==mass and ratio==1
 normalization.append({'mass':str(mass),'area':str(area),'weak_ratio':str(ratio)})
assert F(1)/(2*2)*(1+F(2,4))**2==F(9,16)<1

# Finite-box membership against direct exact radius-event witnesses.
P,Q,R,level=F(4),F(3),F(1),F(1)
base=[(F(0),F(0)),(F(1,2),F(1,4))]
def ceil(q):return -((-q).__floor__())
copies=[]
for i,(x,y) in enumerate(base):
 for j in range(ceil((-R-x)/P),((P+R-x)/P).__floor__()+1):
  for k in range(ceil((-R-y)/Q),((Q+R-y)/Q).__floor__()+1):copies.append((x+j*P,y+k*Q,i))
xs=sorted({p[0] for p in copies});ys=sorted({p[1] for p in copies});boxes=[]
for a,b in product(xs,xs):
 if a>b or b-a>2*R:continue
 for c,d in product(ys,ys):
  if c>d or d-c>2*R:continue
  captured=[p for p in copies if a<=p[0]<=b and c<=p[1]<=d]
  if captured:boxes.append((a,b,c,d,[sum(p[2]==i for p in captured) for i in range(len(base))]))
assert boxes and all(sum(counts)>0 for *_,counts in boxes)
centers=list(product([F(0),F(1,8),F(1,4),F(1,2),F(1),F(31,8)],
                     [F(0),F(1,8),F(1,4),F(1,2),F(23,8)]))
centers +=[(F(0),F(249,1000)),(F(0),F(251,1000))]
checks=0;closed_boundary={}
for weights in [[F(1,16),F(1,4)],[F(1,4),F(1,4)],[F(1,8),F(1,8)]]:
 # All weights in compact [1/16,1/4]^2 have G_0(1)<=5/16.
 assert sum(weights)/(P*Q)*(1+P/(2*R))*(1+Q/(2*R))<=F(5,16)<level
 for x,y in centers:
  distances=sorted({max(abs(x-a),abs(y-b)) for a,b,i in copies})
  direct=any(d<R and sum(weights[i] for a,b,i in copies if max(abs(x-a),abs(y-b))<=d)>4*level*d*d for d in distances)
  union=False
  for a,b,c,d,counts in boxes:
   mass=sum(F(n)*v for n,v in zip(counts,weights));assert mass>=F(1,16)
   distance=max(abs(x-a),abs(x-b),abs(y-c),abs(y-d))
   if distance<R and 4*level*distance**2<mass:union=True;break
  assert direct==union;checks+=1
  if weights==[F(1,4),F(1,4)] and x==0 and y in [F(249,1000),F(1,4),F(251,1000)]:closed_boundary[str(y)]=direct
assert closed_boundary=={'1/4':False,'249/1000':True,'251/1000':False}
result={'scope':'R19 symbolic formulas and bounded rational membership only; no project candidate cell enumeration',
 'python':platform.python_version(),'sympy':s.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'dt_formula_residual':'0','diagonal_sqrt_product_includes_linear_mass':str(expanded),
 'cross_sqrt_product_expansion':str(cross),'activation_one_sided':['0',str(activation_right)],
 'persistent_wall_ray_derivative':str(tangent),'transverse_union_wall_one_sided':[str(transverse_minus),str(transverse_plus)],
 'cap_mechanism':{'lambda':str(cap_lam),'R':str(cap_R),'mass_at_wall':str(cap_mass),'left_area_derivative':str(cap_left),'right_area_derivative':str(cap_right),'scope':'endpoint min mechanism; not a family claimed to satisfy uniform tail on both sides'},
 'general_weak_ratio_derivative':str(Jprime),'nonconserved_area_only_counterexample':normalization,'counterexample_uniform_tail_at_R2':'9/16',
 'tiny_periodic_box_coverage':{'periods':[str(P),str(Q)],'R':str(R),'copies':len(copies),'all_nonempty_coordinate_boxes':len(boxes),'exact_center_membership_comparisons':checks,'strict_closed_boundary_centers':closed_boundary,'uniform_tail_for_full_mass_polytope':'<=5/16<1'},
 'cpu_seconds':time.process_time()-start}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
