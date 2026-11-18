"""B1(.430,.0175): certify 11 continuous dual exclusions and 3 all-pose witnesses."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb, arb_mat, ctx

PINS={
 'certify_ah_released_tangent.py':'f43d04412b1d234095f792aa87890ce47f64242f695a17905a3d4493af349b5b',
 'certify_ah_two_line_cap.py':'3d3781ceccaf14b9db4c0bb096e6d24a98c825c0a5f1eadbe6df639a57c86e9b',
 'ah_root_geometry_at_angle.py':'129ad9329f5d9b43f91d486a71926d69ecfc85ebcb62d354627b4d18712fbf57',
 'certify_gibbs_area_refined.py':'b3e0c0200e60927ee5a76d157149586fa6195eb2dee550eeddf7902239caa38b',
 'rational_trig.py':'97e0076b911f51a9c532e14cacb9fb6e3ae7e68707a3c5af77048f9c6623e87f'}
NEGATIVE=(49,50,51,52,53,57,58,59,60,61,62)
POSITIVE={54:('0.139129526671','0.401440263156'),55:('0.139128831434','0.400301757169'),56:('0.139129936470','0.402111340199')}
UNIT_BINS=((0,69),(627,1209),(1809,2391))
UNKNOWN=(69,1209,1808,2391)
SELECT=(1,0,0,0,0,0,0,1,1,1)

def digest(p):return sha256(p.read_bytes()).hexdigest()
def rat(q):
 q=F(q);return arb(q.numerator)/q.denominator
def dot(u,v):return u[0]*v[0]+u[1]*v[1]
def perp(u):return (-u[1],u[0])
def rotate(u,t):
 sn,cs=t.sin_cos();return (cs*u[0]-sn*u[1],sn*u[0]+cs*u[1])
def unit(t):
 sn,cs=t.sin_cos();return (cs,sn)
def bound(x):return str(x.upper())

class Geometry:
 def __init__(self,comp):
  assert all(digest(comp/name)==expected for name,expected in PINS.items())
  sys.path.insert(0,str(comp))
  from certify_ah_released_tangent import load_geometry
  self.ns=load_geometry(comp);ctx.prec=192;self.pi=arb.pi();self.A=self.ns['A'];self.root3=arb(3).sqrt()
  z,D,d=(self.ns[k] for k in ('z','D','direct'))
  self.normal=tuple(-x for x in perp(rotate(d,rat(F(7,400)))))
  pivot=tuple(x+rat(F(43,100))*y for x,y in zip(z,D))
  self.height=dot(self.normal,pivot)
  assert dot(self.normal,self.normal).contains(1)
 def qa(self,u):return ((-u[0]+self.root3*u[1])/2,(self.root3*u[0]+u[1])/2)
 def tangent(self,j):
  pts,ders=self.ns['point_jets'](self.A(F(j,128)))
  q=self.qa(pts[2]);dq=self.qa(ders[2]);length=dot(dq,dq).sqrt();assert length>0
  u=tuple(-x/length for x in perp(dq));v=self.qa(u)
  assert dot(u,u).contains(1) and dot(v,v).contains(1)
  return v,dot(q,u)

def positive_upper(v):return max(arb(0),v.upper())

def dual(g,j,raw):
 pi=g.pi;theta=61*pi/1800;den=(3*theta).sin();assert den>0
 w0,wt,wf=map(rat,raw['w']);nu=list(map(rat,raw['nu']));assert w0>0 and wt>0 and wf>0
 v,h=g.tangent(j);sources=[];constant=wf/(2*den)-w0*g.height-wt*h;pc=arb(0)
 def add(u,w):
  nonlocal constant,pc
  angle=arb.atan2(u[1],u[0]);sgn=1
  if angle<0:angle+=2*pi
  if angle>pi:angle-=pi;sgn=-1
  assert 0<angle<pi
  sources.append((angle,w*sgn));constant+=w*(rat(F(1-sgn,2))+u[0]/2);pc+=w*u[1]
 add(g.normal,w0);add(v,wt)
 add(unit(theta),-wf/den);add(unit(theta+pi/3),wf/den);add(unit(theta+2*pi/3),-wf/den)
 add(unit(pi/3),-nu[0]);add(unit(2*pi/3),-nu[1]);constant+=nu[0]/2+nu[1]/2+nu[2]
 edges=[arb(0)]+sorted((a for a,w in sources),key=lambda a:float(a.mid()))+[pi]
 assert all(b>a for a,b in zip(edges,edges[1:]))
 integral=arb(0);segments=8192
 for left,right in zip(edges,edges[1:]):
  mid=(left+right)/2;aa=-nu[2];bb=-nu[3]
  for angle,w in sources:
   if angle>mid:
    sn,cs=angle.sin_cos();aa-=w*cs;bb+=w*sn
  mass=(abs(aa)+abs(bb)).upper();step=(right-left)/segments
  prev=positive_upper(aa*left.sin()+bb*left.cos())
  for i in range(1,segments+1):
   x=left+step*i;nxt=positive_upper(aa*x.sin()+bb*x.cos());integral+=step*(prev+nxt)/2;prev=nxt
  integral+=mass*(right-left)*step**2/8
 upper=constant+abs(pc)+integral
 assert upper<0
 return {'j':j,'dual_upper':bound(upper),'weights':raw['w'],'equality_multipliers':raw['nu'],'subdivisions':segments}

def counterexample(g,j,initial):
 pi=g.pi;N=2880;p,fixed=map(arb,initial)
 def edge(k):return pi*k/N
 def integ(a,l,r):
  left=edge(l);right=edge(r)
  if a<=left:return arb(0)
  if a<right:right=a
  return (a-right).cos()-(a-left).cos()
 def sm(l,r):return edge(l).cos()-edge(r).cos()
 def cm(l,r):return edge(r).sin()-edge(l).sin()
 def base(a):return a.cos()/2+p*a.sin()+sum(integ(a,l,r) for l,r in UNIT_BINS)+fixed*integ(a,626,627)
 mat=arb_mat([[integ(a,k,k+1) for k in UNKNOWN] for a in (pi/3,2*pi/3)]+[[sm(k,k+1) for k in UNKNOWN],[cm(k,k+1) for k in UNKNOWN]])
 rhs=arb_mat([[arb(1)/2-base(pi/3)],[arb(1)/2-base(2*pi/3)],
  [1-sum(sm(l,r) for l,r in UNIT_BINS)-fixed*sm(626,627)],[-sum(cm(l,r) for l,r in UNIT_BINS)-fixed*cm(626,627)]])
 solution=mat.solve(rhs);rho={k:solution[i,0] for i,k in enumerate(UNKNOWN)};rho[626]=fixed
 assert all(0<v<1 for v in rho.values())
 def rho_at(i):return rho[i] if i in rho else arb(int(any(l<=i<r for l,r in UNIT_BINS)))
 def h(a):return a.cos()/2+p*a.sin()+sum(integ(a,l,r) for l,r in UNIT_BINS)+sum(v*integ(a,k,k+1) for k,v in rho.items())
 def hdir(v):
  angle=arb.atan2(v[1],v[0])
  if angle<0:angle+=2*pi
  if angle>pi:return 1-h(angle-pi)
  return h(angle)
 v,height=g.tangent(j);ah=hdir(g.normal)-g.height;reflected=hdir(v)-height;assert ah>0 and reflected>0
 sigma=rat(F(15497,10000))*pi/180;others=[];i=0
 for k in range(0,360,60):
  for sign in (1,-1):
   if (k,sign) in ((0,1),(240,-1)):continue
   choice=3 if SELECT[i]==0 else 5;i+=1
   u=unit(pi/6+sigma+choice*pi/3);u=rotate((u[0],sign*u[1]),-k*pi/180)
   violation=hdir(u)-arb(1)/2;assert violation>0
   others.append({'rotation_degrees':-k,'reflection':sign,'cut':'m'+str(choice),'violation_lower':str(violation.lower())})
 assert len(others)==10
 sub=16;step=pi/(N*sub);sn,cs=step.sin_cos();hs=[arb(1)/2];gs=[p]
 for i in range(N):
  r=rho_at(i)
  for _ in range(sub):
   h0,g0=hs[-1],gs[-1];hs.append(r+(h0-r)*cs+g0*sn);gs.append(-(h0-r)*sn+g0*cs)
 L=N*sub//3;err=6*step**2/8
 def fv(i):return hs[i]-hs[i+L]+hs[i+2*L]-arb(1)/2
 def fd(i):return gs[i]-gs[i+L]+gs[i+2*L]
 assert fd(0)+6*step<0 and fd(L)-6*step>0
 interior=max(fv(i).upper() for i in range(1,L));assert interior < -err
 return {'j':j,'p':initial[0],'rho626':initial[1],'rho_solution':{str(k):str(rho[k]) for k in UNKNOWN},
  'AH_violation_lower':str(ah.lower()),'QA_tangent_violation_lower':str(reflected.lower()),'other_poses':others,
  'F_interior_grid_upper':str(interior),'F_interpolation_error_upper':str(err.upper()),
  'F_prime_first_upper':str(fd(0).upper()),'F_prime_last_lower':str(fd(L).lower()),'F_subdivisions':L}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('computation',type=Path);ap.add_argument('output',type=Path);args=ap.parse_args()
 assert __debug__
 params_path=Path(__file__).with_name('dual-parameters.json');parameters=json.loads(params_path.read_text());assert sorted(map(int,parameters))==list(NEGATIVE)
 geometry=Geometry(args.computation.resolve())
 proof=[dual(geometry,j,parameters[str(j)]) for j in NEGATIVE]
 witnesses=[counterexample(geometry,j,POSITIVE[j]) for j in POSITIVE]
 result={'schema':'B1-curvature-classification-v1','sigma_degrees':'1.5497','AH_pivot':'43/100','AH_delta_radians':'7/400','tangent_parameter':'j/128','F_required':'F(61*pi/1800)<=0','negative_indices':list(NEGATIVE),'counterexample_indices':list(POSITIVE),'geometry_source_sha256':PINS,'dual_parameters_sha256':digest(params_path),'certifier_sha256':digest(Path(__file__)),'continuous_duals':proof,'counterexamples':witnesses}
 args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'output':str(args.output),'dual_certified':len(proof),'counterexamples':len(witnesses)}))
if __name__=='__main__':main()
