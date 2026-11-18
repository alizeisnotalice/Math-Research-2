'''Replay 11 dual intervals and one continuous all-pose obstruction interval.'''
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import sys,json,argparse
from flint import arb,arb_series,arb_mat,ctx
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from certify_b1_classification import Geometry,dual,rat,dot,perp,POSITIVE,UNKNOWN,UNIT_BINS,counterexample,PINS
ap=argparse.ArgumentParser();ap.add_argument('computation',type=Path);ap.add_argument('output',type=Path);args=ap.parse_args()
G=Geometry(args.computation.resolve());ctx.prec=192;ctx.cap=3;pi=arb.pi();root3=arb(3).sqrt();R=1/root3

def norm(v):return dot(v,v).sqrt()
def norm_up(v):return sum((abs(x).upper()**2 for x in v),arb(0)).sqrt().upper()
def jets(w):
 y=arb_series([w,1,0]);a,b=G.ns['points'](y)[0][2];q=G.qa((a,b));p=(q[0][0],q[1][0]);v=(q[0][1],q[1][1]);acc=(2*q[0][2],2*q[1][2]);speed=norm(v);assert speed>0
 u=tuple(-x/speed for x in perp(v));h=dot(v,acc);du=tuple(-x for x in perp(tuple(acc[i]/speed-v[i]*h/speed**3 for i in range(2))))
 return p,v,u,du

params=json.loads((P/'dual-parameters.json').read_text());radii={49:F(4,100000),50:F(25,1000000),51:F(15,1000000),52:F(6,1000000),53:F(4,10000000),57:F(15,10000000),58:F(8,1000000),59:F(18,1000000),60:F(3,100000),61:F(45,1000000),62:F(65,1000000)}
dual_rows=[]
for j,r in radii.items():
 center=F(j,128);lo,hi=center-r,center+r
 rate=arb(0)
 for k in range(128):
  a=lo+(hi-lo)*k/128;b=lo+(hi-lo)*(k+1)/128
  p,v,u,du=jets(rat(a).union(rat(b)));value=(R+norm_up(p))*norm_up(du);rate=max(rate,value.upper())
 record=dual(G,j,params[str(j)]);margin=-arb(record['dual_upper']);weight=rat(F(params[str(j)]['w'][1]));remain=margin-weight*rate*rat(r)
 assert remain>0
 dual_rows.append({'j':j,'lo':str(lo),'hi':str(hi),'rate_upper':str(rate.upper()),'remaining_margin_lower':str(remain.lower())})

# Reconstruct and reverify the explicit j56 body, including F and 12 original pose conditions.
witness=counterexample(G,56,POSITIVE[56]);N=2880;p,fixed=map(arb,POSITIVE[56]);ones=UNIT_BINS;unknown=UNKNOWN
def edge(k):return pi*k/N
def integ(a,l,r):
 left=edge(l);right=edge(r)
 if a<=left:return arb(0)
 if a<right:right=a
 return (a-right).cos()-(a-left).cos()
def d_integ(a,l,r):
 left=edge(l);right=edge(r)
 if a<=left:return arb(0)
 if a<right:right=a
 return (a-left).sin()-(a-right).sin()
def ms(l,r):return edge(l).cos()-edge(r).cos()
def mc(l,r):return edge(r).sin()-edge(l).sin()
def base(a):return a.cos()/2+p*a.sin()+sum(integ(a,l,r) for l,r in ones)+fixed*integ(a,626,627)
mat=arb_mat([[integ(a,k,k+1) for k in unknown] for a in (pi/3,2*pi/3)]+[[ms(k,k+1) for k in unknown],[mc(k,k+1) for k in unknown]])
b=arb_mat([[arb(1)/2-base(pi/3)],[arb(1)/2-base(2*pi/3)],[1-sum(ms(l,r) for l,r in ones)-fixed*ms(626,627)],[-sum(mc(l,r) for l,r in ones)-fixed*mc(626,627)]])
s=mat.solve(b);rho={k:s[i,0] for i,k in enumerate(unknown)};rho[626]=fixed;assert all(0<x<1 for x in rho.values())
def hp(a):return -a.sin()/2+p*a.cos()+sum(d_integ(a,l,r) for l,r in ones)+sum(z*d_integ(a,k,k+1) for k,z in rho.items())
def h(a):return base(a)-fixed*integ(a,626,627)+sum(z*integ(a,k,k+1) for k,z in rho.items())
def support_and_point(v):
 angle=arb.atan2(v[1],v[0]);
 if angle<0:angle+=2*pi
 sign=1
 if angle>pi:angle-=pi;sign=-1
 assert 0<angle<pi
 for l,r in (*ones,*[(k,k+1) for k in rho]):
  for idx in (l,r):assert angle<edge(idx) or angle>edge(idx)
 hv=h(angle);dv=hp(angle)
 if sign<0:hv=1-hv;dv=-dv
 Jv=perp(v);return hv,(hv*v[0]+dv*Jv[0],hv*v[1]+dv*Jv[1])

def test_cell(a,b):
 mid=(a+b)/2;radius=(b-a)/2;qm,_,um,_=jets(rat(mid));qb,_,_,dub=jets(rat(a).union(rat(b)));v=G.qa(um);hv,point=support_and_point(v);dist=hv-dot(qm,um);xp=G.qa(point);dx=(xp[0]-qb[0],xp[1]-qb[1]);U=norm_up(dub);loss=U*(norm_up(dx)+U*rat(radius))*rat(radius);return dist-loss

lo,hi=F(4192,10000),F(4383,10000);stack=[(lo+(hi-lo)*k/512,lo+(hi-lo)*(k+1)/512,0) for k in range(512)];count=0;deep=0;worst=None
while stack:
 a,b,d=stack.pop();x=test_cell(a,b);count+=1;deep=max(deep,d)
 if x>0:
  y=x.lower();worst=y if worst is None or y<worst else worst
  continue
 assert d<20,(a,b,x)
 mid=(a+b)/2;stack.append((mid,b,d+1));stack.append((a,mid,d+1))
result={'dual_intervals':dual_rows,'obstruction':{'lo':str(lo),'hi':str(hi),'body_j':56,'covered_binary_cells':count,'max_depth':deep,'minimum_margin_lower':str(worst),'witness':witness},'geometry_source_sha256':PINS,'dual_parameters_sha256':sha256((P/'dual-parameters.json').read_bytes()).hexdigest(),'base_certifier_sha256':sha256((P/'certify_b1_classification.py').read_bytes()).hexdigest(),'extension_certifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'note':'each interval is a closed certificate; gaps are not decided'}
args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('dual',len(dual_rows),'obstruction cells',count,'depth',deep,'margin',worst)
