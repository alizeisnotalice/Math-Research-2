from pathlib import Path
import sys,json,hashlib,math,time,subprocess
sys.path.insert(0,str(Path(__file__).resolve().parent))
# Run the complete angle certificate before consuming its stated lower bounds.
subprocess.run([sys.executable, str(Path(__file__).resolve().parent/'verify_cert.py'), sys.argv[1]], check=True)
from fractions import Fraction as F
# Import the same pinned supports and hosts as the angular verifier.
from cross_common import q, arb, pi, rot, axis, hf, hp, Host, C, A, dot, ctx
D=Path(__file__).resolve().parent
pose155=json.loads((C/'free-p30-reoptimized-pose.json').read_text());pose175=json.loads((C/'seven-point-families-pose.json').read_text())['pose'];hosts={'r155':Host(pose155,False),'r175':Host(pose175,True)}
# Rigorous host outer radius from 64 rational directions. Width of adjacent angular gaps is interval checked.
K=64; seed=[]
for i in range(K):
 a=2*math.pi*i/K;s=1 if math.cos(a)>=0 else -1;b=(a-(0 if s==1 else math.pi)+math.pi)%(2*math.pi)-math.pi;seed.append((s,round(1e12*math.tan(b/2))))
uv=[axis(v)[0] for v in seed]; assert all(u[0]*v[1]-u[1]*v[0]>0 for u,v in zip(uv,uv[1:]+uv[:1]));limit=(2*pi/64+A('.0001')).cos();assert all(q(u[0]*v[0]+u[1]*v[1])>limit for u,v in zip(uv,uv[1:]+uv[:1]));hcap=A('.59')*(pi/64+A('.00005')).cos()
for name,h in hosts.items():
 for z in seed:assert h.row(z)[2]<hcap,(name,z,h.row(z)[2],hcap)
# Disk radius .5 belongs to both hosts. Triangle-cap bound is rational at each asserted escape delta.
deltas={'r155-r151':F('.00154117'),'r155-r216':F('.00149681'),'r175-r151':F('.00145235'),'r175-r216':F('.00141054')}
# For the disk of radius 1/2, the two tangents yield an exterior triangle
# of area d^2/(2 sqrt((b+d)^2-1/4)); 8/25 bounds that square root.
assert all((F(59,100)+d)**2-F(1,4)<F(8,25)**2 for d in deltas.values())
low={key:d*d/F(16,25) for key,d in deltas.items()}
# Exact halfspace intersection of rational halfspaces and square [-2,2]^2.
def clip(poly,u,b):
 out=[]
 for p,r in zip(poly,poly[1:]+poly[:1]):
  fp=u[0]*p[0]+u[1]*p[1]-b;fr=u[0]*r[0]+u[1]*r[1]-b
  if fp<=0:out.append(p)
  if (fp<0 and fr>0) or (fp>0 and fr<0):
   t=fp/(fp-fr);out.append((p[0]+t*(r[0]-p[0]),p[1]+t*(r[1]-p[1])))
 return out
def ceiling(x,scale=10**12):return F(int((x*scale).upper().ceil().unique_fmpz()),scale)
# 1024 directions; seed rational gives exactly unit normal.
N=1024;seeds=[]
for i in range(N):
 a=2*math.pi*i/N;s=1 if math.cos(a)>=0 else -1;b=(a-(0 if s==1 else math.pi)+math.pi)%(2*math.pi)-math.pi;seeds.append((s,round(10**12*math.tan(b/2))))
areas=[];search=json.loads((D/'area-search-final.json').read_text())['search'];assert len(search)==4 and {(r['host'],r['body']) for r in search}=={(h,b) for h in ('r155','r175') for b in ('r151','r216')};assert all(r['sign'] in (-1,1) and len(r['pose'])==3 and all(math.isfinite(x) for x in r['pose']) for r in search);orig_lo={'r155':F('.83751202649'),'r175':F('.8378949872364058')}
for row in search:
 name=row['host'];body=row['body'];sign=row['sign']; a,x,y=[q(F(format(t,'.12f'))) for t in row['pose']];h=hosts[name];f=hp if body=='r151' else hf
 poly=[(F(-2),F(-2)),(F(2),F(-2)),(F(2),F(2)),(F(-2),F(2))]
 for z in seeds:
  u,v,H=h.row(z); w=rot(-a,v);w=(w[0],sign*w[1]); height=(f(w)+x*v[0]+y*v[1]).upper();b=ceiling(max(H,height));poly=clip(poly,u,b);assert len(poly)>=3
 # Cover every rational vertex by its 10^-12 square; integer hull contains the full polygon.
 scale=10**12;points=set()
 for px,py in poly:
  ax=(px*scale).__floor__();ay=(py*scale).__floor__()
  for dx in (0,1):
   for dy in (0,1):points.add((ax+dx,ay+dy))
 def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
 pts=sorted(points);lo=[];hi=[]
 for z in pts:
  while len(lo)>=2 and cross(lo[-2],lo[-1],z)<=0:lo.pop()
  lo.append(z)
 for z in reversed(pts):
  while len(hi)>=2 and cross(hi[-2],hi[-1],z)<=0:hi.pop()
  hi.append(z)
 hull=lo[:-1]+hi[:-1];twice=sum(x0*y1-y0*x1 for (x0,y0),(x1,y1) in zip(hull,hull[1:]+hull[:1]));area=F(twice,2*scale*scale);assert area>0
 key=f'{name}-{body}';entry=dict(host=name,input=body,reflection=sign,pose=[format(t,'.12f') for t in row['pose']],upper_area=str(area),upper_decimal_ceil=str(ceiling(q(area),10**12)),added_upper_decimal_ceil=str(ceiling(q(area-orig_lo[name]),10**12)),added_lower_fraction=str(low[key]),added_lower_decimal_floor=str(F(int(low[key]*10**12),10**12)),vertices=len(poly));areas.append(entry);print(key,'upper area',entry['upper_decimal_ceil'],'added interval',entry['added_lower_decimal_floor'],entry['added_upper_decimal_ceil'],'vertices',len(poly),flush=True)
obj=dict(radius_host_upper='.59',host_radius_method='64 rational directions and angular gap',input_radius='.58',area_brackets=areas,host_area_lower={'r155':str(orig_lo['r155']),'r175':str(orig_lo['r175'])},directions=N,precision_bits=ctx.prec,sha256_certificate=hashlib.sha256((D/'certificate-labels.json').read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256((D/'verify_cert.py').read_bytes()).hexdigest(),area_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(D/'area-result.json').write_text(json.dumps(obj,indent=2)+'\n')
