"""Arb certificates for disk, Reuleaux 3 and 5 via finite directional duals."""
from flint import arb,arb_mat,ctx
from fractions import Fraction as F
from pathlib import Path
import gzip,json,sys,math,hashlib,time
ctx.prec=160;p=Path(__file__).resolve().parent;candidate_bytes=gzip.decompress((p/'directional-candidates.json.gz').read_bytes());data=json.loads(candidate_bytes);n=data['n'];N=data['angle_n'];assert n>0 and N>0 and [r['i'] for r in data['results']]==list(range(N+1));A=lambda x:arb(str(x));pi=arb.pi();tau=pi/A(n);sinc=tau.sin()/tau;r=A('0.5');R3=1/(2*(pi/6).cos());R5=1/(2*(pi/10).cos());
vec=lambda t:(t.cos(),t.sin());dot=lambda u,v:u[0]*v[0]+u[1]*v[1];dirs=[vec(pi*A(F(2*j+1,n))) for j in range(n)];v3=[tuple(R3*z for z in vec(pi*A(F(2*j,3)))) for j in range(3)];v5=[tuple(R5*z for z in vec(pi*A(F(2*j,5)))) for j in range(5)];minlow=None;rows=[];start=__import__('time').time()
for rec in data['results']:
 i=rec['i'];l3=[];l5=[];ks=[];kind=[]
 for j,row in enumerate(rec['lambda']):
  a,b=sum(row[:3]),sum(row[3:]);sa=F(str(round(a,12)));sb=F(str(round(b,12)));s=a+b
  ka=max(range(3),key=lambda k:row[k]);kb=max(range(5),key=lambda k:row[k+3]);ks.append((ka,kb))
  if s<1e-8:sa=sb=F(0)
  elif s>1-1e-8:
   if a<1e-8:sa,sb=F(0),F(1)
   elif b<1e-8:sa,sb=F(1),F(0)
   else:sb=1-sa;kind.append((j,'m'))
  else:
   assert (a<1e-8)!=(b<1e-8),(i,j,a,b)
   if a<1e-8:sa=F(0);kind.append((j,'b'))
   else:sb=F(0);kind.append((j,'a'))
  assert 0<=sa<=1 and 0<=sb<=1 and sa+sb<=1,(i,j)
  l3.append(A(sa));l5.append(A(sb))
 assert len(kind)==4,(i,kind)
 res=[sum((l3[j]*dirs[j][k] for j in range(n)),A(0)) for k in range(2)]+[sum((l5[j]*dirs[j][k] for j in range(n)),A(0)) for k in range(2)]
 columns=[]
 for j,k in kind:
  x,y=dirs[j];columns.append([x if k in ('a','m') else A(0),y if k in ('a','m') else A(0),x if k=='b' else -x if k=='m' else A(0),y if k=='b' else -y if k=='m' else A(0)])
 matrix=arb_mat([[columns[j][i] for j in range(4)] for i in range(4)]);rhs=arb_mat([[-v] for v in res]);adjust=matrix.solve(rhs)
 for t,(j,k) in enumerate(kind):
  v=adjust[t,0]
  if k=='a':l3[j]+=v;assert l3[j]>0 and l3[j]<1 and l5[j].is_zero()
  elif k=='b':l5[j]+=v;assert l5[j]>0 and l5[j]<1 and l3[j].is_zero()
  else:l3[j]+=v;l5[j]=1-l3[j];assert l3[j]>0 and l3[j]<1
 # λ·normal = 0 for both translations by the exact interval-linear solve.
 C=A(0);X=A(0);Y=A(0)
 for j,(ka,kb) in enumerate(ks):
  u=dirs[j];v=v5[kb];C+=l3[j]*dot(u,v3[ka])-(l3[j]+l5[j])*r/sinc;X+=l5[j]*dot(u,v);Y+=l5[j]*dot(u,(-v[1],v[0]))
 center=pi*A(F(i,15*N));half=pi*A(F(1,30*N));co,si=center.cos(),center.sin();mid=X*co+Y*si;der=-X*si+Y*co;cosloss=(half.cos()-1).union(A(0));sinbox=(-half.sin()).union(half.sin());bound=C+mid+mid*cosloss+der*sinbox
 value=(pi/4+pi/(2*A(n))*sinc*bound).lower();rows.append(str(value));minlow=value if minlow is None or value<minlow else minlow
 if i%20==0:print(i,'current minimum',minlow,'elapsed',time.time()-start,flush=True)
assert minlow>A('.83256418')
result={'n':n,'angle_n':N,'precision_bits':ctx.prec,'lower':str(minlow),'positive_all':all(arb(v)>A('.8325') for v in rows),'per_cell_lowers':rows,'input_sha256':hashlib.sha256(candidate_bytes).hexdigest(),'scope':'all placements of unit disk, regular Reuleaux triangle and pentagon'};Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print('global lower',minlow,'all > .8325',result['positive_all'],flush=True)
