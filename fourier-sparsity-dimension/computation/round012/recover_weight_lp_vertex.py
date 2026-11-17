import math,random,json,time,hashlib,sys,argparse
from pathlib import Path
import scipy
from fractions import Fraction
import sympy as sp
from scipy.optimize import linprog
N=64;D=16;W=list(range(8,33,2))
def K(j,i):return sum((-1)**r*math.comb(i,r)*math.comb(N-i,j-r) for r in range(max(0,j-(N-i)),min(i,j)+1))
def vec(j):return [K(j,w)+(0 if w==32 else K(j,64-w)) for w in W],K(j,0)+K(j,64)
def sysm():
 A=[];b=[];E=[[2 if w<32 else 1 for w in W]];e=[(1<<D)-2];c,k=vec(2);E+=[c];e+=[-k]
 for j in range(4,33,2):
  c,k=vec(j);A.append([-z for z in c]);b.append(k-(1<<D)*abs(K(j,28)))
  if j>=8:
   sg=1 if j%4==0 else -1;r=[-z for z in c];r[W.index(j)]+=1<<(D+1);A.append(r);b.append(k+sg*(1<<D)*K(j,28))
 assert len(A)==28 and len(E)==2
 return A,b,E,e
def verify(x):
 # Independent direct enumerator form, intentionally not reusing sysm inequalities.
 full={0:Fraction(1),64:Fraction(1)}
 for w,z in zip(W,x):full[w]=z;full[64-w]=z
 if any(z<0 for z in x) or sum(full.values())!=1<<D:return False
 def B(j):return sum(full[i]*K(j,i) for i in full)/Fraction(1<<D)
 if B(2)!=0:return False
 for j in range(0,33,2):
  if B(j)<abs(K(j,28)):return False
  if j>=8:
   signed=(B(j)+(K(j,28) if j%4==0 else -K(j,28)))/2
   if signed<full[j]:return False
 return True
def main():
 parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True);args=parser.parse_args()
 t=time.time();A,b,E,e=sysm();rng=random.Random(2026111712);ans=None
 for trial in range(100):
  res=linprog([rng.randint(-100,100) for _ in W],A_ub=A,b_ub=b,A_eq=E,b_eq=e,bounds=[(0,None)]*len(W),method='highs')
  if not res.success:continue
  active=[i for i,(r,q) in enumerate(zip(A,b)) if abs(sum(v*z for v,z in zip(r,res.x))-q)<1e-5]
  active += [28+i for i,z in enumerate(res.x) if abs(z)<1e-7]
  rows=E[:];rhs=e[:]
  for i in active:
   r,q=(A[i],b[i]) if i<28 else ([1 if k==i-28 else 0 for k in range(13)],0)
   if sp.Matrix(rows+[r]).rank()>len(rows):rows.append(r);rhs.append(q)
   if len(rows)==13:break
  if len(rows)==13:
   x=list(sp.Matrix(rows).inv()*sp.Matrix(rhs));x=[Fraction(int(z.p),int(z.q)) for z in x]
   if verify(x):ans=x;break
 if ans is None:raise RuntimeError('No rational feasible vertex recovered')
 o={'status':'passed' if ans else 'no_exact_vertex_recovered','weights':W,'rational_A_weights':None if ans is None else [str(x) for x in ans],'fraction_verified':bool(ans),'trials':trial+1,'source_sha256':hashlib.sha256(open(__file__,'rb').read()).hexdigest(),'elapsed_seconds':time.time()-t,'python':sys.version,'scipy':scipy.__version__,'sympy':sp.__version__,'seed':2026111712,'scope':'Continuous rational weight-enumerator relaxation only; not an integer code.'};open(args.out,'w').write(json.dumps(o,indent=2)+'\n');print(json.dumps({'status':o['status'],'trials':o['trials'],'elapsed':o['elapsed_seconds']}))
if __name__=='__main__':main()
