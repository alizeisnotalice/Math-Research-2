"""Independent Fraction/O_K verifier: never imports the generating SymPy code."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,math

def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):
 m,n=x;u,v=y;return m*u+n*v,m*v+n*u+n*v

def main():
 if not __debug__:raise RuntimeError('no -O')
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert not args.out.exists();data=json.loads(args.input.read_text());assert data['D']==5 and data['basis']=='omega=(1+sqrt(5))/2'
 A,B=[list(map(tuple,data[k])) for k in ('A','B')];assert len(A)==len(B)==6 and A[0]==B[0]==(1,0) and A[5]==B[5]==(0,0);assert add(A[1],B[1])==(1,0)
 lo,hi=[F(*x) for x in data['sqrt5_interval']];assert lo>0 and lo*lo<5<hi*hi
 bounds=[];minimal=[]
 def interval(z):
  m,n=z;u=(2*m+n+n*(lo if n>=0 else hi))/2;v=(2*m+n+n*(hi if n>=0 else lo))/2;return u,v
 for label,poly in [('A',A),('B',B)]:
  for i,z in enumerate(poly):
   assert all(type(x) is int for x in z)
   if z in ((0,0),(1,0)):continue
   low,up=interval(z);assert 0<low<=up<1;bounds.append([label,i,str(low),str(up)])
   m,n=z;T,N=2*m+n,m*m+m*n-n*n;assert add(add(mul(z,z),(-T*z[0],-T*z[1])),(N,0))==(0,0)
   assert n!=0 and math.isqrt(5*n*n)**2!=5*n*n;minimal.append([label,i,[N,-T,1]])
 C=[(0,0)]*11
 for i,x in enumerate(A):
  for j,y in enumerate(B):C[i+j]=add(C[i+j],mul(x,y))
 assert C[:6]==[(1,0)]*6 and C[6]==tuple(data['expected_C6'])
 low,up=interval(C[6]);assert 0<low<=up<1
 out={'status':'passed_exact_integral_local_prefix_only','C':C,'C6_bounds':[str(low),str(up)],'coefficient_bounds':bounds,'minimal_polynomials_ascending':minimal,'first_nonbinary_index':6,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'scope':data['scope']};t=args.out.with_suffix('.partial');t.write_text(json.dumps(out,indent=2)+'\n');os.replace(t,args.out);print(json.dumps(out))
if __name__=='__main__':main()
