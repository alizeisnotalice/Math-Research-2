"""Exact finite checks of the tensor pair and the shear realization, r=0..3."""
from pathlib import Path
import json,hashlib,time,platform,sys
P=Path(__file__).resolve().parent

def fwht(a):
 a=a[:];step=1
 while step<len(a):
  for base in range(0,len(a),2*step):
   for i in range(base,base+step):a[i],a[i+step]=a[i]+a[i+step],a[i]-a[i+step]
  step*=2
 return a

def main():
 start=time.perf_counter();raw=(P/'mm_completion_result.json').read_bytes();sol=json.loads(raw)['solution'];b=int(sol['f_hex'],16);bp=int(sol['f_plus_h_hex'],16);h=b^bp
 records=[]
 for r in range(4):
  n=8+2*r;q=1<<(n//2);mask=(1<<r)-1;F=[];H=[]
  for z in range(1<<n):
   x=z&255;u=(z>>8)&mask;v=z>>(8+r)
   F.append((b>>x&1)^((u&v).bit_count()%2));H.append((h>>x&1)*int(v==0))
  G=[f^hval for f,hval in zip(F,H)];WF=fwht([1-2*v for v in F]);WG=fwht([1-2*v for v in G]);assert all(abs(v)==q for v in WF+WG)
  assert sum(H)==2*q
  # Coordinates z=(y,t,s), with y occupying n low bits.
  shear=[];shearT=[]
  for z in range(1<<(n+2)):
   y=z&((1<<n)-1);t=(z>>n)&1;s=(z>>(n+1))&1
   shear.append(F[y]^(t*(s^1^H[y])))
   shearT.append(F[y]^((t^s)*(s^1^H[y])))
  A=[x^y for x,y in zip(shear,shearT)]
  assert sum(A)==4*q
  assert all(abs(v)==2*q for v in fwht([1-2*v for v in shear])+fwht([1-2*v for v in shearT]))
  assert all(A[z]==(((z>>(n+1))&1)*H[z&((1<<n)-1)]) for z in range(len(A)))
  records.append(dict(r=r,tensor_dimension=n,tensor_q=q,tensor_switch_weight=sum(H),tensor_both_bent=True,shear_dimension=n+2,shear_switch_weight=sum(A),shear_both_bent=True,shear_difference_formula_verified=True))
 out=dict(records=records,elapsed_seconds=time.perf_counter()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(raw).hexdigest(),python=sys.version,platform=platform.platform(),scope='finite r=0..3 checks; all-r claims require the analytic proof')
 (P/'tensor_shear_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
