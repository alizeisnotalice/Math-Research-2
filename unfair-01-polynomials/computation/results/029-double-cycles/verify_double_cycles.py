"""Independent integer verification of H8 double-fiber Frobenius inputs.

Every polynomial product, irreducibility claim and squarefree test is
recomputed over a prime field. Source success flags are ignored.
"""
from pathlib import Path
from math import lcm
import argparse,hashlib,json,platform,time
from finite_fields import _prime,_irreducible_mod,_mul_mod,_gcd_mod,req
H=[1,1,1,0,0,0,0,1,1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def subst(h,k):
 a=[0]*(k*(len(h)-1)+1)
 for i,c in enumerate(h):a[k*i]=c
 return a
def factorization(rows,target,prime):
 req(type(rows)is list and len(rows)>0,'factor list');out=[1];degrees=[]
 for r in rows:
  req(type(r['multiplicity'])is int and r['multiplicity']==1,'squarefree multiplicity')
  a=r['ascending'];req(type(a)is list and all(type(v)is int for v in a),'integer coefficients');a=[v%prime for v in a]
  req(type(r['degree'])is int and r['degree']==len(a)-1 and a[-1]==1,'monic degree')
  _irreducible_mod(a,prime);out=_mul_mod(out,a,prime);degrees.append(len(a)-1)
 req(out==[v%prime for v in target],'full factor product')
 return out,degrees
def verify(data):
 req(data['H_ascending']==H,'fixed H8');out=[]
 for p,q in ((3,5),(5,3)):
  r=data['p'+str(p)];ell=r['ell'];req(_prime(ell) and ell%p==1 and r['p']==p,'prime and fiber order')
  h,dh=factorization(r['H'],H,ell);ap,dp=factorization(r['H_xp'],subst(H,p),ell);aq,dq=factorization(r['H_xq'],subst(H,q),ell)
  c=_mul_mod(ap,aq,ell);derivative=[i*c[i]%ell for i in range(1,len(c))]
  req(_gcd_mod(c,derivative,ell)==[1],'whole C squarefree')
  req(p not in dh,'mother has no p-cycle')
  req(dp.count(p)==2 and all(d%p for d in dp if d!=p) and all(d%p for d in dq),'only two p-cycles across C')
  M=lcm(*([d for d in dp if d!=p]+dq));req(M%p!=0 and type(r['M'])is int and r['M']==M,'correct whole-C killing exponent')
  out.append({'p':p,'q':q,'ell':ell,'M':M,'H_degrees':dh,'H_xp_degrees':dp,'H_xq_degrees':dq,'other_C_cycles_fixed_by_power_M':True,'two_remaining_cycles_are_in_distinct_fixed_mother_fibers':True})
 return out
def main():
 a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--record',type=Path,required=True);z=a.parse_args();start=time.monotonic();v=verify(json.loads(z.input.read_text()));rec={'status':'exact_integer_double_cycle_inputs_verified','certificates':v,'input_sha256':sha(z.input),'verifier_sha256':sha(__file__),'kernel_sha256':sha(Path(__file__).with_name('finite_fields.py')),'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'scope':'Frobenius cycle premises only; root-sum and full-fiber reduction require the separate proof.'};z.record.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
if __name__=='__main__':main()
