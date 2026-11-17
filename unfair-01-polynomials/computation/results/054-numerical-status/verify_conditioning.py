"""Actual floating conditioning regression, with an independently checked exact fair control."""
import argparse,json,hashlib,os,sys
from pathlib import Path
import numpy as np

def conv(a,b):
 out=[0.]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out

def group_poly(roots):
 p=[1+0j]
 for r,i in roots:
  z=complex(r,i);q=[0j]*(len(p)+1)
  for j,v in enumerate(p):q[j]-=z*v;q[j+1]+=v
  p=q
 p=[v/p[0] for v in p];assert max(abs(v.imag) for v in p)<1e-7;return [v.real for v in p]
def balanced(polys):
 if not polys:return [1.]
 ordered=sorted([np.asarray(p) for p in polys],key=lambda p:p[1] if len(p)>1 else 0.);a=[]
 for i in range((len(ordered)+1)//2):
  a.append(ordered[i])
  if len(ordered)-1-i!=i:a.append(ordered[-1-i])
 while len(a)>1:a=[np.convolve(a[i],a[i+1]) if i+1<len(a) else a[i] for i in range(0,len(a),2)]
 return a[0].tolist()
def score(a,b):return max(0.,-min(a),-min(b),max(a)-1,max(b)-1)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert not args.out.exists();x=json.loads(args.input.read_text());task=x['task0'];A=[int(i in x['exact_A_support']) for i in range(3)];B=[int(i in x['exact_B_support']) for i in range(61)];assert conv(A,B)==task['C']
 bits=[0]+[int(c) for c in format(x['mask'],f'0{task["groups"]-1}b')];polys=[[group_poly(g) for bit,g in zip(bits,task['roots']) if bit==side] for side in (1,0)];old=[];new=[]
 for factors in polys:
  p=[1.]
  for f in factors:p=conv(p,f)
  old.append(p);new.append(balanced(factors))
 oldscore,newscore=score(*old),score(*new);olderr=max(abs(v-c) for v,c in zip(conv(*old),task['C']));newerr=max(abs(v-c) for v,c in zip(conv(*new),task['C']))
 assert oldscore>1e-3 and newscore<1e-8 and newerr<1e-8
 out={'status':'passed_floating_conditioning_regression','exact_control_product_verified':True,'old_sequential_score':oldscore,'old_product_residual':olderr,'new_balanced_score':newscore,'new_product_residual':newerr,'numpy':np.__version__,'python':sys.version,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Shows floating reconstruction instability on this control only. Near numerical samples remain unknown until exact certification.'};t=args.out.with_suffix('.partial');t.write_text(json.dumps(out,indent=2)+'\n');os.replace(t,args.out);print(json.dumps(out))
if __name__=='__main__':main()
