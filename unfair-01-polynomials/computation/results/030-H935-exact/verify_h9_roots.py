"""Exact Fraction/Rouché validation for H9's 72-root cover (labels audited separately)."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
from exact_complex import cj,pf,mother,lift,dis
G=[1,0,0,-1,0,1];PS={5:[1,1,1,1,1],15:[1,-1,0,1,-1,1,0,-1,1],25:[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}
def req(x,s):
 if not x:raise ValueError(s)
def D(x):
 z,r=cj(x['disk']['center']),pf(x['disk']['radius']);req(r>0,'radius');return z,r
def sep(a):
 for i,(z,r) in enumerate(a):
  for w,s in a[:i]:req(dis(z,r,w,s),'overlap')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(b):
 S=int(b['scale']);req(b['format']=='H9-35-rational-rouche-v1' and b['G5_coefficients_ascending']==G and len(b['mothers'])==5 and len(b['groups'])==37,'schema')
 moms={};md=[]
 for x in b['mothers']:
  i=x['index'];z,r=D(x);L,R=mother(G,z,r,S);req(type(i)is int and i not in moms and L<R and pf(x['lhs'])==L and pf(x['rhs'])==R,'G5 mother');moms[i]=(z,r);md.append((z,r))
 req(set(moms)==set(range(5)),'mother index set');req(moms[0][0][1]==0 and moms[0][0][0]+moms[0][1]<0,'negative real mother')
 for up,down in ((1,3),(2,4)):
  z,r=moms[up];w,s=moms[down];req(z[1]>r and w[0]==z[0] and w[1]==-z[1] and s==r,'exact conjugate mother pair')
 req(moms[1][0][0]+moms[1][1] < moms[2][0][0]-moms[2][1], 'upper mother Re interval ordering');sep(md)
 roots=[];cyc=Counter();liftc=Counter();seen=set();labels=set();quad=line=0
 for x in b['groups']:
  i=x['index'];req(type(i)is int and i not in seen and 0<=i<37,'group index');seen.add(i);z,r=D(x);source=x['source'];req(source in ('cyclotomic','G5_lift'),'group source')
  if source=='cyclotomic':
   n,j=x['n'],x['j'];req(n in PS and type(j)is int and 1<=j<n/2 and __import__('math').gcd(n,j)==1 and (n,j) not in labels and x['kind']=='quadratic','cyclo label');labels.add((n,j));L,R=mother(PS[n],z,r,S);req(L<R and pf(x['lhs'])==L and pf(x['rhs'])==R,'cyclo Rouche');cyc[n]+=1;quad+=1;roots.extend(((z,r),((z[0],-z[1]),r)))
  else:
   par,q=x['parent_pair'],x['exponent'];req(type(par)is int and par in (0,1,2) and q in (3,5),'lift label');L,R=lift(*moms[par],z,r,S,q);req(L<R and pf(x['lhs'])==L and pf(x['rhs'])==R,'lift Rouche');liftc[par,q]+=1
   if x['kind']=='linear':req(par==0 and z[1]==0 and z[0]+r<0,'real lift');line+=1;roots.append((z,r))
   else:req(x['kind']=='quadratic' and abs(z[1])>r,'quad lift');quad+=1;roots.extend(((z,r),((z[0],-z[1]),r)))
 req(seen==set(range(37)) and cyc=={5:2,15:4,25:10} and labels=={(5,1),(5,2),(15,1),(15,2),(15,4),(15,7),(25,1),(25,2),(25,3),(25,4),(25,6),(25,7),(25,8),(25,9),(25,11),(25,12)},'cyclo complete labels')
 req(liftc=={(0,3):2,(0,5):3,(1,3):3,(1,5):5,(2,3):3,(2,5):5} and quad==35 and line==2 and len(roots)==72,'lift fibers/counts');sep(roots)
 return {'status':'H9_72_roots_exact_fraction_rouche_passed','G5_mother_disks':5,'Phi5_mother_quadratic_groups':2,'cyclotomic_quadratic_groups':16,'G5_lift_quadratic_groups':19,'G5_lift_linear_groups':2,'target_root_disks':72,'cyclotomic_label_scope':'Labels exhaust primitive (n,j) upper-half slots and are not used as an exact angle proof.'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args();o=verify(json.loads(a.balls.read_text()));o.update(balls_sha256=sha(a.balls),verifier_sha256=sha(Path(__file__)),helper_sha256=sha(Path(__file__).with_name('exact_complex.py')));a.record.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o))
if __name__=='__main__':main()
