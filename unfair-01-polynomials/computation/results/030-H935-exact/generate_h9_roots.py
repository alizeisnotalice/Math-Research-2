"""Generate rational-disk input for H9(x^3)H9(x^5), independently of search."""
import hashlib,json,mpmath as mp
from fractions import Fraction
from pathlib import Path
from exact_complex import mother,lift

HERE=Path(__file__).resolve().parent; S=10**160; MR=Fraction(1,10**70); LR=Fraction(1,10**65)
G=[1,0,0,-1,0,1]
PHI={5:[1,1,1,1,1],15:[1,-1,0,1,-1,1,0,-1,1],25:[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}
def frac(x):return Fraction(mp.nstr(x,110))
def disk(z,r):return {'center':{'re':f'{frac(mp.re(z)).numerator}/{frac(mp.re(z)).denominator}','im':f'{frac(mp.im(z)).numerator}/{frac(mp.im(z)).denominator}'},'radius':f'{r.numerator}/{r.denominator}'}
def rat(x):return f'{x.numerator}/{x.denominator}'
def main():
 mp.mp.dps=120
 roots=mp.polyroots(list(reversed(G)),maxsteps=1000,error=False)
 neg=next(z for z in roots if abs(mp.im(z))<mp.mpf('1e-90'))
 ups=sorted([z for z in roots if mp.im(z)>0],key=lambda z:mp.re(z)); assert len(ups)==2
 parent_reps=[neg]+ups; mothers=[]
 for ident,z in enumerate([neg,ups[0],ups[1],mp.conj(ups[0]),mp.conj(ups[1])]):
  c=(frac(mp.re(z)),frac(mp.im(z)));L,R=mother(G,c,MR,S);assert L<R;mothers.append({'index':ident,'disk':disk(z,MR),'lhs':rat(L),'rhs':rat(R)})
 groups=[]
 # 16 fixed cyclotomic upper-half quadratic groups.
 for n in (5,15,25):
  for j in range(1,(n+1)//2):
   if __import__('math').gcd(j,n)!=1:continue
   z=mp.e**(2j*mp.pi*j/n);c=(frac(mp.re(z)),frac(mp.im(z)));L,R=mother(PHI[n],c,LR,S);assert L<R
   groups.append({'index':len(groups),'source':'cyclotomic','n':n,'j':j,'kind':'quadratic','disk':disk(z,LR),'lhs':rat(L),'rhs':rat(R)})
 assert len(groups)==16
 # G5 lift groups, with one representative for each conjugate pair.
 for parent,z in enumerate(parent_reps):
  c=(frac(mp.re(z)),frac(mp.im(z)))
  for q in (3,5):
   lifts=mp.polyroots([1]+[0]*(q-1)+[-z],maxsteps=1000,error=False)
   selected=[]
   if parent==0:
    selected += [('linear',next(w for w in lifts if abs(mp.im(w))<mp.mpf('1e-90')))]
    selected += [('quadratic',w) for w in lifts if mp.im(w)>0]
   else:selected += [('quadratic',w) for w in lifts]
   assert len(selected)==(1+(q-1)//2 if parent==0 else q)
   for kind,w in selected:
    cc=(frac(mp.re(w)),frac(mp.im(w)));L,R=lift(c,MR,cc,LR,S,q);assert L<R
    groups.append({'index':len(groups),'source':'G5_lift','parent_pair':parent,'exponent':q,'kind':kind,'disk':disk(w,LR),'lhs':rat(L),'rhs':rat(R)})
 assert len(groups)==37
 out={'format':'H9-35-rational-rouche-v1','G5_coefficients_ascending':G,'scale':str(S),'mother_radius':'1/10^70','lift_radius':'1/10^65','mothers':mothers,'groups':groups,'coverage':{'mother_disks':5,'cyclotomic_quadratic_groups':16,'G5_lift_quadratic_groups':19,'G5_lift_linear_groups':2,'target_root_disks':72}}
 p=HERE/'H9-35-root-balls.json';p.write_text(json.dumps(out,indent=2)+'\n');print(hashlib.sha256(p.read_bytes()).hexdigest())
if __name__=='__main__':main()
