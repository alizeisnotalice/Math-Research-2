import json,hashlib,mpmath as mp
from fractions import Fraction
from pathlib import Path
from root_ball_common import complex_to_json,complex_from_json,parse_fraction,polynomial_rouche_bounds,lift_rouche_bounds,disk_disjoint,disk_nonreal
H=[1,1,1,0,0,0,0,1,1];S=10**160;MR=Fraction(1,10**70);LR=Fraction(1,10**65)
def f(x):return Fraction(mp.nstr(x,105))
def d(z,r):return {'center':complex_to_json((f(mp.re(z)),f(mp.im(z)))),'radius':f'{r.numerator}/{r.denominator}'}
def main():
 mp.mp.dps=110;r=mp.polyroots(list(reversed(H)),maxsteps=1000,error=False);assert not [z for z in r if abs(mp.im(z))<mp.mpf('1e-90')]
 parents=sorted([z for z in r if mp.im(z)>0],key=lambda z:mp.re(z));assert len(parents)==4;moms=[];gs=[]
 for b,z in enumerate(parents):
  c=(f(mp.re(z)),f(mp.im(z)));L,R=polynomial_rouche_bounds(H,c,MR,S);assert L<R and disk_nonreal(c,MR);moms.append({'block':b,'disk':d(z,MR),'lhs':f'{L.numerator}/{L.denominator}','rhs':f'{R.numerator}/{R.denominator}'})
  for q,fam in ((3,'a'),(5,'b')):
   xs=sorted(mp.polyroots([1]+[0]*(q-1)+[-z],maxsteps=1000,error=False),key=lambda w:mp.arg(w))
   for slot,w in enumerate(xs):
    cc=(f(mp.re(w)),f(mp.im(w)));L,R=lift_rouche_bounds(c,MR,cc,LR,S,q);assert L<R and disk_nonreal(cc,LR);gs.append({'index':len(gs),'block':b,'family':fam,'exponent':q,'slot':slot,'kind':'quadratic','disk':d(w,LR),'lhs':f'{L.numerator}/{L.denominator}','rhs':f'{R.numerator}/{R.denominator}'})
 assert len(gs)==32;ds=[]
 for g in gs:
  c=complex_from_json(g['disk']['center']);rr=parse_fraction(g['disk']['radius']);ds += [(c,rr),((c[0],-c[1]),rr)]
 assert len(ds)==64
 for i,(c,rr) in enumerate(ds):
  for z,t in ds[:i]:assert disk_disjoint(c,rr,z,t)
 o={'format':'H8-a3b5-rational-rouche-v1','H':H,'a':3,'b':5,'mp_dps':110,'scale':str(S),'mother_radius':'1/10^70','lift_radius':'1/10^65','mothers':moms,'groups':gs,'coverage':{'mother_pairs':4,'quadratic_groups':32,'root_disks':64,'pairwise_disjoint':True},'order':'mothers increasing real part; each a/b lift increasing argument'};Path('H8-a3b5-root-balls.json').write_text(json.dumps(o,indent=2)+'\n');print(hashlib.sha256(Path('H8-a3b5-root-balls.json').read_bytes()).hexdigest())
if __name__=='__main__':main()
