import json,hashlib,mpmath as mp
from fractions import Fraction
from pathlib import Path
from root_ball_common import complex_to_json,polynomial_rouche_bounds,lift_rouche_bounds,disk_disjoint,complex_from_json,parse_fraction
H=[1,1,1,0,0,0,1,1];S=10**160;MR=Fraction(1,10**70);LR=Fraction(1,10**65)
def fr(x):return Fraction(mp.nstr(x,105))
def dj(z,r):return {'center':complex_to_json((fr(mp.re(z)),fr(mp.im(z)))),'radius':f'{r.numerator}/{r.denominator}'}
def main():
 mp.mp.dps=110; rr=mp.polyroots(list(reversed(H)),maxsteps=1000,error=False); real=[z for z in rr if abs(mp.im(z))<mp.mpf('1e-90')][0]; pos=sorted([z for z in rr if mp.im(z)>0],key=lambda z:mp.re(z));moms=[]; groups=[]
 for b,z in enumerate(pos+[real]):
  c=(fr(mp.re(z)),fr(mp.im(z)));L,R=polynomial_rouche_bounds(H,c,MR,S);assert L<R;moms.append({'block':b,'disk':dj(z,MR),'lhs':f'{L.numerator}/{L.denominator}','rhs':f'{R.numerator}/{R.denominator}'})
  for q,fam in ((3,'a'),(5,'b')):
   xs=mp.polyroots([1]+[0]*(q-1)+[-z],maxsteps=1000,error=False)
   # Nonreal parents: one representative per lift (its conjugate is implicit).
   # Real parent: exact schema order real lift, then positive-imag representatives.
   keep=sorted(xs,key=lambda w:mp.arg(w)) if b<3 else sorted([w for w in xs if abs(mp.im(w))<mp.mpf('1e-90')]+[w for w in xs if mp.im(w)>mp.mpf('1e-90')],key=lambda w:(0 if abs(mp.im(w))<mp.mpf('1e-90') else 1,mp.arg(w)))
   for slot,w in enumerate(keep):
    cc=(fr(mp.re(w)),fr(mp.im(w)));L,R=lift_rouche_bounds(c,MR,cc,LR,S,q);assert L<R
    kind='linear' if b==3 and abs(mp.im(w))<mp.mpf('1e-90') else 'quadratic';groups.append({'index':len(groups),'block':b,'family':fam,'slot':slot,'kind':kind,'disk':dj(w,LR),'lhs':f'{L.numerator}/{L.denominator}','rhs':f'{R.numerator}/{R.denominator}'})
 assert len(groups)==29 and [g['kind'] for g in groups].count('linear')==2
 disks=[]
 for g in groups:
  c=complex_from_json(g['disk']['center']);r=parse_fraction(g['disk']['radius']);disks.append((c,r));
  if g['kind']=='quadratic':disks.append(((c[0],-c[1]),r))
 assert len(disks)==56
 for i,(c,r) in enumerate(disks):
  for d,t in disks[:i]:assert disk_disjoint(c,r,d,t)
 out={'format':'H7-a3b5-mixed-rational-rouche-v1','H':H,'a':3,'b':5,'mp_dps':110,'scale':str(S),'mother_radius':'1/10^70','lift_radius':'1/10^65','mothers':moms,'groups':groups,'coverage':{'root_disks':56,'quadratic_groups':27,'real_linear_groups':2,'pairwise_disjoint':True},'schema_order':'complex blocks 0..2: a3 then b5; real block: a linear,a quadratic,b linear,b quadratic,b quadratic'}
 Path('H7-a3b5-root-balls.json').write_text(json.dumps(out,indent=2)+'\n');print(hashlib.sha256(Path('H7-a3b5-root-balls.json').read_bytes()).hexdigest())
if __name__=='__main__':main()
