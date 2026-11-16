import json,hashlib,mpmath as mp
from fractions import Fraction
from pathlib import Path
from root_ball_common import qcomplex,complex_to_json,fraction_string,polynomial_rouche_bounds,lift_rouche_bounds,disk_disjoint,disk_nonreal
H=[1,1,1,0,0,0,1,1];S=10**160;MR=Fraction(1,10**70);LR=Fraction(1,10**65)
def fr(x):return Fraction(mp.nstr(x,105))
def disk(z,r):return {'center':complex_to_json((fr(mp.re(z)),fr(mp.im(z)))),'radius':fraction_string(r)}
def main():
 mp.mp.dps=110; roots=mp.polyroots(list(reversed(H)),maxsteps=1000,error=False); real=[z for z in roots if abs(mp.im(z))<mp.mpf('1e-90')]; pos=sorted([z for z in roots if mp.im(z)>0],key=lambda z:mp.re(z));assert len(real)==1 and len(pos)==3
 mothers=[];lifts=[]
 for b,z in enumerate(pos+[real[0]]):
  c=(fr(mp.re(z)),fr(mp.im(z)));L,R=polynomial_rouche_bounds(H,c,MR,S);assert L<R
  mothers.append({'block':b,'kind':'complex_mother' if b<3 else 'real_mother','disk':disk(z,MR),'lhs':fraction_string(L),'rhs':fraction_string(R)})
  for q,fam in ((3,'a'),(4,'b')):
   zs=mp.polyroots([1]+[0]*(q-1)+[-z],maxsteps=1000,error=False)
   for j,w in enumerate(zs):
    if b==3 and mp.im(w)<-mp.mpf('1e-90'): continue
    cc=(fr(mp.re(w)),fr(mp.im(w)));L,R=lift_rouche_bounds(c,MR,cc,LR,S,q);assert L<R
    typ='linear_real' if b==3 and q==3 and abs(mp.im(w))<mp.mpf('1e-90') else 'quadratic'
    lifts.append({'source_group_index':len(lifts),'block':b,'family':fam,'exponent':q,'lift_index':j,'kind':typ,'disk':disk(w,LR),'lhs':fraction_string(L),'rhs':fraction_string(R)})
 # target root disks: nonreal entries represent their conjugate pair, real one itself
 ds=[]
 for x in lifts:
  z=x['disk']['center'];from root_ball_common import complex_from_json,parse_fraction
  c=complex_from_json(z);r=parse_fraction(x['disk']['radius']);ds.append((c,r))
  if x['kind']=='quadratic':ds.append(((c[0],-c[1]),r))
 assert len(lifts)==25 and len(ds)==49
 for i,(c,r) in enumerate(ds):
  for d,t in ds[:i]:assert disk_disjoint(c,r,d,t)
 out={'format':'H7-a3b4-mixed-rational-rouche-v1','H':H,'mp_dps':110,'scale':str(S),'mothers':mothers,'lifts':lifts,'coverage':{'root_disks':49,'quadratic_groups':24,'real_linear_groups':1,'pairwise_disjoint':True}}
 Path('H7-a3b4-mixed-root-balls.json').write_text(json.dumps(out,indent=2)+'\n');print(hashlib.sha256(Path('H7-a3b4-mixed-root-balls.json').read_bytes()).hexdigest())
if __name__=='__main__':main()
