#!/usr/bin/env python3
"""Independent full endpoint wall census and four direct vertical unions, no coverage tree."""
import argparse,gzip,hashlib,json,os,sys
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
sys.dont_write_bytecode=True

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def E(e):
 a,v,s,q,b=e;a,v,q,b=map(F,(a,v,q,b));s=int(s)
 if not s:return a,v,0,F(0),F(0)
 if not b:
  n,d=isqrt(q.numerator),isqrt(q.denominator)
  if n*n==q.numerator and d*d==q.denominator:return a+s*F(n,d),v,0,F(0),F(0)
 return a,v,s,q,b

def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
 p=[x*y for x in a for y in b];return min(p),max(p)
def root(q,bits=256):
 assert q>0;D=1<<bits;n=isqrt(q.numerator*D*D//q.denominator)
 assert n*n*q.denominator<=q.numerator*D*D<(n+1)*(n+1)*q.denominator
 return F(n,D),F(n,D) if n*n*q.denominator==q.numerator*D*D else F(n+1,D)

def run(slabpath,leftpath,rightpath):
 slab=json.loads(gzip.decompress(Path(slabpath).read_bytes()));left=json.loads(Path(leftpath).read_text());right=json.loads(Path(rightpath).read_text());assert len(slab['boxes'])==26327
 assert sha(slabpath)==left['slab_file_sha256']==right['slab_file_sha256']
 assert left['wall_t_bounds']==right['wall_t_bounds'] and left['wall_equation']==right['wall_equation']
 eq=left['wall_equation'];A,B=E(eq['A_endpoint']),E(eq['B_endpoint']);t0=F(slab['t0']);bounds=list(map(F,left['wall_t_bounds']))
 # Independently derive a root enclosure directly from A/B wall equation.
 C=A[0]-B[0];assert A[2]==-1 and A[4]==0 and B[2]==1 and B[4]<0
 ra=root(A[3],336);low,high=C-ra[1],C-ra[0];assert low>0
 w=sorted(((low*low-B[3])/B[4]+t0,(high*high-B[3])/B[4]+t0))
 assert w[0]<=bounds[1] and bounds[0]<=w[1]
 # Use union of independent and producer enclosures, rounded outwards to 256 dyadics.
 D=1<<256;tl=min(w[0],bounds[0]);th=max(w[1],bounds[1]);tl=F(tl.numerator*D//tl.denominator,D);th=F(-((-th.numerator*D)//th.denominator),D)
 @lru_cache(None)
 def val(e):
  a,v,s,q,b=e
  if not s:return a,a
  lo,hi=sorted((q+b*(tl-t0),q+b*(th-t0)));assert lo>0
  rl,rh=root(lo)[0],root(hi)[1]
  return (a+rl,a+rh) if s>0 else (a-rh,a-rl)
 @lru_cache(None)
 def der(e):
  a,v,s,q,b=e
  if not s:return F(0),F(0)
  lo,hi=sorted((q+b*(tl-t0),q+b*(th-t0)));rl,rh=root(lo)[0],root(hi)[1];c=s*b/2
  return (c/rh,c/rl) if c>=0 else (c/rl,c/rh)
 boxes=[tuple(map(E,b)) for b in slab['boxes']];fl,fr,fb,ft=map(F,slab['domain']);constant=lambda a:(a,F(0),0,F(0),F(0));X={constant(fl),constant(fr)};Y={constant(fb),constant(ft)}
 for l,r,a,b in boxes:X.update((l,r));Y.update((a,b))
 def ordered(es,axis):
  seq=sorted(es,key=lambda e:sum(val(e)));exceptions=[]
  for a,b in zip(seq,seq[1:]):
   if val(a)[1]>=val(b)[0]:
    assert axis=='x' and {a,b}=={A,B},('unresolved endpoint equality',axis,a,b)
    exceptions.append((a,b))
  assert len(exceptions)==(1 if axis=='x' else 0)
  return seq
 xs=ordered(X,'x');ys=ordered(Y,'y');yrank={e:i for i,e in enumerate(ys)};yvals=[val(e) for e in ys];xwall=val(A)
 # Independent mass-only jet/order: A_t=0 and B_t<0 => A-B increasing.
 da,db=der(A),der(B);assert da==(0,0) and db[1]<0
 def relative(e):
  if e==A:return 'A'
  if e==B:return 'B'
  v=val(e)
  if v[1]<xwall[0]:return -1
  assert v[0]>xwall[1],('noncritical endpoint touching critical wall',e)
  return 1
 C0=set();events={'A':(set(),set()),'B':(set(),set())};ranges=[]
 for i,(l,r,a,b) in enumerate(boxes):
  ll,rr=relative(l),relative(r);bottom,top=max(yrank[a],yrank[constant(fb)]),min(yrank[b],yrank[constant(ft)]);assert bottom<top;ranges.append((bottom,top))
  if ll==-1 and rr in (1,'A','B'):C0.add(i)
  for e,kind in [(ll,0),(rr,1)]:
   if e in ('A','B'):events[e][kind].add(i)
 def apply(S,event):
  starts,ends=events[event];assert not starts&S and ends<=S
  return (S|starts)-ends
 CA=apply(C0,'A');CB=apply(C0,'B');C1=apply(CA,'B');assert C1==apply(CB,'A')
 def union(S):
  ranges_sorted=sorted(ranges[i] for i in S);merged=[]
  for a,b in ranges_sorted:
   if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
   else:merged.append((a,b))
  total=(F(0),F(0))
  for a,b in merged:total=add(total,sub(yvals[b],yvals[a]))
  return total,len(merged)
 lengths=[union(S) for S in [C0,CA,CB,C1]];L0,LA,LB,L1=[r[0] for r in lengths];coef=sub(add(LA,LB),add(L0,L1));jump=mul(sub(da,db),coef);betajump=(-coef[1],-coef[0]);assert jump[1]<0 and betajump[0]>0
 directjump=sub(tuple(map(F,right['A_t_wall_bounds'])),tuple(map(F,left['A_t_wall_bounds'])));directbeta=sub(tuple(map(F,right['beta_wall_bounds'])),tuple(map(F,left['beta_wall_bounds'])))
 assert jump[0]<=directjump[1] and directjump[0]<=jump[1]
 assert betajump[0]<=directbeta[1] and directbeta[0]<=betajump[1]
 return {'status':'passed','scope':'complete raw endpoint wall census and four direct vertical union states; no full x coverage-tree replay or omitted bbox rescan','slab_sha256':sha(slabpath),'left_output_sha256':sha(leftpath),'right_output_sha256':sha(rightpath),'candidate_count':len(boxes),'x_endpoint_count':len(X),'y_endpoint_count':len(Y),'nonpersistent_x_events':1,'nonpersistent_y_events':0,'persistent_endpoint_duplicates_grouped_by_full_expression':True,'critical_groups':{k:{'starts':len(v[0]),'ends':len(v[1])} for k,v in events.items()},'state_active_counts':[len(S) for S in [C0,CA,CB,C1]],'state_merged_counts':[v[1] for v in lengths],'four_union_length_bounds':[[str(q) for q in v[0]] for v in lengths],'C_bounds':list(map(str,coef)),'A_t_jump_bounds':list(map(str,jump)),'beta_jump_bounds':list(map(str,betajump)),'overlaps_full_sweep_derivative_intervals':True,'wall_t_bounds_independent':list(map(str,w)),'all_radicands_positive':True,'A_t_jump_strict_negative':True,'beta_jump_strict_positive':True,'audit_code_sha256':sha(__file__),'python':sys.version,'executable':sys.executable}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--slab',required=True);p.add_argument('--left',required=True);p.add_argument('--right',required=True);p.add_argument('--output',required=True);a=p.parse_args();assert os.path.ismount('/Volumes/data');out=run(a.slab,a.left,a.right);Path(a.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['four_union_length_bounds','wall_t_bounds_independent','A_t_jump_bounds','beta_jump_bounds','C_bounds']},indent=2));print('jump_display',[float(F(q)) for q in out['A_t_jump_bounds']]);print('beta_display',[float(F(q)) for q in out['beta_jump_bounds']])
