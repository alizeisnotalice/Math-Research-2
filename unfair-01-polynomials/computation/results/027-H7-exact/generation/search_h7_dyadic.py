"""H7 mixed-group searcher: outward dyadic intervals; separate verifier required."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
import argparse,json,hashlib,time,os
from root_ball_common import complex_from_json,parse_fraction,abs_bounds,polynomial_rouche_bounds,lift_rouche_bounds

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom(p,x):
 tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
def floorf(x):return x.numerator//x.denominator
def ceilf(x):return -((-x.numerator)//x.denominator)
def mul(p,q,D):
 a=[0]*(len(p)+len(q)-1);b=[0]*len(a)
 for i,(l,h) in enumerate(p):
  for j,(s,t) in enumerate(q):
   vals=[l*s,l*t,h*s,h*t];a[i+j]+=min(vals);b[i+j]+=max(vals)
 return [(l//D,-((-h)//D)) for l,h in zip(a,b)]
def prepare(b,bits):
 S=int(b['scale']);D=1<<bits;M={m['block']:m for m in b['mothers']};f={}
 for m in M.values():
  z=complex_from_json(m['disk']['center']);r=parse_fraction(m['disk']['radius']);assert r>0
  l,h=polynomial_rouche_bounds(b['H'],z,r,S);assert l<h
 for g in b['lifts']:
  z=complex_from_json(g['disk']['center']);r=parse_fraction(g['disk']['radius']);m=M[g['block']];c=complex_from_json(m['disk']['center']);t=parse_fraction(m['disk']['radius'])
  l,h=lift_rouche_bounds(c,t,z,r,S,g['exponent']);assert l<h
  if g['kind']=='linear_real':assert z[1]==0 and z[0]+r<0;iv=[(-z[0]-r,-z[0]+r),(F(1),F(1))]
  else:
   assert abs(z[1])>r;l,h=abs_bounds(z,S);iv=[(max(F(0),l-r)**2,(h+r)**2),(-2*(z[0]+r),-2*(z[0]-r)),(F(1),F(1))]
  i=g['source_group_index'];assert i not in f;f[i]=[(floorf(l*D),ceilf(h*D)) for l,h in iv]
 return f

def build(b):
 by={(k,f):[g['source_group_index'] for g in b['lifts'] if g['block']==k and g['family']==f] for k in range(4) for f in ['a','b']}
 lin=next(g['source_group_index'] for g in b['lifts'] if g['kind']=='linear_real');aq=next(i for i in by[3,'a'] if i!=lin);all_a=sum(1<<i for k in range(4) for i in by[k,'a']);tasks=[]
 for E in [4,6,8,10,12]:
  choices=[]
  for k in range(3):
   rows=[]
   for u in range(4):
    for v in range(5):
     if 4*u+3*v!=E:continue
     for aa in combinations(by[k,'a'],u):
      for bb in combinations(by[k,'b'],v):rows.append(sum(1<<i for i in aa+bb))
   if E==12 and k==0:rows=[sum(1<<i for i in by[0,'a'])]
   choices.append(rows)
  real=[]
  for eps in [0,1]:
   for eta in [0,1]:
    for w in range(3):
     if 4*eps+8*eta+6*w!=E:continue
     for ids in combinations(by[3,'b'],w):real.append(sum(1<<i for i in ids)+(eps<<lin)+(eta<<aq))
  for parts in product(*choices,real):tasks.append((E,parts[0]|parts[1]|parts[2]|parts[3]))
 assert len(tasks)==len(set(tasks))==12158;return sorted(tasks),all_a

def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--parent-masks',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--bits',type=int,default=128);p.add_argument('--smoke',action='store_true');a=p.parse_args();assert Path('/Volumes/data').is_mount() and 32<=a.bits<=512
 assert not a.out.exists() or not any(a.out.iterdir());a.out.mkdir(parents=True,exist_ok=True);b=json.loads(a.balls.read_text());f=prepare(b,a.bits);tasks,standard=build(b);assert tasks==[tuple(x) for x in json.loads(a.parent_masks.read_text())];D=1<<a.bits
 if a.smoke:tasks=[next(t for t in tasks if t[0]==E) for E in [4,6,8,10]]+[(12,standard),next(t for t in tasks if t[0]==12 and t[1]!=standard)]
 rec={'status':'running','pid':os.getpid(),'bits':a.bits,'arithmetic':'integer directed dyadic interval convolution','balls_sha256':sha(a.balls),'source_sha256':sha(__file__),'parent_masks_sha256':sha(a.parent_masks),'complete_necessary_mask_set_matched':12158,'requested':len(tasks),'processed':0,'negative_certificates':[],'unknown':[]};start=time.monotonic()
 try:
  for E,mask in tasks:
   found=False
   for side,use in [('A',mask),('B',((1<<25)-1)^mask)]:
    q=[(D,D)]
    for i in range(25):
     if use>>i&1:q=mul(q,f[i],D)
    bad=next(((k,l,h) for k,(l,h) in enumerate(q) if h<0),None)
    if bad:
     k,l,h=bad;rec['negative_certificates'].append({'E':E,'mask':mask,'side':side,'index':k,'lo_n':str(l),'hi_n':str(h)});found=True;break
   if not found:rec['unknown'].append({'E':E,'mask':mask,'standard':mask==standard})
   rec['processed']+=1
   if rec['processed']%100==0:rec['elapsed_seconds']=time.monotonic()-start;atom(a.out/'checkpoint.json',rec)
  rec['status']='smoke_complete' if a.smoke else 'complete'
 except BaseException as e:rec['status']='partial_failed_or_interrupted';rec['exception']=type(e).__name__;raise
 finally:
  rec['elapsed_seconds']=time.monotonic()-start;atom(a.out/'checkpoint.json',rec)
  if rec['status'] in ('complete','smoke_complete'):atom(a.out/'certificate.json',rec)
  print(json.dumps({k:rec[k] for k in ['status','pid','processed','elapsed_seconds']}|{'negative':len(rec['negative_certificates']),'unknown':len(rec['unknown'])}),flush=True)
if __name__=='__main__':main()
