"""Independent searcher for H7(3,5): 128bit endpoint intervals, compact witnesses."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
import argparse,json,hashlib,time,os
from root_ball_common import complex_from_json,parse_fraction,abs_bounds,polynomial_rouche_bounds,lift_rouche_bounds

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,separators=(',',':'))+'\n');t.replace(p)
def mul(p,q,D):
 lo=[0]*(len(p)+len(q)-1);hi=[0]*len(lo)
 for i,(a,b) in enumerate(p):
  for j,(c,d) in enumerate(q):
   v=(a*c,a*d,b*c,b*d);lo[i+j]+=min(v);hi[i+j]+=max(v)
 return [(a//D,-((-b)//D)) for a,b in zip(lo,hi)]
def prepare(b,D):
 scale=int(b['scale']);m={x['block']:x for x in b['mothers']};fac={}
 for x in m.values():
  z=complex_from_json(x['disk']['center']);r=parse_fraction(x['disk']['radius']);l,h=polynomial_rouche_bounds(b['H'],z,r,scale);assert r>0 and l<h
 for g in b['groups']:
  z=complex_from_json(g['disk']['center']);r=parse_fraction(g['disk']['radius']);x=m[g['block']];c=complex_from_json(x['disk']['center']);s=parse_fraction(x['disk']['radius']);q=3 if g['family']=='a' else 5;l,h=lift_rouche_bounds(c,s,z,r,scale,q);assert r>0 and l<h
  if g['kind']=='linear':assert z[1]==0 and z[0]+r<0;v=[(-z[0]-r,-z[0]+r),(F(1),F(1))]
  else:
   assert abs(z[1])>r;l,h=abs_bounds(z,scale);v=[(max(F(0),l-r)**2,(h+r)**2),(-2*(z[0]+r),-2*(z[0]-r)),(F(1),F(1))]
  i=g['index'];assert i not in fac;fac[i]=[(int((l*D).__floor__()),int((h*D).__ceil__())) for l,h in v]
 return fac

def build(b):
 by={(k,f):[g['index'] for g in b['groups'] if g['block']==k and g['family']==f] for k in range(4) for f in ('a','b')};a=[g['index'] for g in b['groups'] if g['block']==3 and g['family']=='a' and g['kind']=='linear'][0];aq=next(i for i in by[3,'a'] if i!=a);bl=[g['index'] for g in b['groups'] if g['block']==3 and g['family']=='b' and g['kind']=='linear'][0];bq=[i for i in by[3,'b'] if i!=bl];out=[]
 for E in range(1,16):
  choices=[]
  for k in range(3):
   masks=[]
   for u in range(4):
    for v in range(6):
     if 5*u+3*v!=E:continue
     for aa in combinations(by[k,'a'],u):
      for bb in combinations(by[k,'b'],v):masks.append(sum(1<<i for i in aa+bb))
   if E==15 and k==0:masks=[sum(1<<i for i in by[0,'a'])]
   choices.append(masks)
  real=[]
  for ea,eta,eb,w in product(range(2),range(2),range(2),range(3)):
   if 5*ea+10*eta+3*eb+6*w!=E:continue
   for ids in combinations(bq,w):real.append((ea<<a)+(eta<<aq)+(eb<<bl)+sum(1<<i for i in ids))
  out.extend((E,x|y|z|w) for x,y,z,w in product(*choices,real))
 assert len(out)==len(set(out))==119062;return sorted(out),sum(1<<g['index'] for g in b['groups'] if g['family']=='a')
def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--parent-masks',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--smoke',action='store_true');a=p.parse_args();assert Path('/Volumes/data').is_mount();assert not a.out.exists() or not any(a.out.iterdir());a.out.mkdir(parents=True,exist_ok=True)
 b=json.loads(a.balls.read_text());tasks,standard=build(b);assert tasks==[tuple(x) for x in json.loads(a.parent_masks.read_text())];digest=hashlib.sha256(json.dumps(tasks,separators=(',',':')).encode()).hexdigest();D=1<<128;fac=prepare(b,D)
 if a.smoke:
  tasks=[next(t for t in tasks if t[0]==E) for E in sorted({E for E,m in tasks}) if E<15]+[(15,standard),next(t for t in tasks if t[0]==15 and t[1]!=standard)]
 rec={'format':'H7-a3b5-compact-ball-witness-v1','status':'running','pid':os.getpid(),'bits':128,'balls_sha256':sha(a.balls),'canonical_sha256':digest,'search_source_sha256':sha(__file__),'parent_mask_file_sha256':sha(a.parent_masks),'arithmetic':'integer directed 128bit endpoint intervals','complete_necessary_mask_set_matched':119062,'requested':len(tasks),'processed':0,'witnesses':[],'unknown':[]};start=time.monotonic();negative=0;minimum=None
 try:
  for E,mask in tasks:
   witness=None
   for side,use in [(0,mask),(1,((1<<29)-1)^mask)]:
    q=[(D,D)]
    for i in range(29):
     if use>>i&1:q=mul(q,fac[i],D)
    v=next(((k,h) for k,(l,h) in enumerate(q) if h<0),None)
    if v:
     k,h=v;witness=[side,k];negative+=1;minimum=min(minimum,-h) if minimum is not None else -h;break
   if witness is None:
    witness=[-1,-1];rec['unknown'].append({'position':rec['processed'],'E':E,'mask':mask,'standard':mask==standard})
   rec['witnesses'].append(witness);rec['processed']+=1
   if rec['processed']%1000==0:rec['elapsed_seconds']=time.monotonic()-start;atom(a.out/'checkpoint.json',rec)
  rec['status']='smoke_complete' if a.smoke else 'complete'
 except BaseException as e:rec['status']='partial_failed_or_interrupted';rec['exception']=type(e).__name__;raise
 finally:
  rec.update(elapsed_seconds=time.monotonic()-start,negative=negative,minimum_search_negative_margin=str(F(minimum,D)) if minimum else None);atom(a.out/'checkpoint.json',rec)
  if rec['status'] in ('complete','smoke_complete'):atom(a.out/'certificate.json',rec)
  print(json.dumps({k:rec[k] for k in ['status','pid','processed','negative','elapsed_seconds']}|{'unknown':rec['unknown']}),flush=True)
if __name__=='__main__':main()
