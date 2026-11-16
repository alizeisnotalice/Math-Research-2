"""H8(3,5) search suggestions: reversed monic factors, exact dyadic boxes.

This is a searcher, not an independent coverage/certificate verifier.
Each 2-byte row stores side, offset from the highest degree; 255,255 is
reserved for the standard mask. Every other unresolved case is explicit.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
import argparse, hashlib, json, os, struct, time
from root_ball_common import complex_from_json, parse_fraction, abs_bounds, polynomial_rouche_bounds, lift_rouche_bounds

H=[1,1,1,0,0,0,0,1,1]; ALL=(1<<32)-1; STD=0x07070707
COUNTS={3:625,5:81,6:10000,8:50625,9:10000,10:81,11:810000,12:625,13:50625,14:810000,15:8}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom(p,v):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def mul(p,q,D,cutoff):
 n=min(len(p)+len(q)-1,cutoff+1);lo=[0]*n;hi=[0]*n
 for i,(a,b) in enumerate(p):
  for j,(c,d) in enumerate(q):
   if i+j>=n:break
   v=(a*c,a*d,b*c,b*d);lo[i+j]+=min(v);hi[i+j]+=max(v)
 return [(a//D,-((-b)//D)) for a,b in zip(lo,hi)]
def options():
 out={}
 for u in range(4):
  for v in range(6):
   for aa in combinations(range(3),u):
    for bb in combinations(range(3,8),v):out.setdefault(5*u+3*v,[]).append(sum(1<<i for i in aa+bb))
 return {e:sorted(a) for e,a in out.items()}
def expected():
 o=options()
 assert {e:len(o[e])**4 for e in range(1,15) if e in o}=={e:n for e,n in COUNTS.items() if e<15}
 for E in range(1,16):
  choices=[o.get(E,[])]*4
  if E==15:choices=[[7]]+choices[1:]
  for b in product(*choices):yield E,sum(x<<(8*i) for i,x in enumerate(b))
def prepare(b,D):
 assert b['H']==H and b['a']==3 and b['b']==5 and len(b['groups'])==32
 S=int(b['scale']);m={x['block']:x for x in b['mothers']};fac={}
 assert set(m)==set(range(4))
 for x in m.values():
  z=complex_from_json(x['disk']['center']);r=parse_fraction(x['disk']['radius']);l,h=polynomial_rouche_bounds(H,z,r,S);assert z[1]>r>0 and l<h
 for g in b['groups']:
  i=g['index'];k=i//8;off=i%8;fam='a' if off<3 else 'b';slot=off if off<3 else off-3
  assert g['block']==k and g['family']==fam and g['slot']==slot and g['kind']=='quadratic'
  z=complex_from_json(g['disk']['center']);r=parse_fraction(g['disk']['radius']);x=m[k];c=complex_from_json(x['disk']['center']);s=parse_fraction(x['disk']['radius']);q=3 if fam=='a' else 5
  l,h=lift_rouche_bounds(c,s,z,r,S,q);assert abs(z[1])>r>0 and l<h
  l,h=abs_bounds(z,S);v=[(F(1),F(1)),(-2*(z[0]+r),-2*(z[0]-r)),(max(F(0),l-r)**2,(h+r)**2)]
  assert i not in fac;fac[i]=[(int((l*D).__floor__()),int((h*D).__ceil__())) for l,h in v]
 assert set(fac)==set(range(32))
 block=[]
 for k in range(4):
  row={0:[(D,D)]}
  for bits in range(1,256):
   bit=(bits&-bits).bit_length()-1;row[bits]=mul(row[bits^(1<<bit)],fac[8*k+bit],D,16)
  block.append(row)
 return block
def top_one(block,m):
 lo=hi=0
 for k in range(4):
  p=block[k][m>>(8*k)&255]
  if len(p)>1:l,h=p[1];lo+=l;hi+=h
 return lo,hi
def slow(block,m,D):
 degree=2*m.bit_count()
 for cutoff in (4,8,16,32,64):
  p=[(D,D)]
  for k in range(4):p=mul(p,block[k][m>>(8*k)&255],D,min(cutoff,degree))
  for j,(l,h) in enumerate(p):
   if h<0:return j
  if cutoff>=degree:break
 return None
def main():
 a=argparse.ArgumentParser();a.add_argument('--balls',type=Path,required=True);a.add_argument('--ledger',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--smoke',action='store_true');args=a.parse_args()
 assert Path('/Volumes/data').is_mount();assert not args.out.exists() or not any(args.out.iterdir());args.out.mkdir(parents=True,exist_ok=True)
 raw=args.ledger.read_bytes();assert len(raw)==4*sum(COUNTS.values())
 for n,(E,m) in enumerate(expected()):assert struct.unpack_from('<I',raw,4*n)[0]==m
 assert n+1==1742670
 D=1<<128;block=prepare(json.loads(args.balls.read_text()),D);selected=None
 if args.smoke:
  selected=set();cursor=0
  for E,count in COUNTS.items():selected.update((cursor,cursor+count-1));cursor+=count
  selected.add(next(n for n,(_,m) in enumerate(expected()) if m==STD))
 meta={'format':'H8-a3b5-top-offset-u8-v1','status':'running','bits':128,'balls_sha256':sha(args.balls),'ledger_sha256':sha(args.ledger),'searcher_sha256':sha(__file__),'root_helper_sha256':sha(Path(__file__).with_name('root_ball_common.py')),'canonical_count':1742670,'processed':0,'negative':0,'standard':0,'unknown':[],'offset_counts':{},'smoke_positions':sorted(selected) if selected else None,'pid':os.getpid()};start=time.monotonic();dest=args.out/'witnesses.u8';hist={}
 try:
  with dest.open('wb') as f:
   for n,(E,m) in enumerate(expected()):
    if selected is not None and n not in selected:continue
    found=None
    for side,use in [(0,m),(1,ALL^m)]:
     if top_one(block,use)[1]<0:found=(side,1);break
    if found is None and m!=STD:
     for side,use in [(0,m),(1,ALL^m)]:
      j=slow(block,use,D)
      if j is not None:found=(side,j);break
    if found is None:
     found=(255,255)
     if m==STD:meta['standard']+=1
     else:meta['unknown'].append({'position':n,'E':E,'mask':m})
    else:
     meta['negative']+=1;hist[str(found[1])]=hist.get(str(found[1]),0)+1
    f.write(bytes(found));meta['processed']+=1
    if meta['processed']%10000==0:
     f.flush();meta.update(elapsed_seconds=time.monotonic()-start,offset_counts=hist);atom(args.out/'checkpoint.json',meta)
  meta['status']='smoke_complete' if args.smoke else 'complete'
 except BaseException as e:meta['status']='partial_failed_or_interrupted';meta['exception']=type(e).__name__;raise
 finally:
  meta.update(elapsed_seconds=time.monotonic()-start,offset_counts=hist,witnesses_sha256=sha(dest),witness_bytes=dest.stat().st_size);atom(args.out/'checkpoint.json',meta)
  if meta['status'] in ('complete','smoke_complete'):atom(args.out/'certificate.json',meta)
  print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
