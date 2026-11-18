from fractions import Fraction as F
from pathlib import Path
import sys,json,hashlib
import numpy as np
from flint import arb,ctx
ctx.prec=160
comp=Path(sys.argv[1]).resolve();treepath=Path(sys.argv[2]).resolve();scale=F(sys.argv[5]);out=Path(sys.argv[6]);sys.path.insert(0,str(comp))
from certify_polygon_083731 import ORDERS,Q,ar,verified_table
from certify_tomar_six_high import SourceTree
tree=SourceTree(treepath);source=json.loads((comp/'tomar-six-high-input.json').read_text());M=32768;target=(Q*scale)**2;pi=arb.pi();result=[]
for level in (32,64):
 path=Path(sys.argv[3 if level==32 else 4]);counts,points,stamp=verified_table(path,tree,level,source['tables'][str(level)]);offsets=np.r_[0,np.cumsum(counts)]
 for i,n in enumerate(ORDERS):
  body=points[:,offsets[i]:offsets[i+1]].astype(np.int64);assert np.abs(body).max()<Q
  if n==0:
   assert np.all(body==body[0]);d=np.square(body[0]).sum(axis=1);jidx=int(np.argmax(d));peak=4*int(d[jidx]);key=[0,0,jidx];condition='scaled-disc-radius'
  else:
   fixed=n==3
   if fixed:assert np.all(body==body[0])
   end=F(0) if fixed else F(tree.root['hi'][3*i-4]);amplitude=Q*ar(scale)/(2*(pi/(2*n)).cos());peak=-1;key=[];condition='scaled-unit-disc'
   for j in range(1 if fixed else 1025):
    for k in range(n):
     t=ar(end*j/1024)+2*pi*k/n;sy,cx=t.sin_cos();c=[amplitude*cx,amplitude*sy];whole=np.array([round(float(v.mid())) for v in c],dtype=np.int64);frac=np.array([round(float(((v-int(w))*M).mid())) for v,w in zip(c,whole)],dtype=np.int64)
     assert max(abs(whole))<Q and all(abs(v-int(w))<1 and abs((v-int(w))*M-int(f))<1 for v,w,f in zip(c,whole,frac))
     d=body[j]-whole;linear=2*(np.abs(d)-d*frac).sum(axis=1);upper=np.square(d).sum(axis=1)-((-linear)//M)+2;idx=int(np.argmax(upper));v=int(upper[idx])
     if v>peak:peak=v;key=[j,k,idx]
  gap=target-peak;result.append({'table':level,'order':n,'check':condition,'peak_integer_upper':peak,'gap_exact':str(gap),'positive':gap>0,'argmax_row_circle_point':key,'point':list(map(int,body[key[0],key[2]]))});print(level,n,'gap',str(gap),'positive',gap>0,'argmax',key,flush=True)
out.write_text(json.dumps({'scale':str(scale),'target':str(target),'precision_bits':ctx.prec,'results':result,'all_passed':all(x['positive'] for x in result),'inputs':{str(i):hashlib.sha256(Path(sys.argv[i]).read_bytes()).hexdigest() for i in (3,4)}},indent=2)+'\n')
