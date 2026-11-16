import itertools,json,subprocess,hashlib,os
from pathlib import Path
B=Path(__file__).parent; n=5;N=25
nb=[]
for v in range(N):
 x,y=v%n,v//n;s=set()
 for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
  if 0<=X<n and 0<=Y<n:s.add(Y*n+X)
 nb.append(s)
records=[]
for a in range(N):
 for b in nb[a]:
  if a>b:continue
  ca=nb[a]-{b};cb=nb[b]-{a}
  for P in itertools.combinations(ca|cb,3):
   S=set(P)
   if any((nb[x]&S)-{x} for x in S) or set.intersection(*(nb[x] for x in S)):continue
   for h,o,H,L in ((a,b,ca,cb),(b,a,cb,ca)):
    if len(S&H)==2 and len(S&L)==1:
     records.append((h,o,sum(1<<x for x in S)))
records=sorted(set(records));assert len(records)==172,len(records)
out=B/"n5-all-checkpoints";out.mkdir(exist_ok=True)
src=hashlib.sha256((B/"helper_dp.cpp").read_bytes()).hexdigest()
for i,(a,b,p) in enumerate(records):
 fn=out/f"{i:03d}.json"
 if fn.exists():
  try:
   z=json.loads(fn.read_text())
   if z.get("record_index")==i and z.get("source")==src:continue
  except: pass
 tmp=out/f"{i:03d}.tmp"
 q=subprocess.run([str(B/"helper_dp"),"5",str(a),str(b),str(p)],capture_output=True,text=True,check=True)
 z=json.loads(q.stdout);z.update(record_index=i,source=src,args=[a,b,p])
 tmp.write_text(json.dumps(z,separators=(',',':')));os.replace(tmp,fn)
 print(i+1,len(records),flush=True)
