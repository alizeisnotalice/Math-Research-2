#!/usr/bin/env python3
"""Small exact n=5 prefix and sampled-word structural gate."""
from functools import lru_cache
import hashlib, json
from pathlib import Path
import subprocess, sys

binary, reference, out = map(Path,sys.argv[1:4])
if out.exists(): raise SystemExit("fresh output directory required")
out.mkdir(parents=True)
n=5
adj=[sum(1<<u for u in range(n*n) if abs(u%n-v%n)+abs(u//n-v//n)==1) for v in range(n*n)]
def region(h):return sum(1<<(y*n+x) for y,c in enumerate(h) for x in range(int(c)))
def neigh(S,H):
    z=0
    while S:
        b=S&-S;S-=b;z|=adj[b.bit_length()-1]
    return z&H
def safe_roots(H, side):
    A=0
    for y in range(n):
      for x in range(n):
        if (side==0 and x>=2) or (side==1 and x<3): continue
        v=y*n+x
        if H>>v&1 and not(adj[v]&~H):A|=1<<v
    return A
def direct(H,S):
  @lru_cache(None)
  def f(T):
    if T==H:return 1
    B=neigh(T,H)&~T
    return sum(f(T|(1<<v)) for v in range(n*n) if B>>v&1)
  return f(S)
def valid_event(word,h):
  H=region(h); V=(1<<(n*n))-1; K=V^H; position={v:i for i,v in enumerate(word)}
  for G,A in ((H,safe_roots(H,0)),(K,safe_roots(K,1))):
    peaks=[v for v in range(n*n) if G>>v&1 and all(position[v]<position[u] for u in range(n*n) if (adj[v]>>u&1) and (G>>u&1))]
    if len(peaks)!=1 or not(A>>peaks[0]&1):return False
  peaks=[v for v in range(n*n) if all(position[v]<position[u] for u in range(n*n) if adj[v]>>u&1)]
  return len(peaks)==2
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

ref={r['curve']:r for r in json.loads(reference.read_text())['rows']}
cases=[('22222',1),('22232',1),('23232',1|(1<<5)),('33333',1)]
receipts=[]
for i,(h,S) in enumerate(cases):
  target=out/f'{h}-{S:x}.json'
  p=subprocess.Popen([str(binary),'5',h,str(target),'--prefix-mask',f'{S:x}','--samples','5','--seed',str(9100+i)])
  rc=p.wait(); got=json.loads(target.read_text()) if rc==0 else None
  H=region(h); want=direct(H,S); row=ref[h]
  ok=(got is not None and int(got['prefix_completion'])==want and got['left_rooted_count']==row['left_rooted_count'] and got['right_rooted_count']==row['right_rooted_count'] and got['count']==row['count'] and all(valid_event(w,h) for w in got['samples']))
  receipts.append({'curve':h,'prefix_mask':f'{S:x}','pid':p.pid,'returncode':rc,'waited_reaped':True,'matched':ok,'output_sha256':sha(target) if target.exists() else None})
  if not ok:raise SystemExit(f'n5 gate mismatch {h} {S:x}')
audit={'schema':'exact-mixture-n5-small-prefix-gate-v1','cases':len(cases),'direct_prefix_completion_checked':True,'sampled_words_structurally_checked':True,'binary_sha256':sha(binary),'reference_sha256':sha(reference),'receipts':receipts}
(out/'audit.json').write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in audit.items() if k!='receipts'},sort_keys=True))
