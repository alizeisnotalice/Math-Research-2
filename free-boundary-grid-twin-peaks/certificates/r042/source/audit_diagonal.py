"""Independent n5 single-event counts by induced-region subset DP and shuffle."""
import json,math,sys
from pathlib import Path
p=Path(__file__).resolve().parent.parent;lipschitz='--lipschitz' in sys.argv
d=json.loads((p/('n5-lipschitz-ledger.json' if lipschitz else 'n5-ledger.json')).read_text());n=d['n'];N=n*n
V=set(range(N));adj=[{u for u in V if abs(u%n-v%n)+abs(u//n-v//n)==1} for v in V]
diagonal={j['curves'][0]:int(j['count']) for j in d['jobs'] if j['curves'][0]==j['curves'][1]}
def rooted_total(H,allowed):
 if not allowed:return 0
 vertices=sorted(H);idx={v:i for i,v in enumerate(vertices)};m=len(vertices)
 nb=[sum(1<<idx[u] for u in adj[v]&H) for v in vertices]
 f=[0]*(1<<m)
 for v in allowed:f[1<<idx[v]]=1
 for S in range(1,1<<m):
  if S&(S-1)==0:continue
  q=S;z=0
  while q:
   bit=q&-q;q^=bit;T=S^bit;i=bit.bit_length()-1
   if nb[i]&T:z+=f[T]
  f[S]=z
 return f[-1]
rows=[]
for h in d['curves']:
 L={v for v in V if v%n<int(h[v//n])};R=V-L;band=(n+3)//4
 AL={v for v in L if v%n<band and adj[v]<=L};AR={v for v in R if v%n>=n-band and adj[v]<=R}
 a=rooted_total(L,AL);b=rooted_total(R,AR);value=math.comb(N,len(L))*a*b
 # The four graph symmetries reduce the saved diagonal ledger.
 comp=''.join(str(n-int(x)) for x in h);key=min(h,h[::-1],comp,comp[::-1])
 assert value==diagonal[key],h
 rows.append(dict(curve=h,left_vertices=len(L),left_rooted_count=str(a),right_rooted_count=str(b),count=str(value)))
print(json.dumps(dict(independent_single_events=len(rows),S1=str(sum(int(r['count']) for r in rows)),method='induced-region integer DP times binomial shuffle; no multi-cut full-grid DP',rows=rows),indent=2))
