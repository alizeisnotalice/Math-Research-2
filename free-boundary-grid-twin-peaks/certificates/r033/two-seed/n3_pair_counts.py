import itertools,json
n=3;N=9
def adj(v,w):return abs(v%n-w%n)+abs(v//n-w//n)==1
out={(p,r):0 for p in range(N) for r in range(N) if p!=r and not adj(p,r)}
for w in itertools.permutations(range(N)):
 p,r=w[:2]
 if (p,r) not in out:continue
 seen={p,r};ok=True
 for v in w[2:]:
  if not any(adj(v,u) for u in seen):ok=False;break
  seen.add(v)
 if ok:out[p,r]+=1
rows=[{'p':p,'r':r,'d':abs(p%n-r%n)+abs(p//n-r//n),'B':v} for (p,r),v in out.items()]
print(json.dumps({'n':n,'pairs':rows,'min_by_distance':{d:min(x['B'] for x in rows if x['d']==d) for d in range(2,5)},'max_by_distance':{d:max(x['B'] for x in rows if x['d']==d) for d in range(2,5)}},indent=2))
