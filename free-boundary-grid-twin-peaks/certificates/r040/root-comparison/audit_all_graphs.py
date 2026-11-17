import itertools,json
rows=[]
for n in range(1,6):
 edges=list(itertools.combinations(range(n),2));graphs=words=promotions=twopeak_words=0
 for em in range(1<<len(edges)):
  adj=[0]*n
  for k,(u,v) in enumerate(edges):
   if em>>k&1:adj[u]|=1<<v;adj[v]|=1<<u
  seen=1
  while True:
   nxt=seen
   for v in range(n):
    if seen>>v&1:nxt|=adj[v]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:continue
  graphs+=1;counts=[0]*n;twopeak=0;buckets={}
  for w in itertools.permutations(range(n)):
   mask=1<<w[0];births=1
   for v in w[1:]:
    if not(adj[v]&mask):births+=1
    mask|=1<<v
   if births==2:
    twopeak+=1;twopeak_words+=1;comps=[1<<w[0]]
    for v in w[1:]:
     touched=[i for i,c in enumerate(comps) if c&adj[v]]
     if len(touched)==2:
      q=sum(c.bit_count() for c in comps);low=next(c for c in comps if not(c&(1<<w[0]))).bit_count();buckets[q,low]=buckets.get((q,low),0)+1;break
     if not touched:comps.append(1<<v)
     else:comps[touched[0]]|=1<<v
   if births!=1:continue
   words+=1;counts[w[0]]+=1;pos={v:i for i,v in enumerate(w)}
   for x in range(n):
    chain=[x]
    while chain[-1]!=w[0]:
     v=chain[-1];chain.append(min((u for u in range(n) if adj[v]>>u&1 and pos[u]<pos[v]),key=lambda u:pos[u]))
    path=set(chain);eta=chain+[v for v in w if v not in path];P={pos[v] for v in chain};assert 0 in P
    mask=1<<eta[0]
    for v in eta[1:]:assert adj[v]&mask;mask|=1<<v
    d=len(P);rev=iter(reversed(eta[:d]));other=iter(eta[d:]);decoded=tuple(next(rev) if i in P else next(other) for i in range(n));assert decoded==w
    promotions+=1
  A1=sum(counts)
  assert A1<=2**(n-1)*min(counts)
  from math import comb
  for (q,low),value in buckets.items():assert value<=A1*2**(low-1)*comb(q-1,low)
  assert 4*twopeak<=A1*(3**(n-1)-2*n+1)
 rows.append(dict(n=n,connected_labelled_graphs=graphs,onepeak_words=words,promotion_left_inverses=promotions,twopeak_words=twopeak_words,all_root_and_branch_bounds_pass=True))
# Sharpness on paths: a(P_n)=2^(n-1), endpoint count1.
for n in range(1,17):
 from math import comb
 assert sum(comb(n-1,j) for j in range(n))==2**(n-1)
print(json.dumps(rows,indent=2))
