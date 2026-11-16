import importlib.util,itertools,json,hashlib,random,math,pathlib,time
s=pathlib.Path('/Volumes/data/CodexTemp/T2/a022-tree-top2-scaling-r024/tree_peak_dp.py')
spec=importlib.util.spec_from_file_location('tree',s);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sha=hashlib.sha256(s.read_bytes()).hexdigest(); rng=random.Random(240116);rows=[];start=time.time()
def graph(n):
 seq=[rng.randrange(n) for _ in range(n-2)];deg=[1]*n
 for v in seq:deg[v]+=1
 adj=[set() for _ in range(n)]
 for v in seq:
  u=next(i for i in range(n) if deg[i]==1);adj[u].add(v);adj[v].add(u);deg[u]-=1;deg[v]-=1
 u,v=[i for i in range(n) if deg[i]==1];adj[u].add(v);adj[v].add(u);return adj
for n in range(3,9):
 for trial in range(2):
  adj=graph(n);a1=a2=b=0
  for lab in itertools.permutations(range(n)):
   peaks=[v for v in range(n) if all(lab[v]>lab[w] for w in adj[v])]
   if len(peaks)==1:a1+=1
   if len(peaks)==2:
    a2+=1;b+=sorted(lab[v] for v in peaks)==[n-2,n-1]
  aa,bb,poly=m.counts(adj);assert(aa,bb)==(a2,b)
  assert poly[1]*math.factorial(n)==a1
  for root in [0,n-1]:assert m.peak_probability_poly(adj,root)==poly
  hook=0
  for root in range(n):
   sizes={}
   def walk(v,p):
    sizes[v]=1+sum(walk(w,v) for w in adj[v] if w!=p);return sizes[v]
   walk(root,-1);h=math.factorial(n-1)
   den=math.prod(v for k,v in sizes.items() if k!=root);assert h%den==0;hook+=h//den
  assert hook==a1
  rows.append({'n':n,'trial':trial,'edges':[[v,w] for v in range(n) for w in adj[v] if v<w],'A1':a1,'A2':a2,'Btop2':b})
assert hashlib.sha256(s.read_bytes()).hexdigest()==sha
out={'status':'passed','method':'independent all-label permutations and tree hook-length A1, plus message root invariance','seed':240116,'source_sha256':sha,'cases':rows,'wall_seconds':time.time()-start}
pathlib.Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
