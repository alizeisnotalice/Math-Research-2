import json,sys
from pathlib import Path
def nb(v,n):
 x,y=v%n,v//n;return [q for q,o in ((v-1,x>0),(v+1,x+1<n),(v-n,y>0),(v+n,y+1<n)) if o]
def peaks(a,n):return [v for v in range(len(a))if all(a[q]<a[v]for q in nb(v,n))]
def merge(a,n):
 p=list(range(len(a)));on=[0]*len(a)
 def f(x):
  while p[x]!=x:p[x]=p[p[x]];x=p[x]
  return x
 for v in sorted(range(len(a)),key=a.__getitem__,reverse=True):
  on[v]=1;r={f(q)for q in nb(v,n)if on[q]}
  if len(r)>1:return a[v]
  for q in r:p[q]=v
def one(z):
 n=z['n'];a=z['labels'];N=n*n;assert sorted(a)==list(range(1,N+1));p=peaks(a,n);assert len(p)==2;d=abs(p[0]%n-p[1]%n)+abs(p[0]//n-p[1]//n);v=min(p,key=a.__getitem__);M=b=0
 for x in nb(v,n):
  M=max(M,a[x])
  if all(q==v or a[q]<a[x]for q in nb(x,n)):b=max(b,a[x])
 assert d==z['d'] and d>=n/2;assert (M,b,a[v])==(z['M'],z['b'],z['t']);beta=(M-b)/(a[v]-1);assert abs(beta-z['beta'])<1e-12
 # deterministic steepest-ascent basins and interface edge count
 root={}
 for w in range(N):
  path=[];u=w
  while u not in root:
   path.append(u);hi=[q for q in nb(u,n)if a[q]>a[u]]
   if not hi:root[u]=u;break
   u=max(hi,key=a.__getitem__)
  for q in path:root[q]=root[u]
 edge=sum(root[x]!=root[y]for x in range(N)for y in nb(x,n)if x<y)
 return {'n':n,'distance':d,'beta':beta,'M':M,'b':b,'t':a[v],'first_merge_label':merge(a,n),'interface_edges':edge,'basin_sizes':[sum(root[x]==q for x in range(N))for q in p]}
inp=Path(sys.argv[1]);out=Path(sys.argv[2]);rows=[one(json.loads(x))for x in inp.read_text().splitlines()if 'labels'in x];out.write_text(json.dumps({'count':len(rows),'rows':rows},indent=2)+'\n');print(json.dumps({'count':len(rows),'n':sorted(set(x['n']for x in rows))}))
