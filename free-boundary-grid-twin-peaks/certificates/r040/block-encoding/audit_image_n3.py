from itertools import permutations
from fractions import Fraction as F
from math import comb
import json
n=3;N=9;ALL=(1<<N)-1
adj=[]
for v in range(N):
 x,y=v%n,v//n;adj.append(sum(1<<(Y*n+X) for X,Y in [(x-1,y),(x+1,y),(x,y-1),(x,y+1)] if 0<=X<n and 0<=Y<n))
# Independent last-vertex recursion for all induced rooted one-peak counts.
a=[[0]*N for _ in range(1<<N)]
for mask in sorted(range(1,1<<N),key=int.bit_count):
 if mask&(mask-1)==0:a[mask][mask.bit_length()-1]=1;continue
 for root in range(N):
  if not(mask>>root&1):continue
  for last in range(N):
   if last!=root and mask>>last&1 and adj[last]&(mask^(1<<last)):a[mask][root]+=a[mask^(1<<last)][root]
D=[[ (r%n-p%n)**2+(r//n-p//n)**2 for p in range(N)] for r in range(N)]
far=[[abs(r%n-p%n)+abs(r//n-p//n)>n for p in range(N)] for r in range(N)]
direct={};images={};A1=0;A2=0;Iwords=0
for word in permutations(range(N)):
 seen=0;roots=[]
 for v in word:
  if not(adj[v]&seen):roots.append(v)
  seen|=1<<v
 if len(roots)==2:
  A2+=1;comps=[1<<word[0]]
  for q,v in enumerate(word[1:],1):
   hit=[i for i,c in enumerate(comps) if c&adj[v]]
   if len(hit)==2:
    low=next(c for c in comps if not(c&(1<<word[0])));b=low.bit_count();key=(q,b)
    anchor=adj[v]&low;x=(anchor&-anchor).bit_length()-1
    row=direct.setdefault(key,[0,0,0,0]);row[0]+=1;row[1]+=D[roots[0]][roots[1]];row[2]+=far[roots[0]][roots[1]];row[3]+=D[roots[0]][x];break
   if not hit:comps.append(1<<v)
   else:comps[hit[0]]|=1<<v
 if len(roots)!=1:continue
 A1+=1;r=word[0];pm=[0];unions=[0]
 for v in word:pm.append(pm[-1]|(1<<v));unions.append(unions[-1]|adj[v])
 for q in range(2,N):
  for b in range(1,q):
   i=q-b+1;B=pm[q+1]^pm[i];w=word[i-1];x=word[i]
   if unions[i-1]&B:continue
   anchor=adj[w]&B
   if not anchor or x!=(anchor&-anchor).bit_length()-1:continue
   seenB=1<<x;ok=True
   for v in word[i+1:q+1]:
    if not(adj[v]&seenB):ok=False;break
    seenB|=1<<v
   if not ok:continue
   den=a[B][x];assert den>0 and sum(a[B])<=2**(b-1)*den
   row=images.setdefault((q,b),[F(),F(),F(),F()]);row[0]+=F(sum(a[B]),den);row[1]+=F(sum(a[B][p]*D[r][p] for p in range(N)),den);row[2]+=F(sum(a[B][p]*far[r][p] for p in range(N)),den);row[3]+=F(sum(a[B])*D[r][x],den);Iwords+=1
   for p in range(N):
    if B>>p&1:assert (D[r][p]-D[r][x])**2<=8*(n-1)**2*(b-1)**2
assert A1==29568==sum(a[ALL]) and A2==143808
rows={}
for q in range(2,N):
 for b in range(1,q):
  key=(q,b);left=direct.get(key,[0,0,0,0]);im=images.get(key,[F(),F(),F(),F()]);right=[v*comb(q-1,b) for v in im];assert right==left,(key,left,right)
  rows['q=%d,b=%d'%key]=dict(direct_count=left[0],direct_euclidean_squared_sum=left[1],direct_Manhattan_far_count=left[2],direct_entry_euclidean_squared_sum=left[3],image_root_ratio_sums=[str(x) for x in im],all_exact_identities_pass=True)
for k in range(1,N-1):
 subset=[v for (q,b),v in direct.items() if b<=k];count=sum(v[0] for v in subset);diff=sum(v[1]-v[3] for v in subset);assert diff*diff<=8*(n-1)**2*(k-1)**2*count*count
Dsum=sum(r[1] for r in direct.values());Fsum=sum(r[2] for r in direct.values())
out=dict(A1=A1,A2=A2,image_words_across_parameters=Iwords,R=str(F(A2,A1)),Q=str(F(Dsum,N*A2)),total_euclidean_squared_sum=Dsum,Manhattan_far_threshold='rho>3',Manhattan_far_count=Fsum,rows=rows)
print(json.dumps(out,indent=2))
