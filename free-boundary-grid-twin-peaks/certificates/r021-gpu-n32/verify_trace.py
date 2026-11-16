import json, math, sys

M=(1<<64)-1;S=0x9e3779b97f4a7c15;A=0xbf58476d1ce4e5b9;B=0x94d049bb133111eb
def mix(x):
 x=((x^(x>>30))*A)&M;x=((x^(x>>27))*B)&M;return x^(x>>31)
def draw(x,m):
 assert m>0
 t=(-m)&M;t%=m
 while 1:
  x=(x+S)&M;z=mix(x)
  if z>=t:return z%m,x
def key(s,r,c,d):return s^(r<<32)^((c*S)&M)^((d*0xd1b54a32d192ed03)&M)
def rebuild(w,n,t):
 seen=set();bd=set();roots=[]
 for q,v in enumerate(w,1):
  if v in seen:return None
  birth=v not in bd
  if birth:roots.append(v)
  if birth != (q==1 or q==t):return None
  seen.add(v);x=v%n;y=v//n
  for z in ([v-1]if x else [])+([v+1]if x+1<n else[])+([v-n]if y else[])+([v+n]if y+1<n else[]):
   if z not in seen:bd.add(z)
  bd.discard(v)
 return roots,seen,bd
def allowed(w,n,t,r):
 _,seen,bd=rebuild(w,n,t);free=set(range(n*n))-seen
 return sorted(free if r==1 else (free-bd if r==t else bd))
def close(a,b):
 a=float(a);b=float(b)
 if a==0.0 or b==0.0:return a==b
 return abs(a-b)<=5e-12*max(abs(a),abs(b))

path=sys.argv[1]
x=json.load(open(path));n=x['n'];t=x['tau'];P=x['P'];L=x['L'];seed=int(x['seed_hex'],16);N=n*n
assert x['implementation']=='r021-n32-rank-parallel-terminal-inverse-fibre' and x['method']=='strat'
assert x['base_source_sha256']=='6aa1f829dbe97c141781c7440c0f99ef62fa966cbe2c665fab7432a48b2ec9f8'
assert 2<=n<=32 and P>0 and L>=0
old=[[]for _ in range(P)]
rep=json.loads(open(path.replace('.trace.json','')).read().splitlines()[1])
assert rep['tau']==t and rep['replicate']==0 and int(rep['seed_hex'],16)==seed
rows=x['rows'];assert [z['rank']for z in rows]==list(range(1,N+1))
dead=False;totals=[]
for row in rows:
 r=row['rank'];natural=[len(allowed(w,n,t,r))for w in old]
 ww=[0]*P if dead else natural
 assert len(row['allowed'])==P and len(row['prefix'])==P and row['allowed']==ww
 groups=[sum(ww[j:j+256]) for j in range(0,P,256)];cdf=[];s=0
 for g in groups:s+=g;cdf.append(str(s))
 assert row['cdf']==cdf
 total=sum(ww);totals.append(total)
 if row.get('terminal'):
  assert r==N and len(row['state'])==P and len(row['moment_groups'])==len(groups)
  assert len(row['root_F_groups'])==len(groups) and len(row['root_G_groups'])==len(groups)
  if dead or total==0:
   assert row['moment_groups']==['0']*len(groups)
   assert row['root_F_groups']==['0']*len(groups) and row['root_G_groups']==['0']*len(groups)
   dead=True
   continue
  rawGroups=[];rawFGroups=[];rawGGroups=[]
  for i,w in enumerate(old):
   assert len(w)==r-1 and row['prefix'][i]==w
   roots,seen,bd=rebuild(w,n,t);o=[0]*32;b=[0]*32
   for v in seen:o[v>>5]|=1<<(v&31)
   for v in bd:b[v>>5]|=1<<(v&31)
   assert row['state'][i]==[o,b,roots[0],roots[1] if len(roots)>1 else 0]
   aa=allowed(w,n,t,r);q,_=draw(key(seed,r,i,2),len(aa));v=aa[q]
   r2=v if r==t else (roots[1] if len(roots)>1 else 0)
   if i%256==0: rawGroups.append(0);rawFGroups.append(0);rawGGroups.append(0)
   rawGroups[-1]+=len(aa)*((roots[0]%n-r2%n)**2+(roots[0]//n-r2//n)**2)
   if t<=N:
    xp=2*(roots[0]%n)-(n-1);yp=2*(roots[0]//n)-(n-1)
    xq=2*(r2%n)-(n-1);yq=2*(r2//n)-(n-1)
    rawFGroups[-1]+=len(aa)*(xp*xp+yp*yp+xq*xq+yq*yq)
    rawGGroups[-1]+=len(aa)*(xp*xq+yp*yq)
  assert row['moment_groups']==[str(z) for z in rawGroups]
  assert row['root_F_groups']==[str(z) for z in rawFGroups]
  assert row['root_G_groups']==[str(z) for z in rawGGroups]
  raw=sum(rawGroups)
  rawF=sum(rawFGroups);rawG=sum(rawGGroups)
  if t>N:
   assert rep['U'] is None and rep['U_over_Z'] is None
   assert all(rep[z] is None for z in ('D','F','G','F_over_Z','G_over_Z','V','C','F_minus_2G_minus_4D'))
  else:
   assert rawF-2*rawG==4*raw
   assert close(rep['U_over_Z'],raw/(total*N))
   assert close(rep['F_over_Z'],rawF/total) and close(rep['G_over_Z'],rawG/total)
   assert close(rep['V'],rawF/(8*total)) and close(rep['C'],rawG/(4*total))
  continue
 if total==0:
  assert row['parents']==[4294967295]*P
  dead=True
  continue
 assert len(row['parents'])==P and len(row['state'])==P
 for i in range(P):
  pa=row['parents'][i];after=row['prefix'][i];bits=row['state'][i]
  assert len(after)==r and len(bits)==4 and len(bits[0])==32 and len(bits[1])==32
  u0,_=draw(key(seed,r,i,1),total);u=(i*total+u0)//P;c=0
  for j,z in enumerate(ww):
   c+=z
   if u<c:assert pa==j;break
  aa=allowed(old[pa],n,t,r);q,_=draw(key(seed,r,i,2),len(aa));w=old[pa]+[aa[q]]
  for st in range(L):
   z=key(seed,r,i,st+3);ii,z=draw(z,r);v,z=draw(z,N);y=w[:]
   if v in y:j=y.index(v);y[ii],y[j]=y[j],y[ii]
   else:y[ii]=v
   if rebuild(y,n,t):w=y
  assert after==w
  roots,seen,bd=rebuild(w,n,t);o=[0]*32;b=[0]*32
  for v in seen:o[v>>5]|=1<<(v&31)
  for v in bd:b[v>>5]|=1<<(v&31)
  assert bits==[o,b,roots[0],roots[1] if len(roots)>1 else 0]
  # Equality with the independently rebuilt 32-word state is the high-word
  # audit.  The n>16 extension must not fabricate zero high words.
 old=row['prefix']
if 0 in totals:
 assert all(z==0 for z in totals[totals.index(0):])
 assert rep['dead'] is True and float(rep['Z'])==0.0
 assert rep['logZ'] is None and rep['U'] is None and rep['U_over_Z'] is None
 assert all(rep[z] is None for z in ('D','F','G','F_over_Z','G_over_Z','V','C','F_minus_2G_minus_4D'))
else:
 ell=sum(math.log(z/(P*(N-r+1))) for r,z in enumerate(totals,1))
 assert rep['dead'] is False and close(rep['logZ'],ell) and close(rep['Z'],math.exp(ell))
 if t<=N:
  raw=sum(map(int,rows[-1]['moment_groups']))
  assert close(rep['U'],math.exp(ell)*raw/(totals[-1]*N))
  if t<=N:
   rawF=sum(map(int,rows[-1]['root_F_groups']));rawG=sum(map(int,rows[-1]['root_G_groups']))
   d=math.exp(ell)*raw/totals[-1];f=math.exp(ell)*rawF/totals[-1];g=math.exp(ell)*rawG/totals[-1]
   assert close(rep['D'],d) and close(rep['F'],f) and close(rep['G'],g)
   assert close(rep['F']-2*rep['G'],4*rep['D'])
   assert abs(float(rep['F_minus_2G_minus_4D']))<=5e-12*max(abs(float(rep['F'])),abs(2*float(rep['G'])),abs(4*float(rep['D'])))
print(json.dumps({'passed':True,'n':n,'tau':t,'L':L,'dead':dead,'logZ_checked':True}))
