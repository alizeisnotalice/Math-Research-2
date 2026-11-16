#!/usr/bin/env python3
import json,sys
MASK=(1<<64)-1; STEP=0x9e3779b97f4a7c15; C1=0xbf58476d1ce4e5b9; C2=0x94d049bb133111eb
def mix(x): x=((x^(x>>30))*C1)&MASK; x=((x^(x>>27))*C2)&MASK; return x^(x>>31)
def draw(x,m):
 t=(-m)&MASK; t%=m
 while 1:
  x=(x+STEP)&MASK; z=mix(x)
  if z>=t:return z%m,x
def key(s,r,c,d):return s^(r<<32)^((c*STEP)&MASK)^((d*0xd1b54a32d192ed03)&MASK)
def rebuild(w,n,tau):
 if not w:return []
 seen=set(); bd=set(); roots=[]
 for q,v in enumerate(w,1):
  if v in seen or not 0<=v<n*n:return None
  birth=v not in bd
  if birth:roots.append(v)
  if birth != (q==1 or q==tau):return None
  seen.add(v);x=v%n;y=v//n
  for z in ([v-1] if x else [])+([v+1] if x+1<n else [])+([v-n] if y else [])+([v+n] if y+1<n else []):
   if z not in seen:bd.add(z)
  bd.discard(v)
 return roots if len(roots)==(1 if len(w)<tau else 2) else None
def allowed(w,n,tau,r):
 roots=rebuild(w,n,tau);assert roots is not None
 seen=set(w);bd=set()
 for v in w:
  x=v%n;y=v//n
  for z in ([v-1] if x else [])+([v+1] if x+1<n else [])+([v-n] if y else [])+([v+n] if y+1<n else []):
   if z not in seen:bd.add(z)
 free=set(range(n*n))-seen
 return sorted(free if r==1 else (free-bd if r==tau else bd))
def main(p):
 d=json.load(open(p));n=d['n'];tau=d['tau'];P=d['P'];seed=int(d['seed_hex'],16);L=d['L'];old=[[] for _ in range(P)];checks=0
 assert d['method'].startswith('stratified')
 assert [x['rank'] for x in d['traces']]==list(range(1,len(d['traces'])+1)) and d['traces'][-1]['rank']<=n*n
 for row in d['traces']:
  r=row['rank']
  assert len(row['states'])==P and len(row['roots'])==P and len(row['occ'])==P and len(row['bd'])==P
  for w,o,b,rr in zip(row['states'],row['occ'],row['bd'],row['roots']):
   roots=rebuild(w,n,tau);assert roots is not None and rr==roots+[0]*(2-len(roots))
   oo=[0]*8;bb=[0]*8
   for v in w:oo[v>>5]|=1<<(v&31)
   seen=set(w)
   for v in w:
    x=v%n;y=v//n
    for z in ([v-1] if x else [])+([v+1] if x+1<n else [])+([v-n] if y else [])+([v+n] if y+1<n else []):
     if z not in seen:bb[z>>5]|=1<<(z&31)
   assert o==oo and b==bb
  if row.get('dead'):
   ww=[len(allowed(w,n,tau,r)) for w in old];assert row['allowed']==ww and row['cdf']==[sum(ww[:i+1]) for i in range(P)] and not any(ww);break
  if row.get('terminal'):
   ww=[len(allowed(w,n,tau,r)) for w in old];assert row['allowed']==ww and row['cdf']==[sum(ww[:i+1]) for i in range(P)]
   for i,w in enumerate(old):
    a=allowed(w,n,tau,r); x,_=draw(key(seed,r,i,2),len(a)); assert row['states'][i]==w+[a[x]];checks+=1
   break
  weights=[len(allowed(w,n,tau,r)) for w in old]; total=sum(weights);assert row['allowed']==weights and row['cdf']==[sum(weights[:i+1]) for i in range(P)]
  for i,(pa,v,after,moves) in enumerate(zip(row['parents'],row['child'],row['states'],row['moves'])):
   raw,_=draw(key(seed,r,i,1),total);u=(i*total+raw)//P; acc=0
   for j,z in enumerate(weights):
    acc+=z
    if u<acc:assert pa==j;break
   aa=allowed(old[pa],n,tau,r);q,_=draw(key(seed,r,i,2),len(aa));assert v==aa[q];w=old[pa]+[v]
   assert len(moves)==L
   for step,rec in enumerate(moves):
    rng=key(seed,r,i,step+3);ii,rng=draw(rng,r);vv,rng=draw(rng,n*n); y=w[:];
    if vv in y:j=y.index(vv);y[ii],y[j]=y[j],y[ii]
    else:y[ii]=vv
    ok=rebuild(y,n,tau) is not None
    if rec != [ii,vv,int(ok)]: raise AssertionError(('move',r,i,step,rec,[ii,vv,int(ok)],w))
    if ok:w=y
   assert after==w and rebuild(after,n,tau) is not None;checks+=1
  old=row['states']
 print(json.dumps({'file':p,'passed':True,'checks':checks,'n':n,'tau':tau,'L':L}))
if __name__=='__main__':main(sys.argv[1])
