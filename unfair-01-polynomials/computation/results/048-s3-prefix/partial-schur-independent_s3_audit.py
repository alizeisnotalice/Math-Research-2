#!/usr/bin/env python3
"""Independent three-TN replay: descending Schur and direct P*Pbar convolution."""
import gzip,hashlib,json,math,time
from collections import Counter
from pathlib import Path
ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round048-parent');PREFIX=ROOT/'s3-prefix6-independent.json';OUT=Path(__file__).parent/'result.json';PTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3))
def bd(p,q,T,D,a,b):
 r=2*a-b*abs(2*p+q*T);return r>0 and b*b*q*q*D<r*r
def ts():
 z=[]
 for T in range(-16,17):
  for N in range(-71,72):
   D=T*T-4*N
   if D>0 and math.isqrt(D)**2!=D and bd(0,1,T,D,203,24):z.append((T,N,D))
 return z
def schur_desc(d):
 d=list(d)
 while len(d)>1:
  g=math.gcd(*d);g=-g if d[0]<0 else g;d=[x//g for x in d]
  if abs(d[-1])>=d[0]:return len(d)-1
  lead,last=d[0],d[-1];n=len(d)-1;d=[lead*d[j]-last*d[n-j] for j in range(n)]
 return None
def ppbar(c,d,e,f,T,N):
 # Direct K=Q(alpha) multiplication, alpha^2=T*alpha-N.
 F=(1,0,c,0,d,e,f,1);P=[(x,1 if i==3 else 0) for i,x in enumerate(F)]
 Q=[(u+v*T,-v) for u,v in P];H=[(0,0)]*15
 for i,(u,v) in enumerate(P):
  for j,(r,s) in enumerate(Q):
   a,b=H[i+j];H[i+j]=(a+u*r-v*s*N,b+u*s+v*r+v*s*T)
 assert all(v==0 for u,v in H);return [u for u,v in H]
def one(tr,allow):
 T,N,D=tr;C=Counter(TN=1);codes=bytearray();sur=[];tests=[]
 for p,q in PTS:
  w=[p**i*q**(7-i) for i in range(8)];tests.append((2*(w[0]+w[7])+T*w[3],2*w[2],2*w[4],2*w[5],2*w[6],D*w[3]**2))
 for c in range(-6,7):
  C['c']+=1;es=[e for e in range(-29,30) if bd(-5*e,5*c,T,D,3107,32)];permit=(T,N,D,c) in allow
  for d in range(-52,53):
   if abs(2*c*c-4*d)*16>=787:continue
   C['d']+=1;C['e']+=len(es);ge=[e for e in es if bd(7*(c*e-1),7*(d-c*c),T,D,49283,128)];C['s7']+=len(ge);fs=[f for f in range(-8,9) if bd(-3*N-2*c**3+6*c*d-6*f,3*T,T,D,12355,64)];C['f']+=len(ge)*len(fs)
   for e in ge:
    for f in fs:
     if not permit:C['prefix_rejected']+=1;codes.append(0);continue
     C['prefix_allowed']+=1;code=None
     for j,(a,bc,dd,ee,ff,V) in enumerate(tests):
      L=a+bc*c+dd*d+ee*e+ff*f
      if L<=0 or L*L<=V:code=1+j;C[f'axis_{j}']+=1;break
     if code is None:
      C['positive_axis']+=1;H=ppbar(c,d,e,f,T,N);z=schur_desc([(v<<i) for i,v in enumerate(H)][::-1])
      if z is not None:code=20+z;C['outer_rejected']+=1
      else:
       z=schur_desc([H[j]<<(14-j) for j in range(15)])
       if z is not None:code=40+z;C['inner_rejected']+=1
       else:code=255;C['annulus_survived']+=1;sur.append([T,N,D,c,d,e,f])
     codes.append(code)
 return dict(C),bytes(codes),sur
def main():
 p=json.loads(PREFIX.read_text());allow={tuple(x[:4]) for x in p['states']};all=ts();assert len(all)==669;out={'indices':{},'prefix_sha256':hashlib.sha256(PREFIX.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};t=time.monotonic()
 for i in (0,334,668,109,148,60):
  c,b,s=one(all[i],allow);out['indices'][str(i)]={'triple':all[i],'counts':c,'survivors':s,'raw_code_sha256':hashlib.sha256(b).hexdigest(),'gzip_code_sha256':hashlib.sha256(gzip.compress(b,mtime=0)).hexdigest(),'code_len':len(b)}
 out['seconds']=time.monotonic()-t;OUT.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
if __name__=='__main__':main()
