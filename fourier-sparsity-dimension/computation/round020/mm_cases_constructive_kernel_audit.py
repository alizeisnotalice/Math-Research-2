#!/usr/bin/env python3
"""Independent constructive R20 audit for root's 293 standard-MM cases."""
import argparse,hashlib,json,multiprocessing as mp,os,platform,sys,time
from pathlib import Path
N=256
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def span(bs):
 S={0}
 for b in bs:S|={x^b for x in tuple(S)}
 return S
def walsh(mask):
 a=[1-2*((mask>>x)&1) for x in range(N)];h=1
 while h<N:
  for i in range(0,N,2*h):
   for j in range(i,i+h):
    x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
  h*=2
 return a
def bent(mask):return all(abs(x)==16 for x in walsh(mask))
def diff(A):
 pts=[x for x in range(N) if A>>x&1];D={x^y for x in pts for y in pts}
 return sum(1<<d for d in D)
def basis_subspace(dim,avoid=0):
 out=[]
 for x in range(1,16):
  if rank(out+[x])>len(out) and not (dim==3 and x==avoid):out.append(x)
  if len(out)==dim:return out
 raise AssertionError
def fmask(sigma,h):
 return sum(((((x&15)&sigma[x>>4]).bit_count()&1)^((h>>(x>>4))&1))<<x for x in range(N))
def avoid(A,K):
 D=diff(A);return rank(K)==3 and all(not(D>>z&1) for z in span(K)-{0})
def K_four(sigma,h):
 S=[y for y in range(16) if sigma[y]!=y];d=[y^sigma[y] for y in S];e=[(h>>y)&1 for y in S]
 v=next(v for v in range(16) if all(((v&z).bit_count()&1)==1 for z in d))
 W=next(W for W in (span([a,b]) for a in range(1,16) for b in range(a+1,16) if rank([a,b])==2)
        if not ((set(S[i]^S[j] for i in range(4) for j in range(4))- {0})& (W-{0})))
 wb=basis_subspace(2); # replace with a basis of selected W
 wb=[]
 for z in W:
  if rank(wb+[z])>len(wb):wb.append(z)
  if len(wb)==2:break
 return [v]+[z<<4 for z in wb]
def K_two_one(sigma,h):
 S=[y for y in range(16) if sigma[y]!=y];a,b=S;d=a^b
 full=[y for y in range(16) if sigma[y]==y and (h>>y)&1];assert len(full)==1;c=full[0]
 ea,eb=(h>>a)&1,(h>>b)&1
 lam=next(t for t in range(1,16) if ((t&(c^a)).bit_count()&1)==1 and ((t&(c^b)).bit_count()&1)==1)
 W=[w for w in range(16) if ((w&lam).bit_count()&1)==0];assert d in W
 target=1^ea^eb;u=next(x for x in range(16) if ((x&d).bit_count()&1)==target)
 wb=[d]
 for w in W:
  if rank(wb+[w])>len(wb):wb.append(w)
  if len(wb)==3:break
 return [u|(d<<4),wb[1]<<4,wb[2]<<4]
def K_two_full(sigma,h):
 ys=[y for y in range(16) if (h>>y)&1];assert len(ys)==2;delta=ys[0]^ys[1]
 lam=next(t for t in range(1,16) if ((t&delta).bit_count()&1)==1)
 W=[w for w in range(16) if ((w&lam).bit_count()&1)==0];wb=[]
 for w in W:
  if rank(wb+[w])>len(wb):wb.append(w)
  if len(wb)==3:break
 assert delta not in span(wb)
 return [w<<4 for w in wb]
def audit(row):
 sigma=row['sigma'];h=row['h'];f=fmask(list(range(16)),0);g=fmask(sigma,h);A=f^g
 assert A==int(row['support_hex'],16) and A.bit_count()==32 and bent(f) and bent(g)
 c=sum(sigma[y]!=y for y in range(16));t=sum(sigma[y]==y and (h>>y)&1 for y in range(16));assert c+2*t==4
 K=K_four(sigma,h) if c==4 else K_two_one(sigma,h) if c==2 else K_two_full(sigma,h)
 assert avoid(A,K)
 return {'kind':row['kind'],'support_hex':row['support_hex'],'case_type':[c,t],
         'kernel_basis':K,'kernel_rank':3,'direct_AplusA_avoidance':True,'worker_pid':os.getpid()}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);ap.add_argument('--output',type=Path,default=Path(__file__).with_name('mm-cases-constructive-kernel-audit.json'))
 a=ap.parse_args();assert 1<=a.workers<=12;started=time.monotonic();raw=a.input.read_bytes();doc=json.loads(raw);assert len(doc['cases'])==293
 with mp.get_context('fork').Pool(a.workers) as pool:rows=pool.map(audit,doc['cases'])
 out={'schema':'a28-r020-mm-constructive-kernel-audit-v1','scope':'root 293 standard-coordinate MM representatives only','arithmetic':'exact F2 and integer FWHT','input_path':str(a.input),'input_sha256':hashlib.sha256(raw).hexdigest(),'cases':len(rows),'all_verified':True,'worker_pids':sorted({r['worker_pid'] for r in rows}),'workers_requested':a.workers,'rows':rows,'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-started,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['cases','all_verified','worker_pids','elapsed_seconds']},indent=2))
if __name__=='__main__':main()
