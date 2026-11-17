#!/Users/zzh/.venvs/main/bin/python
"""Exact free-cofactor prefix tree for the post-annulus core6 survivors.

Model: A=G P, B=Q Pbar, C=H G Q, where G,Q are arbitrary integer
polynomials of unrestricted degree.  Each level records every live (G,Q)
prefix; a cap is reported as unknown, never as an exclusion.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json
from math import isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
DEPTH=128; CAP=100_000

def floor_pair(p,q,T,D):
    if q==0: return p
    u=2*p+q*T; r=isqrt(q*q*D)
    return (u+r)//2 if q>0 else (u-r-1)//2

def opts(base,T,D):
    p,q=base; f=floor_pair(p,q,T,D)
    return (-f,1-f) if q==0 else (-f,)

def H(b,T,N,l,d,e):
 return (1,2*b,b*b+T,2*b*T+2*l,2*d+b*b*T+2*b*l+N,
 2*e+2*b*d+2*b*N+l*T,2+2*b*e+d*T+b*b*N+b*l*T+l*l,
 2*b+e*T+b*d*T+2*d*l,T+e*b*T+2*e*l+d*d,b*T+2*l+2*d*e,
 2*d+e*e,2*e,1)

def base(coef,z,k):
    return (sum(coef[j][0]*z[k-j] for j in range(1,min(k,6)+1)),
            sum(coef[j][1]*z[k-j] for j in range(1,min(k,6)+1)))

def ck(h,g,q,k):
    return sum(h[j]*sum(g[i]*q[k-j-i] for i in range(k-j+1))
               for j in range(min(k,len(h)-1)+1))

def run(row):
 b,T,N,D,l,d,e=(row[x] for x in ('b','T','N','delta','ell','d','e'))
 # Coefficients in Z+Z*a, with abar=T-a.
 P=((1,0),(b,0),(0,1),(l,b),(d,0),(e,0),(1,0))
 Pb=((1,0),(b,0),(T,-1),(b*T+l,-b),(d,0),(e,0),(1,0))
 h=H(b,T,N,l,d,e); states=[((1,),(1,))]; levels=[]
 for k in range(1,DEPTH+1):
  nxt=[]; rejected={'forced':0,'binary_C':0}
  for g,q in states:
   go=opts(base(P,g,k),T,D); qo=opts(base(Pb,q,k),T,D)
   if k==1: go=tuple(x for x in go if x==-b); qo=tuple(x for x in qo if x==-b)
   if k==3:
    go=tuple(x for x in go if x==-b*g[2]-l); qo=tuple(x for x in qo if x==-b*q[2]-l)
   if not go or not qo: rejected['forced']+=1; continue
   for gi in go:
    for qi in qo:
     gg=g+(gi,); qq=q+(qi,); c=ck(h,gg,qq,k)
     if c in (0,1): nxt.append((gg,qq))
     else: rejected['binary_C']+=1
  levels.append({'k':k,'live':[[list(g),list(q)] for g,q in nxt], 'rejected':rejected})
  if not nxt: return {'row':row,'status':'DEAD','levels':levels}
  if len(nxt)>CAP: return {'row':row,'status':'UNKNOWN_STATE_CAP','levels':levels}
  states=nxt
 return {'row':row,'status':'UNKNOWN_DEPTH','levels':levels}

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--start',type=int,required=True); ap.add_argument('--stop',type=int,required=True); args=ap.parse_args()
 src=HERE/'core6_de_complete_necessary_screen_summary.json'; data=json.loads(src.read_text()); rows=data['survivors'][args.start:args.stop]
 out=HERE/f'prefix_tree_{args.start}_{args.stop}.jsonl.gz'; summary={'dead':0,'unknown_depth':0,'unknown_cap':0,'rows':len(rows),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
 with gzip.open(out,'wt') as f:
  for r in rows:
   rec=run(r); f.write(json.dumps(rec,separators=(',',':'))+'\n'); summary[{'DEAD':'dead','UNKNOWN_DEPTH':'unknown_depth','UNKNOWN_STATE_CAP':'unknown_cap'}[rec['status']]]+=1; f.flush()
 sp=HERE/f'prefix_summary_{args.start}_{args.stop}.json'; summary['tree_sha256']=hashlib.sha256(out.read_bytes()).hexdigest(); summary['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(); sp.write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary))
if __name__=='__main__': main()
