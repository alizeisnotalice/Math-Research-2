#!/usr/bin/env python3
"""Exact rational prefix feasibility; never a full polynomial certificate."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,time,os

def run(s,t,a,p=2,q=21):
    G={0:Q(1),s:a,t:Q(1),q:Q(1),q+s:a,q+t:Q(1)}
    L=2*q; forbidden={z for j in G for z in (j-p,j+p) if 0<=z<L}
    f=[];h=[]; nodes=0;leaves=0; first=None
    def dfs(i):
        nonlocal nodes,leaves,first
        nodes+=1
        if i==L:
            leaves+=1
            if first is None:first=(f.copy(),h.copy())
            return
        c=sum((g*f[i-j] for j,g in G.items() if 0<j<=i),Q(0))
        if c>1:return
        forced_zero=i in forbidden or (i>=p and h[i-p]==1)
        forced_one=i in G
        choices=(0,) if forced_zero else ((1,) if forced_one or c>0 else (0,1))
        if forced_zero and forced_one:return
        for hi in choices:
            fi=Q(hi)-c
            if fi<0:continue
            f.append(fi);h.append(hi)
            # All F involved in these A/B rows are already known.
            valid=True
            for u,v in ((p,q),(q,p)):
                k=u*i
                for e in range(k,k+u):
                    if e>=u*L:continue
                    coeff=sum((g*f[(e-v*j)//u] for j,g in G.items() if e>=v*j and (e-v*j)%u==0 and (e-v*j)//u<=i),Q(0))
                    if coeff>1:valid=False;break
                if not valid:break
            if valid:dfs(i+1)
            f.pop();h.pop()
    dfs(0)
    return {'p':p,'q':q,'s':s,'t':t,'a':str(a),'prefix_length':L,'nodes':nodes,'leaves':leaves,'first':None if first is None else {'F_prefix':[str(x) for x in first[0]],'H_prefix':first[1]},'scope':'prefix only; missing tail and full F/G factorization'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--s',type=int,required=True);ap.add_argument('--t',type=int,required=True);ap.add_argument('--a',default='1/2');ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data');start=time.monotonic();r=run(args.s,args.t,Q(args.a));r['seconds']=time.monotonic()-start;args.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='first'}))
