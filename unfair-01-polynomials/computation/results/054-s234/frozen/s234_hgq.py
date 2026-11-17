"""Unrestricted signed integer cofactor prefix trees, via integral H(GQ)."""
import argparse,hashlib,importlib.util,json,math,os,time,sys
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('annulus',ROOT/'s234_annulus.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a);m=a.m
def floor_ring(x,D):return a.mfloor(x,D) if hasattr(a,'mfloor') else exact_floor(x,D)
def exact_floor(x,D):
    R,V,den=m.rv(x,D);B=abs(R)+abs(V)*(D+1)+2;lo,hi=-B,B
    while hi-lo>1:
        mid=(lo+hi)//2
        if m.sign(m.add(x,(-mid,0)),D)>=0:lo=mid
        else:hi=mid
    return lo
def choices(x,D):
    n=-floor_ring(x,D)
    return (n,n+1) if x[1]==0 else (n,)
def state_hash(states):
    h=hashlib.sha256()
    for G,Q,W in sorted(states,key=lambda z:(z[0],z[1])):h.update(json.dumps([G,Q],separators=(',',':')).encode()+b'\n')
    return h.hexdigest()
def tree(row,depth=96,cap=25000,seconds=8):
    D,L,am,an,k0,dm,dn,e,f=row;assert D>1 and math.isqrt(D)**2!=D
    P=[(1,0),(L,0),(am,an),m.add(m.scale(L,(am,an)),(k0,0)),(dm,dn),(e,0),(f,0),(1,0)];B=[m.conj(x,D) for x in P];H=a.norm(P,D);states=[((1,),(1,),(1,))];trace=[];start=time.monotonic();status='unknown_depth'
    for k in range(1,depth+1):
        nxt=[];rejected=Counter(forced_endpoint=0,nonbinary_C=0)
        for G,Q,W in states:
            U=(0,0);V=(0,0)
            for j in range(1,min(7,k)+1):U=m.add(U,m.scale(G[k-j],P[j]));V=m.add(V,m.scale(Q[k-j],B[j]))
            internal=sum(G[i]*Q[k-i] for i in range(1,k));R=internal+sum(H[j]*W[k-j] for j in range(1,min(14,k)+1))
            for g in choices(U,D):
                for q in choices(V,D):
                    if k in (1,3) and (m.add(U,(g,0))!=(0,0) or m.add(V,(q,0))!=(0,0)):rejected['forced_endpoint']+=1;continue
                    if R+g+q not in (0,1):rejected['nonbinary_C']+=1;continue
                    nxt.append((G+(g,),Q+(q,),W+(internal+g+q,)))
        assert len({(g,q) for g,q,w in nxt})==len(nxt);states=nxt;trace.append({'k':k,'live':len(states),'state_sha256':state_hash(states),'rejected':dict(rejected)})
        if not states:status='dead_all_prefixes';break
        if len(states)>cap:status='unknown_state_cap';break
        if time.monotonic()-start>seconds:status='unknown_time_limit';break
    out={'row':row,'status':status,'depth_reached':len(trace),'trace':trace}
    if status!='dead_all_prefixes':out['live_prefixes']=[{'G':list(g),'Q':list(q)} for g,q,w in states]
    return out
