"""Directed rectangle proof over the 948 diagonal, with bounded subdivision."""
import argparse,json,os,time
from functools import lru_cache
from pathlib import Path
from certify import D,HI,LO,Box,ZERO,ONE,E,gaussian

N=948
LEFT=D('2.291')
RIGHT=D('3.5')
CUT=D('.125')

@lru_cache(maxsize=16)
def time_cells(j,depth):
    k=N-j
    h=D(1)/(2**depth)
    critical=None if not j else Box(N)/Box(j)
    peak=ONE if j in (0,N) else (Box(j)/Box(N)).pow(j)*(Box(k)/Box(N)).pow(k)
    ans=[]
    def weight(t):
        a=(ZERO-Box(t)).exp()
        return (ZERO-Box(j)*Box(t)).exp()*(ONE-a).pow(k)
    lo=CUT
    wl=weight(lo)
    el=Box(lo).exp()
    while lo<D('1.5'):
        hi=lo+h
        wh=weight(hi)
        eh=Box(hi).exp()
        low=min(wl.a,wh.a)
        max_counted=j in (0,N) or eh.b<=critical.a or el.a>=critical.b
        high=max(wl.b,wh.b) if max_counted else peak.b
        ans.append((lo,hi,Box(low,high)*Box(h)/peak))
        lo,wl,el=hi,wh,eh
    return ans


def check(j,a,b,depth):
    p=D(0);n=D(0)
    for lo,hi,w in time_cells(j,depth):
        if lo<D(1):
            g=gaussian(lo,b).a # b²>1, Gaussian increases with variance on this slab.
            p=LO.add(p,LO.multiply(w.a,g))
        if D('.5')<=lo<D(1):
            va,vb=lo-D('.375'),hi-D('.375')
            xx=a-D(2)
            critical=xx*xx
            vmax=min(max(critical,va),vb)
            g=gaussian(vmax,xx).b
            n=HI.add(n,HI.multiply(w.b,HI.multiply(E.b,g)))
        if lo>=D(1):
            va,vb=lo-D('.875'),hi-D('.875')
            xx=b-D(2)
            g=min(gaussian(va,xx).a,gaussian(vb,xx).a)
            p=LO.add(p,LO.multiply(w.a,LO.multiply(E.a,g)))
    return LO.subtract(1,HI.divide(n,p))


def verify(j,deadline):
    agenda=[(LEFT,RIGHT,0)]; passes=[]; incomplete=[];calls=0
    while agenda:
        a,b,depth=agenda.pop()
        if time.monotonic()>deadline:
            incomplete.extend([(a,b,depth),*agenda]);break
        if abs(j-551)<55:
            timedepth=min(7,5+depth//3)
        else:
            timedepth=min(6,4+depth//3)
        bnd=check(j,a,b,timedepth);calls+=1
        if bnd>0:
            passes.append([str(a),str(b),str(bnd),timedepth]);continue
        if depth==6:
            incomplete.append((a,b,depth));continue
        mid=(a+b)/2
        agenda.extend([(mid,b,depth+1),(a,mid,depth+1)])
    return dict(j=j,k=N-j,calls=calls,passed=passes,unresolved=[[str(a),str(b),d] for a,b,d in incomplete],full=not incomplete)


def main():
    p=argparse.ArgumentParser();p.add_argument('index',type=int);p.add_argument('out',type=Path);p.add_argument('--limit',type=float,default=105)
    a=p.parse_args();start=time.time();cpu=time.process_time(); deadline=time.monotonic()+a.limit
    pairs=list(range(a.index,N+1,15))
    pairs.sort(key=lambda j:-abs(j-551))
    report=[]
    for j in pairs:
        if time.monotonic()>deadline:break
        report.append(verify(j,deadline))
    a.out.write_text(json.dumps(dict(shard=a.index,pid=os.getpid(),start=start,end=time.time(),cpu=time.process_time()-cpu,all_j=pairs,results=report),indent=2)+'\n')

if __name__=='__main__':main()
