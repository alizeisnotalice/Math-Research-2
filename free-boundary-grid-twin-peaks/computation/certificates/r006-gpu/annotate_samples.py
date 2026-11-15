"""Independent CPU audit/annotation of saved Metal-search K=2 labelings.

The metrics describe samples only; no frequency is interpreted probabilistically.
"""
import json
import sys
from pathlib import Path


def nbrs(v,n):
    x,y=v%n,v//n
    if x: yield v-1
    if x+1<n: yield v+1
    if y: yield v-n
    if y+1<n: yield v+n


def peaks(a,n):
    return [v for v,z in enumerate(a) if all(a[w]<z for w in nbrs(v,n))]


def first_merge_label(a,n):
    N=len(a); p=list(range(N)); active=[False]*N
    def find(x):
        while p[x]!=x:
            p[x]=p[p[x]];x=p[x]
        return x
    for v in sorted(range(N),key=a.__getitem__,reverse=True):
        active[v]=True; roots={find(w) for w in nbrs(v,n) if active[w]}
        if len(roots)>=2: return a[v]
        for q in roots: p[q]=v
    return None


def basin_interface(a,n,ps):
    # Follow the deterministic highest-neighbour map to a peak.
    root={}
    for v in range(len(a)):
        path=[];u=v
        while u not in root:
            path.append(u); higher=[w for w in nbrs(u,n) if a[w]>a[u]]
            root[u]=u if not higher else -1
            if root[u]==u: break
            u=max(higher,key=a.__getitem__)
        r=root[u]
        for q in path: root[q]=r
    boundary=sum(root[v]!=root[w] for v in range(len(a)) for w in nbrs(v,n) if v<w)
    sizes=[sum(root[v]==p for v in range(len(a))) for p in ps]
    return boundary,sizes


def beta(a,n,ps):
    low=min(ps,key=a.__getitem__); t=a[low]; M=b=0
    for x in nbrs(low,n):
        M=max(M,a[x])
        if all(z==low or a[z]<a[x] for z in nbrs(x,n)): b=max(b,a[x])
    assert M>=b and M<t
    return (M-b)/(t-1),M,b,t


inp=Path(sys.argv[1]); out=Path(sys.argv[2]); checked=0
with inp.open() as f, out.open("w") as g:
    for line in f:
        z=json.loads(line)
        if z.get("summary") or "labels" not in z: continue
        a=z["labels"];n=z["n"];assert sorted(a)==list(range(1,n*n+1));ps=peaks(a,n);assert len(ps)==2
        d=abs(ps[0]%n-ps[1]%n)+abs(ps[0]//n-ps[1]//n);assert d==z["distance"]
        b,M,bb,t=beta(a,n,ps); edge, sizes=basin_interface(a,n,ps)
        z.update({"cpu_verified":True,"peaks":ps,"first_merge_label":first_merge_label(a,n),
                  "beta":b,"beta_M":M,"beta_b":bb,"low_peak_label":t,
                  "gradient_basin_interface_edges":edge,"gradient_basin_sizes":sizes})
        g.write(json.dumps(z,separators=(",",":"))+"\n"); checked+=1
print(json.dumps({"input":str(inp),"output":str(out),"cpu_verified_samples":checked}))
