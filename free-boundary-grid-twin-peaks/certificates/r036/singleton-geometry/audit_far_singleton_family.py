#!/usr/bin/env python3
"""Strict checks for the explicit far singleton-endpoint family (n=5,6,8,16)."""
import hashlib, json, platform
from pathlib import Path

HERE=Path(__file__).resolve().parent

def grid(n):
    adj=[set() for _ in range(n*n)]
    for y in range(n):
        for x in range(n):
            u=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xx,yy=x+dx,y+dy
                if 0<=xx<n and 0<=yy<n:adj[u].add(yy*n+xx)
    return adj

def comps(seen,adj):
    todo=set(seen);out=[]
    while todo:
        c={todo.pop()};q=list(c)
        while q:
            u=q.pop();hit=adj[u]&todo;todo-=hit;c|=hit;q.extend(hit)
        out.append(c)
    return out

def strict_scan(word,adj):
    seen=set();births=[];merge=None
    for rank,v in enumerate(word,1):
        if not(adj[v]&seen):births.append(v)
        cs=comps(seen,adj) if seen else []
        hit=[c for c in cs if c&adj[v]]
        if merge is None and len(hit)>=2:
            assert len(cs)==len(hit)==2
            merge=(rank,v,hit)
        seen.add(v)
    return births,merge

def family(n):
    def v(x,y):return y*n+x
    r=v(0,0); low=v(n-1,n-1); x=v(n-2,n-1); y=v(n-1,n-2)
    z=v(n-2,n-2); t=v(n-1,n-3)
    A=set(range(n*n))-{low,x,y,z,t}
    adj=grid(n)
    # Deterministic connected order of the high branch, rooted at r.
    order=[r];seen={r}
    while len(order)<len(A):
        frontier=sorted((set().union(*(adj[u] for u in seen))&A)-seen)
        assert frontier
        u=frontier[0];seen.add(u);order.append(u)
    word=(r,low)+tuple(order[1:])+(x,y,z,t)
    assert len(word)==n*n and len(set(word))==n*n
    return word,adj,(r,low,x,y,z,t),len(A)

def main():
    rows=[]
    for n in (5,6,8,16):
        word,adj,data,a=family(n);r,low,x,y,z,t=data
        births,merge=strict_scan(word,adj)
        rank,mv,cs=merge
        assert births==[r,low]
        assert mv==x and rank==a+2
        high=next(c for c in cs if r in c);small=next(c for c in cs if low in c)
        assert len(high)==a and small=={low}
        q=a+1;theta_num=a; theta_den=q
        E_num=2*(n-1)**2; E_den=n*n
        # kappa=q/2 exactly for a singleton side.
        rows.append({"n":n,"births":births,"first_merge_rank":rank,
                     "first_merge_vertex":mv,"branch_sizes":[len(high),len(small)],
                     "q":q,"theta":[theta_num,theta_den],
                     "min_theta":[1,q],"kappa":[q,2],
                     "E":[E_num,E_den]})
    out={"kind":"strict-explicit-far-singleton-family-audit","python":platform.python_version(),
         "rows":rows,"all_checks_passed":True,
         "scope":"checks one explicit top-two family only; no counting or asymptotic mass claim"}
    text=json.dumps(out,indent=2,sort_keys=True)+"\n"
    HERE.joinpath("audit-far-singleton-family.json").write_text(text)
    print(text,end="");print("output_sha256="+hashlib.sha256(text.encode()).hexdigest())
if __name__=="__main__":main()
