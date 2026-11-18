"""Independent integer audit of cavity weights and immediate joins."""
from collections import Counter,defaultdict,deque
from fractions import Fraction
from itertools import permutations
from math import factorial
from pathlib import Path
from random import Random
import argparse
import json
import sys

D=Path(__file__).resolve().parent


def connected_component(mask,start,around):
    seen=start
    while more:=around[seen]&mask&~seen:
        seen|=more
    return seen


def cells(n,grid):
    N=n*n
    adj=grid(n)
    full=(1<<N)-1
    around=[0]*(full+1)
    suffix=[0]*(full+1)
    for s in range(1,full+1):
        b=s&-s
        around[s]=around[s^b]|adj[b.bit_length()-1]
    suffix[full]=1
    for s in range(full-1,0,-1):
        frontier=around[s]&~s
        while frontier:
            b=frontier&-frontier
            suffix[s]+=suffix[s|b]
            frontier^=b
    edges=prefixes=0
    best=defaultdict(Fraction)
    for K in range(1,full+1):
        if connected_component(K,K&-K,around)!=K:
            continue
        prefixes+=1
        U=full&~(K|around[K])
        rest=U
        M=N-K.bit_count()
        dist=[0]*N
        flood=K
        depth=0
        while more:=around[flood]&~flood:
            depth+=1
            flood|=more
            for v in range(N):
                if (more>>v)&1:
                    dist[v]=depth
        while rest:
            C=connected_component(U,rest&-rest,around)
            rest^=C
            b=C.bit_count()
            for q in range(N):
                if not (C>>q)&1:
                    continue
                for u in range(q+1,N):
                    if ((C&adj[q])>>u)&1:
                        weight=suffix[K|(1<<q)|(1<<u)]
                        ratio=Fraction(weight*M*(M-1),suffix[K])
                        ell=min(dist[q],b)
                        path_bound=3*factorial(ell)*(2*(b+1))**(ell+1)
                        assert 1<=ratio<=min(factorial(b+1),path_bound)
                        best[b]=max(best[b],ratio)
                        edges+=1
    return {'n':n,'connected_prefix_sets':prefixes,'internal_good_edges':edges,
            'largest_edge_multiplier_by_cavity_size':{b:str(v) for b,v in sorted(best.items())}}


def all_bad(grid,rho,roots,trace):
    adj=grid(3)
    orig=Counter()
    outputs=defaultdict(Counter)
    for word in permutations(range(9)):
        rs=roots(word,adj)
        if len(rs)==1:
            orig['one_peak']+=1
        if len(rs)!=2:
            continue
        p,q=rs
        seq,j0,K,U=trace(word,q,adj)
        u=seq[j0]
        if u in U:
            continue
        out=seq[:j0+1]+(q,)+seq[j0+1:]
        assert roots(out,adj)==(p,) and rho(u,q,3)==1
        slot=word.index(q)
        assert 1<=slot<=j0
        restored=seq[:slot]+(q,)+seq[slot:]
        assert restored==word
        outputs[j0+1][out]+=1
    report={}
    for j,obs in sorted(outputs.items()):
        assert len(obs)<=orig['one_peak'] and max(obs.values())<=j-1
        report[j]={'bad_words':sum(obs.values()),'one_peak_images':len(obs),
                   'largest_inverse':max(obs.values())}
    return {'n':3,'one_peak':orig['one_peak'],'by_first_contact_rank':report}


def random_box():
    rng=Random(2026111854)
    count=0
    for n in (16,24,40):
        world={(x,y) for x in range(n) for y in range(n)}
        def boundary(S):
            return {v for x,y in S for v in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if v in world-S}
        for _ in range(100):
            p=(rng.randrange(n),rng.randrange(n))
            K={p}
            for _ in range(rng.randrange(n//4)):
                K.add(rng.choice(sorted(boundary(K))))
            assert len(K)<=n//4
            left=max(0,min(x for x,y in K)-1)
            right=min(n-1,max(x for x,y in K)+1)
            bottom=max(0,min(y for x,y in K)-1)
            top=min(n-1,max(y for x,y in K)+1)
            outside={v for v in world if not (left<=v[0]<=right and bottom<=v[1]<=top)}
            root=tuple(0 if a>=n//2 else n-1 for a in p)
            assert root in outside and sum(abs(a-b) for a,b in zip(p,root))>=n-1
            U=world-K-boundary(K)
            assert outside<=U
            found={root}
            todo=deque([root])
            while todo:
                x,y=todo.popleft()
                for v in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
                    if v in outside and v not in found:
                        found.add(v)
                        todo.append(v)
            assert found==outside
            count+=1
    return count


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('base',type=Path)
    args=p.parse_args()
    sys.path.insert(0,str(args.base))
    from audit_distance_layers import grid,rho,roots
    from verify_contact import trace
    print(json.dumps({'cavity_weights':[cells(n,grid) for n in (3,4)],
                      'bad_contact_words':all_bad(grid,rho,roots,trace),
                      'larger_box_checks':random_box()},indent=2))
