"""Check the weighted peeling construction on deterministic finite affine systems."""
from fractions import Fraction as Q
from itertools import combinations
import random,json,platform

def span(vs):
    S={0}
    for v in vs:
        if v not in S:S|={x^v for x in S}
    return S

def d(S):return len(S).bit_length()-1

def multiplicity(fam,N):
    return max((sum(all(not ((x^b)&v).bit_count()%2 for v in W) for _,b,W in fam)
                for x in range(N)),default=0)

def verify(fam,N):
    orig=sum((list(W) for _,_,W in fam),[])
    full=d(span(orig)); L=multiplicity(fam,N)
    old={i:d(W) for i,_,W in fam};increments=dict.fromkeys(old,0)
    active=list(fam);total_rank=full
    while total_rank>0:
        J=list(active)
        for t in list(J):
            if d(span(v for q in J if q!=t for v in q[2]))==total_rank:J.remove(t)
        new=[]
        for i,b,W in J:
            U=W&span(v for j,_,Z in J if j!=i for v in Z)
            p=d(W)-d(U);assert p>=1
            increments[i]+=p;new.append((i,b,U))
        assert multiplicity(new,N)<=multiplicity(active,N)<=L
        active=new;total_rank=d(span(v for _,_,W in active for v in W))
    assert sum(increments.values())==full
    cost=sum((Q((1+x)**2,1<<old[i]) for i,x in increments.items()),Q())
    assert cost<=7*L
    return cost/L if L else Q()

def main():
    N=8
    spaces={tuple(sorted(span(s))) for m in range(4) for s in combinations(range(1,N),m)}
    flats=[]
    for W in sorted(spaces):
        W=set(W);done=set()
        for b in range(N):
            pattern=tuple((v&b).bit_count()%2 for v in sorted(W))
            if pattern not in done:
                done.add(pattern);flats.append((len(flats),b,W))
    rng=random.Random(6907)
    runs=[]
    for _ in range(2000):
        fam=rng.sample(flats,rng.randrange(1,17))
        runs.append(verify(fam,N))
    print(json.dumps(dict(python=platform.python_version(),seed=6907,ambient=3,affine_flats=len(flats),
      sampled_families=len(runs),max_witnessed_cost_over_overlap=str(max(runs)),analytic_cap=7),indent=2))
if __name__=="__main__":main()
