"""python3 -B check_joint_cut.py. Integer and Fraction arithmetic only."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, isqrt
import platform
import json


def parity(x):
    return x.bit_count() & 1


def greedy(labels):
    pivots = {}
    selected = []
    for v in labels:
        a = v
        while a and a.bit_length() in pivots:
            a ^= pivots[a.bit_length()]
        if a:
            pivots[a.bit_length()] = a
            selected.append(v)
    return selected


def rank(labels):
    return len(greedy(labels))


def aff(labels):
    return rank(2*v+1 for v in labels)


def walsh(values):
    out = list(values)
    w = 1
    while w < len(out):
        for base in range(0,len(out),2*w):
            for j in range(base,base+w):
                x,y=out[j],out[j+w]
                out[j],out[j+w]=x+y,x-y
        w *= 2
    return out


def spectrum(t,n):
    return tuple(j for j,v in enumerate(walsh(1-2*((t>>x)&1) for x in range(2**n))) if v)


def sharp_width(l,M):
    if l == 1:
        return Q(M)
    k=l//2
    zero=comb(l,k) if l%2==0 else 2*comb(2*k,k-1)
    return M*(1-Q(zero,2**l))


@lru_cache(None)
def geometry(supports):
    M=len(supports)
    labels=sorted(set().union(*supports))
    freq={v:tuple(i for i,S in enumerate(supports) if v in S) for v in labels}
    w={v:Q(1,len(a)) for v,a in freq.items()}
    widths=sum(sharp_width(len(a),M) if aff(a)==len(a) else Q(M,len(a)) for a in freq.values())
    W=sum(w.values())
    means=[sum(w[v] for v in S)/len(S) for S in supports]
    ranks=[aff(S) for S in supports]
    costs=[sum(w[x//2] for x in greedy(2*v+1 for v in sorted(S,key=w.__getitem__))) for S in supports]
    r=aff(labels)
    assert sum(costs)>=r
    return (r,sum(ranks)+M,sum(map(len,supports)),len(labels),M*W,widths,
            sum(a*(R+1) for a,R in zip(means,ranks)),
            sum(a*(R+1)-c for a,R,c in zip(means,ranks,costs)))


def subspace(basis):
    items=[0]
    for v in basis:
        items+= [x^v for x in items]
    return items


def two_dimensional_kernels(n):
    result={}
    for v,z in combinations(range(1,2**n),2):
        result.setdefault(tuple(sorted((v,z,v^z))),(v,z))
    return result


def physical_permutation(n,v,z):
    def answer(x):
        return parity(x&v)+2*parity(x&z)
    private=subspace(greedy(x for x in range(2**n) if answer(x)==0))
    inv=subspace([next(x for x in range(2**n) if answer(x)==j) for j in (1,2)])
    order=[u^b for u in inv for b in private]
    return [[sum(((val>>(x-off))&1)<<i for i,x in enumerate(order) if off<=x<off+8)
             for val in range(1<<min(8,2**n-off))] for off in range(0,2**n,8)]


def exhaustive(n):
    small=2**(n-2)
    N=2**n
    tables=2**N
    child=[spectrum(t,n-2) for t in range(2**small)]
    geometry_table={}
    packages=[]
    for t in range(tables):
        parent=spectrum(t,n)
        R=aff(parent)
        S=tuple(child[(t>>(i*small))%(2**small)] for i in range(4))
        if S not in geometry_table:
            r,U,L,size,unc,sharp,weighted,local=geometry(S)
            roots=sum(isqrt(len(A)) for A in S)
            assert all(isqrt(len(A))**2==len(A) for A in S)
            geometry_table[S]=(r,U,L,size,unc,sharp,weighted,local,roots)
        r,U,L,size,unc,sharp,weighted,local,roots=geometry_table[S]
        assert len(parent)>=sharp>=unc
        numeric = (R+1)**2*roots**2 <= U**2*unc
        improved = (R+1)**2*roots**2 <= U**2*sharp
        aggregate=(U-R-1)*L >= (L-size)*U
        weights=R+1<=weighted
        basis=R-r+1<=local
        test=(aggregate,numeric,improved,weights,basis)
        if any(test):
            assert any((R+1)**2*len(A)<=(aff(A)+1)**2*len(parent) for A in S)
        packages.append((R,len(parent),R>r,test))
    stats=Counter()
    successes=set()
    for v,z in two_dimensional_kernels(n).values():
        picks=physical_permutation(n,v,z)
        seen=bytearray(tables)
        for table in range(tables):
            tid=picks[0][table&255]
            if len(picks)>1:
                tid |= picks[1][table>>8]
            assert not seen[tid]
            seen[tid]=1
            R,s,informative,tests=packages[tid]
            assert packages[table][:2]==(R,s)
            stats['table_subspace_pairs']+=1
            if informative:
                stats['positive_rank_loss']+=1
                for key,ok in zip(('overlap','width','sharp_width','weighted_rank','local_greedy'),tests):
                    if ok:
                        stats[key]+=1
                if any(tests):
                    successes.add(table)
        assert all(seen)
    stats['different_tables_with_certified_rank_reduction']=len(successes)
    return {"n":n,"counts":dict(stats)}


def signed_tree(B):
    support={}
    for i,(a,c) in enumerate(((0,1),(2,4),(3,7),(5,6))):
        normals=[u for u in range(8) if not parity(u&(a^c))]
        start=3+i*(2*B-1)
        for leaf in range(B,2*B):
            path={}
            node=leaf
            while node>1:
                path[start+node//2-1]=node&1
                node//=2
            tag=1<<(start+leaf-1)
            for mask in subspace(1<<p for p in path):
                weight=sum((mask>>p & 1)*v for p,v in path.items())
                for u in normals:
                    support[tag|mask|u]=(-1)**(weight+parity(a&u))
    return support


def cut_sparse(signed,dir_basis):
    K=subspace(dir_basis)
    # The joining basis is converted to coefficient masks, so both arbitrary
    # single directions and the two standard address directions are permitted.
    projected={}
    for v,a in signed.items():
        tag=min(v^k for k in K)
        index=K.index(v^tag)
        values=projected.setdefault(tag,[0]*len(K))
        for b in range(len(K)):
            values[b]+=a*(-1)**parity(index&b)
    return tuple(tuple(v for v,values in projected.items() if values[b]) for b in range(len(K)))


def test_tree(B):
    coeff=signed_tree(B)
    S=tuple(coeff)
    assert len(S)==16*B*B and aff(S)==8*B-1
    if B<=4:
        sums=Counter()
        for (u,x),(v,y) in combinations(coeff.items(),2):
            sums[u^v]+=x*y
        assert not any(sums.values())
    slices=cut_sparse(coeff,(1,4))
    assert all(len(A)==4*B*B and aff(A)==4*B-1 for A in slices)
    counts=Counter(v for A in slices for v in A)
    assert set(counts.values())=={2}
    r=aff(counts)
    assert r==8*B-3
    result={"B":B,"s":len(S),"union_rank":r,"multiplicity":2,"omega":2,"theta":1}
    if B==1:
        g=greedy(v^S[0] for v in S)
        cuts=subspace(g)[1:]
        best=Q(0)
        for q in cuts:
            for A in cut_sparse(coeff,(q,)):
                best=max(best,Q((aff(A)+1)**2,len(A)))
        assert best<4
        result.update({"single_directions":len(cuts),"maximum_single_child_P_squared":str(best)})
    return result


def sharp_cases():
    for l in range(1,13):
        weights=[1]*l
        if l%2 and l>1:
            weights[0]=2
        signs=list(product((-1,1),repeat=l))
        zero=sum(sum(a*x for a,x in zip(weights,row))==0 for row in signs)
        assert len(signs)-zero==sharp_width(l,2**l)
    return 12


if __name__=="__main__":
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "two_query_exhaustion":[exhaustive(n) for n in (3,4)],
                      "four_line_trees":[test_tree(B) for B in (1,2,4,8)],
                      "sharp_independent_activity_sizes":sharp_cases()},indent=2))
