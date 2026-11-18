"""python3 -B check_zero_band.py. Exhaustive two-bit branches on separate data blocks."""
import json
import platform
from collections import Counter
from itertools import combinations, product


def walsh(v):
    w=list(v)
    size=1
    while size<len(w):
        for lower in range(0,len(w),2*size):
            for j in range(lower,lower+size):
                a,b=w[j],w[j+size]
                w[j],w[j+size]=a+b,a-b
        size*=2
    return w


def rank(values):
    basis={}
    for v in values:
        while v and v.bit_length() in basis:
            v^=basis[v.bit_length()]
        if v:
            basis[v.bit_length()]=v
    return len(basis)


def aff(values):
    return rank(2*v+1 for v in values)


def support(rows):
    result=set()
    width=len(rows)
    for v in set().union(*rows):
        for j,a in enumerate(walsh(row.get(v,0) for row in rows)):
            if a:
                result.add(v*width+j)
    return result


def branch_tables(shift):
    out=[]
    for truth in range(1,15):
        w=walsh([1-2*((truth>>i)&1) for i in range(4)])
        out.append({i<<shift:a for i,a in enumerate(w) if a})
    return out


def main():
    counts=Counter()
    As,Bs=branch_tables(2),branch_tables(0)
    for group in product(Bs,repeat=3):
        rB=aff(set().union(*group))
        lines=[support([group[i],group[j]]) for i,j in combinations(range(3),2)]
        delta=sum(map(aff,lines))-2*rB
        assert delta>=1
        for A in As:
            root=support([A,*group])
            left={v for v in root if (v>>2) & 12}
            right=root-left-{0}
            if rank(left)+rank(right)-rank(root)!=2:
                continue
            D=aff(root)+1
            let=len(root)
            cap=4*len(A)*5+2*sum(map(len,lines))
            row_cost=2*(aff(A)+1)+sum(aff(S)+1 for S in lines)
            zero=0 in A and any(0 in B for B in group)
            if zero:
                counts['shared_zero']+=1
                assert delta>=2 and 2*D+3<=row_cost
                assert cap<=5*let+22 and let>=2*D-4 and D>=7
                assert 22*D*D*4<=5*let*(12*D+9)
            else:
                counts['disjoint_zero']+=1
                assert 2*D<=row_cost and cap<=5*let
            # It suffices in the proof to use A and the three B-only lines.
            assert any(8*D*D*len(S)<=17*let*(aff(S)+1)**2 for S in [A,*lines])
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "scope":"all four independently selected nonconstant two-bit Boolean branch tables; two independent data spaces",
                      "table_tuples":len(As)*len(Bs)**3,"accepted_interface_cases":dict(counts)},indent=2))


if __name__=='__main__':
    main()
