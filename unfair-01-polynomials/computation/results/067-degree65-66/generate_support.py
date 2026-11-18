#!/usr/bin/env python3
"""Exact necessary support envelopes for total degrees 65 and 66."""

import argparse
from collections import Counter
from itertools import combinations,product
import json
from pathlib import Path

LEFT=(0,1,3,4)
RIGHT=(0,1,2,5)


def heads(total):
    for a in range(31,total-30):
        b=total-a
        if not abs(a-b)>=2:continue
        if a>b and a-b in LEFT:continue
        if b>a and b-a in RIGHT:continue
        yield a,b


def mask(S):
    return sum(1<<i for i in set(S))


def bits(m):
    while m:
        one=m&-m
        yield one.bit_length()-1
        m-=one


def allowed_closure(a,b,branches):
    total=a+b
    pa,pb=set(LEFT)|{a},set(RIGHT)|{b}
    sa,sb=set(range(8,a+1))|pa,set(range(8,b+1))|pb
    for k,side in enumerate(branches,8):
        if side in "AX":pa.add(k)
        if side in "BX":pb.add(k)
        if side=="A":sb.discard(k)
        if side=="B":sa.discard(k)
    ones=(set(range(31))-{7})|{total}
    zeros={7}
    saturated={a:(a,0),b:(0,b)}
    while True:
        before=(len(pa),len(pb),len(sa),len(sb),len(ones),len(zeros))
        for k,(i0,j0) in saturated.items():
            sa.difference_update(k-j for j in pb if (k-j,j)!=(i0,j0))
            sb.difference_update(k-i for i in pa if (i,k-i)!=(i0,j0))
        for k in zeros:
            sa.difference_update(k-j for j in pb)
            sb.difference_update(k-i for i in pa)
        ones.update(i+j for i in pa for j in pb)
        can_reach={i+j for i in sa for j in sb}
        zeros.update(set(range(total+1))-can_reach)
        if not pa<=sa or not pb<=sb or (ones&zeros):
            return None
        for k in ones:
            pairs=[(i,k-i) for i in sa if k-i in sb]
            if len(pairs)==1:
                i,j=pairs[0]
                pa.add(i);pb.add(j)
        after=(len(pa),len(pb),len(sa),len(sb),len(ones),len(zeros))
        if before==after:
            free=set(range(total+1))-ones-zeros
            return dict(A_positive=sorted(pa),B_positive=sorted(pb),A_allowed=sorted(sa),B_allowed=sorted(sb),
                        A_zero=sorted(set(range(a+1))-sa),B_zero=sorted(set(range(b+1))-sb),
                        one=sorted(ones),zero=sorted(zeros),free=sorted(free),count=1<<len(free))


def direct_envelope(a,b,shape):
    """Independently intersect pairwise endpoint constraints as bit vectors."""
    pa,pb=mask((*LEFT,a)),mask((*RIGHT,b))
    sa=(1<<(a+1))-1
    sb=(1<<(b+1))-1
    sa=sa^mask((2,5,6,7))
    sb=sb^mask((3,4,6,7))
    for k,side in enumerate(shape,8):
        if side=="A":sb&=~(1<<k)
        if side=="B":sa&=~(1<<k)
        if side in "AX":pa|=1<<k
        if side in "BX":pb|=1<<k
    for i,j in product(range(a+1),range(b+1)):
        status=(i+j==a and (i,j)!=(a,0)) or (i+j==b and (i,j)!=(0,b)) or i+j==7
        if not status:continue
        if pa>>i&1:sb&=~(1<<j)
        if pb>>j&1:sa&=~(1<<i)
    allmask=(1<<(a+b+1))-1
    reach=0
    for i in bits(sa):reach|=sb<<i
    reached=0
    for i in bits(pa):reached|=pb<<i
    one=reached|(((1<<31)-1)^(1<<7))|(1<<(a+b))
    zero=(allmask^reach)|(1<<7)
    if pa&~sa or pb&~sb or one&zero:return None
    return tuple(map(tuple,map(bits,(pa,pb,sa,sb,one,zero))))


def minimize(cubes):
    ordered=sorted(set(cubes),key=lambda e:((e[0]|e[1]).bit_count(),e[0],e[1]))
    result=[]
    for x,y in ordered:
        if any(a&x==a and b&y==b for a,b in result):continue
        result.append((x,y))
    return result


def subtract(cube,other):
    """Disjoint partition of cube minus its intersection with other."""
    one,zero=cube
    o,z=other
    if one&z or zero&o:
        return [cube]
    missing=(o|z)&~(one|zero)
    candidates=[]
    for i in bits(missing):
        b=1<<i
        if o&b:
            candidates.append((one,zero|b))
            one|=b
        else:
            candidates.append((one|b,zero))
            zero|=b
    return candidates


def union_cubes(cubes):
    family=minimize(cubes)
    disjoint=[]
    for cube in family:
        pieces=[cube]
        for previous in disjoint:
            pieces=[x for item in pieces for x in subtract(item,previous)]
            if not pieces:break
        disjoint.extend(pieces)
    for (a,b),(u,v) in combinations(disjoint,2):
        assert a&v or b&u
    return disjoint,len(family)


def shannon_count(cubes,n):
    """Independent union cardinality by recursive Boolean restriction."""
    cubes=minimize(cubes)
    if not cubes:return 0
    if cubes[0]==(0,0):return 1<<n
    flags=0
    for o,z in cubes:flags|=o|z
    flag=flags&-flags
    # Cubes constrained to the other value vanish in each child.
    l=[(o,z&~flag) for o,z in cubes if not o&flag]
    r=[(o&~flag,z) for o,z in cubes if not z&flag]
    return shannon_count(l,n-1)+shannon_count(r,n-1)


def one_degree(n):
    per_direction=[]
    all_cubes=[]
    fields=("A_positive","B_positive","A_allowed","B_allowed","one","zero")
    for a,b in heads(n):
        base=allowed_closure(a,b,"???")
        assert base is not None
        rows=[]
        for shape in map("".join,product("ABX",repeat=3)):
            record=allowed_closure(a,b,shape)
            independent=direct_envelope(a,b,shape)
            if record is None:
                assert independent is None
                continue
            assert tuple(tuple(record[k]) for k in fields)==independent
            rows.append(dict(shape=shape,**record))
            all_cubes.append((mask(record["one"]),mask(record["zero"])))
        per_direction.append(dict(degrees=[a,b],base=base,branches=rows,
                                  width_histogram=Counter(len(r["free"]) for r in rows)))
    disjoint,minimum=union_cubes(all_cubes)
    records=[]
    for o,z in disjoint:
        free=[i for i in range(n+1) if not (o|z)>>i&1]
        records.append(dict(one=list(bits(o)),zero=list(bits(z)),free=free,count=1<<len(free)))
    count=sum(row["count"] for row in records)
    assert count==shannon_count(all_cubes,n+1)
    return dict(degree=n,directions=per_direction,maximum_branch_width=max(len(r["free"]) for group in per_direction for r in group["branches"]),
                undominated_cubes=minimum,disjoint_cubes=records,distinct_products_to_check=count)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert not a.out.exists()
    data=dict(status="verified",degrees=[one_degree(n) for n in (65,66)])
    a.out.write_text(json.dumps(data,sort_keys=True,indent=2)+"\n")
    print(json.dumps([dict(degree=r["degree"],distinct_products=r["distinct_products_to_check"],
                           maximum_width=r["maximum_branch_width"],undominated=r["undominated_cubes"],
                           disjoint=len(r["disjoint_cubes"])) for r in data["degrees"]]))


if __name__=="__main__":main()
