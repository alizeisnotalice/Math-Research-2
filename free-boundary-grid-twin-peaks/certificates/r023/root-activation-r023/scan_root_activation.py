#!/usr/bin/env python3
"""Exact completion-count scan of low-root first activation on n=3,4 grids."""
from fractions import Fraction
from decimal import Decimal, getcontext
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
getcontext().prec = 60


def build(n):
    N=n*n; adj=[0]*N
    for v in range(N):
        x,y=v%n,v//n
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            X,Y=x+dx,y+dy
            if 0<=X<n and 0<=Y<n: adj[v]|=1<<(Y*n+X)
    full=(1<<N)-1
    nbr=[0]*(1<<N)
    # Incremental neighbour unions avoid repeated vertex scans in h recurrence.
    for m in range(1,1<<N):
        bit=m & -m; v=bit.bit_length()-1
        nbr[m]=nbr[m^bit]|adj[v]
    h=[0]*(1<<N); h[full]=1
    for size in range(N-1,0,-1):
        for m in range(1,1<<N):
            if m.bit_count()!=size: continue
            b=nbr[m]&~m
            total=0
            while b:
                bit=b&-b; total+=h[m|bit]; b^=bit
            h[m]=total
    return N,adj,full,h


def connected(mask,adj):
    if not mask:return False
    todo=mask&-mask; seen=0
    while todo:
        bit=todo&-todo; todo^=bit; seen|=bit
        v=bit.bit_length()-1
        todo|=adj[v]&mask&~seen
    return seen==mask


def coords(mask,n):
    return [[v%n,v//n] for v in range(n*n) if mask>>v&1]


def describe(record,n,adj,h):
    C,y,num,den=record['C'],record['y'],record['num'],record['den']
    S=C|(1<<y)
    terms=[]
    ny=adj[y]
    while ny:
        bit=ny&-ny; z=bit.bit_length()-1; ny^=bit
        terms.append({"z": [z%n,z//n], "h_S_plus_z": str(h[S|bit]),
                      "also_adjacent_C": bool(adj[z]&C)})
    p=Fraction(num,den); c=C.bit_count()
    return {"C_coordinates":coords(C,n),"y":[y%n,y//n],"C_size":c,
            "h_S":str(den),"activation_numerator":str(num),"p_act":f"{p.numerator}/{p.denominator}",
            "p_times_C_plus_1":f"{(p*(c+1)).numerator}/{(p*(c+1)).denominator}",
            "p_squared_times_C":f"{(p*p*c).numerator}/{(p*p*c).denominator}",
            "neighbour_completion_terms":terms}


def scan(n):
    N,adj,full,h=build(n)
    records=[]
    for C in range(1,1<<N):
        if not connected(C,adj): continue
        c=C.bit_count()
        for y in range(N):
            ybit=1<<y
            if C&ybit or adj[y]&C: continue
            S=C|ybit
            # Every y-neighbour is unrevealed, distinct, and belongs to B(S).
            ny=adj[y]; num=0
            while ny:
                bit=ny&-ny; ny^=bit; num+=h[S|bit]
            den=h[S]
            assert den>0 and 0<num<=den
            records.append({"C":C,"y":y,"num":num,"den":den})
    def p(r): return Fraction(r['num'],r['den'])
    min_p=min(records,key=p)
    min_linear=min(records,key=lambda r:p(r)*(r['C'].bit_count()+1))
    min_sqrt=min(records,key=lambda r:p(r)*p(r)*r['C'].bit_count())
    ps=p(min_sqrt); sqrtvalue=(Decimal(ps.numerator)/Decimal(ps.denominator)*Decimal(min_sqrt['C'].bit_count()).sqrt())
    by_size=[]
    for k in sorted({r['C'].bit_count() for r in records}):
        layer=[r for r in records if r['C'].bit_count()==k]
        layer_p=min(layer,key=p)
        layer_linear=min(layer,key=lambda r:p(r)*(k+1))
        layer_sqrt=min(layer,key=lambda r:p(r)*p(r)*k)
        by_size.append({"C_size":k,"eligible_pairs":len(layer),
                        "minimum_p_act":describe(layer_p,n,adj,h),
                        "minimum_p_times_C_plus_1":describe(layer_linear,n,adj,h),
                        "minimum_p_times_sqrt_C":describe(layer_sqrt,n,adj,h),
                        "minimum_p_times_sqrt_C_decimal":str(
                            (Decimal(p(layer_sqrt).numerator)/Decimal(p(layer_sqrt).denominator)*Decimal(k).sqrt()))})
    out={"n":n,"N":N,"method":"exact integer h(S)=sum_{z in B(S)}h(S+z), all connected C and y nonadjacent to C",
         "eligible_pairs":len(records),
         "minimum_p_act":describe(min_p,n,adj,h),
         "minimum_p_times_C_plus_1":describe(min_linear,n,adj,h),
         "minimum_p_times_sqrt_C":describe(min_sqrt,n,adj,h),
         "minimum_p_times_sqrt_C_decimal":str(sqrtvalue),
         "by_C_size":by_size,
         "checks":{"unique_y_neighbour_sum":True,"all_y_neighbours_unrevealed":True,
                   "all_y_neighbours_are_boundary":True,"all_probabilities_in_0_1":True}}
    return out


def main():
    result={"results":[scan(3),scan(4)]}
    (HERE/'root-activation-scan-n3-n4.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
