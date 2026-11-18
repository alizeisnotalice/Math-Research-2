"""python3 -B check_sharp_loss.py. Exact small spectra, basis certificates and limits."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
import platform


def walsh(row):
    out=list(row)
    w=1
    while w<len(out):
        for s in range(0,len(out),2*w):
            for j in range(s,s+w):
                a,b=out[j],out[j+w]
                out[j],out[j+w]=a+b,a-b
        w*=2
    return out


def rank(xs):
    pivots={}
    for x in xs:
        while x and x.bit_length() in pivots:
            x^=pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()]=x
    return len(pivots)


def powers(vs):
    items=[0]
    for v in vs:
        items.extend([v^b for b in items])
    return items


class Bits:
    def __init__(self):
        self.position=2
    def axes(self,d):
        start=self.position
        self.position+=d
        return [1<<(start+j) for j in range(d)]


def tree_support(leaves,local_dims,bits,compressed):
    nodes=bits.axes(leaves-1)
    rows=set()
    for j,dim in enumerate(local_dims):
        guard=bits.axes(1)[0]
        labels=powers(bits.axes(dim))
        node=leaves+j
        path=[]
        while node>1:
            node//=2
            path.append(nodes[node-1])
        address=[0]+path if compressed else powers(path)
        rows.update(guard|w|v for w in address for v in labels)
    return rows


def exact_parameters(M,seed_n):
    j=M.bit_length()-1
    T=M*M
    m=T.bit_length()-1
    g=M*M
    f=2*(m+j)
    t=(2*T+2)//5
    k=(M+3)//7
    s0=1<<seed_n
    R0=seed_n+1
    E0=R0+1
    K=M*(4*M+4*k)
    RX=4*M+k-1
    L=T//4
    sA=(L*(3*t-2*L))*K*s0
    RA=t*E0+RX-2
    q=T*T*s0
    formula={}
    groups=[(0,1,2,3),*((i,) for i in range(4)),*combinations(range(4),2)]
    for idx in groups:
        if len(idx)==4:
            s,D=4*sA+(24*g+24)*q, RA+3*T*E0+R0+f+4
        elif idx==(0,):
            s,D=sA,RA+1
        elif len(idx)==1:
            s,D=(8*g+2)*q,T*E0+R0+f+3
        elif idx[0]==0:
            s,D=2*sA+2*(8*g+2)*q, RA+T*E0+R0+f+4
        else:
            s,D=(20*g+8)*q,2*T*E0+R0+f+4
        formula[idx]=(s,D)
    return T,t,k,f,formula


def tables(M,seed_n,compressed):
    T,t,k,f,formula=exact_parameters(M,seed_n)
    bits=Bits()
    coarse=bits.axes(2)
    fine=bits.axes(f)
    shared_seed=bits.axes(seed_n)
    X=tree_support(M,[3]*k+[2]*(M-k),bits,compressed)
    gs=tree_support(t,[seed_n]*t,bits,compressed)
    if compressed:
        u,v=min(X),min(gs)
        A={u^b for b in gs}|{v^b for b in X}
    else:
        A={u^v for u in X for v in gs}
    parts=[({0},A)]
    for i in (1,2,3):
        tree=tree_support(T,[seed_n]*T,bits,compressed)
        # For incident edge spots in F2^2, the Fourier support of their
        # complementary two-point line is the annihilator of their difference.
        edge_values={1:(1,2),2:(1,3),3:(2,3)}[i]
        direction=edge_values[0]^edge_values[1]
        window=[u for l,u in enumerate(powers(coarse)) if (l&direction).bit_count()%2==0]
        parts.append(({i},{u^v for u in tree for v in window}))
    for edge in ((1,2),(1,3),(2,3)):
        guard=bits.axes(1)[0]
        fine_cols=[0]+fine if compressed else powers(fine)
        common={guard|v|c|y for v in powers(shared_seed) for c in powers(coarse) for y in fine_cols}
        parts.append((set(edge),common))
    counts={}
    covered=set()
    for _,rows in parts:
        assert not rows.intersection(covered)
        covered.update(rows)
    for idx,(size,D) in formula.items():
        support=set()
        only_A=set()
        for active,rows in parts:
            coefficients=walsh([int(i in active) for i in idx])
            tags={v|h for h,a in enumerate(coefficients) if a for v in rows}
            support.update(tags)
            if active=={0}:
                only_A=tags
        if not compressed:
            assert len(support)==size
        assert rank(2*v+1 for v in support)==D-1,(idx,D)
        if len(idx)==4:
            assert rank(only_A)+rank(support-only_A)-rank(support)==2
        counts[''.join(map(str,idx))]={"s":size,"D":D,"sampled_columns":len(support)}
    return {"M":M,"seed_full_dimension":seed_n,"mode":"sampled actual spanning columns" if compressed else "full actual supports",
            "records":counts}


def local_bound():
    stats=Counter()
    worst=Fraction(0)
    for triple in product(range(-4,5),repeat=3):
        if not any(triple):
            continue
        parent=sum(v!=0 for v in walsh([0,*triple]))
        pair=sum(sum(v!=0 for v in walsh(p)) for p in combinations(triple,2))
        assert 2*pair<=5*parent
        worst=max(worst,Fraction(pair,parent))
        stats[sum(v!=0 for v in triple)]+=1
    assert worst==Fraction(5,2)
    return {"nonzero_inputs_by_occupancy":dict(stats),"max_ratio":str(worst)}


def family_limits():
    cases=[]
    for j in (2,3,4,6,10,16,24):
        M=1<<j
        _,_,_,_,f=exact_parameters(M,0)
        s,D=f[(0,1,2,3)]
        parent=Fraction(D*D,s)
        best=max(Fraction(d*d,size) for key,(size,d) in f.items() if len(key)<4)
        ratio=parent/best
        cases.append({"M":M,"parent_to_best_squared":str(ratio),"decimal":float(ratio)})
    assert abs(ratio-Fraction(17,8))<Fraction(1,10000)
    return cases


if __name__=='__main__':
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "local_bound":local_bound(),"spectra":[tables(2,0,False),tables(4,0,False),tables(2,3,False),tables(8,0,True)],
                      "ratios":family_limits()},indent=2))
