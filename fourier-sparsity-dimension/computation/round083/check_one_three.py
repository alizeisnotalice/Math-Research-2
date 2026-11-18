"""python3 -B check_one_three.py. Integer coefficients and rational comparisons."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform


def pivots(vectors):
    basis = {}
    for x in vectors:
        while x and x.bit_length() in basis:
            x ^= basis[x.bit_length()]
        if x:
            basis[x.bit_length()] = x
    return basis


def rank(vectors):
    return len(pivots(vectors))


def walsh(values):
    out = list(values)
    w = 1
    while w < len(out):
        for base in range(0,len(out),2*w):
            for j in range(base,base+w):
                a,b=out[j],out[j+w]
                out[j],out[j+w]=a+b,a-b
        w *= 2
    return out


class Bits:
    def __init__(self):
        self.pos=2  # reserve low bits for the two interface frequencies
    def new(self,n):
        index=self.pos
        self.pos+=n
        return index


def choice_tree(T, bits, seed, seed_n, compressed):
    # Integers in this block are normalized tree coefficients times T*2**seed_n.
    internal=bits.new(T-1)
    out={}
    for leaf in range(T,2*T):
        guard=1<<bits.new(1)
        offset=bits.new(seed_n)
        masks=[(0,1)]
        node=leaf
        while node>1:
            tag=1<<(internal+node//2-1)
            sign=(-1)**(node&1)
            if compressed:
                masks.append((tag,sign))
            else:
                masks += [(q^tag,s*sign) for q,s in masks]
            node//=2
        for b,c in enumerate(seed):
            if c:
                out.update({guard|(b<<offset)|q:c*s for q,s in masks})
    return out


def branch_rows(T, seed_n, compressed):
    m=T.bit_length()-1
    N=T*T
    D0=2**seed_n
    seed=walsh(([-1]+[1]*(D0-1)) if seed_n else [1])
    common_denom=32*N*D0
    pool=Bits()
    v_base=pool.new(2)
    y_base=pool.new(2*m)
    w=1<<pool.new(1)
    common_seed=pool.new(seed_n)
    octet=pool.new(3)
    rows=[{} for _ in range(4)]
    oct_w=walsh([-1]+[1]*7)
    A=choice_tree(T//2,pool,seed,seed_n,compressed)
    for data,c in A.items():
        for q,z in enumerate(oct_w):
            rows[0][data|(q<<octet)]=c*z*(common_denom//(8*(T//2)*D0))
    coarse=[walsh([int(v not in (1,2)) for v in range(4)]),
            walsh([int(v!=1) for v in range(4)]),walsh([int(v!=2) for v in range(4)])]
    for i in range(3):
        G=choice_tree(T,pool,seed,seed_n,compressed)
        for q,z in enumerate(coarse[i]):
            if z:
                rows[i+1].update({data|(q<<v_base):c*z*(common_denom//(4*T*D0))
                                  for data,c in G.items()})
    fine=(0,*(1<<j for j in range(2*m))) if compressed else range(N)
    for y in fine:
        u=N-2 if y==0 else -2
        v=(-1)**(((y%T)&(y//T)).bit_count()%2)
        for q in range(4):
            for b,z in enumerate(seed):
                if z:
                    tag=w|(b<<common_seed)|(y<<y_base)|(q<<v_base)
                    a=u*z*(-1)**(q&1)*(common_denom//(4*N*D0))
                    c=v*z*(-1)**(q>>1)*(common_denom//(4*T*D0))
                    assert len(set(map(abs,(a+c,a,c))))==3 and a+c!=0
                    rows[1][tag],rows[2][tag],rows[3][tag]=a+c,a,c
    if not compressed:
        assert all(sum(x*x for x in row.values())==common_denom**2 for row in rows)
    return rows,seed


def gather(rows, indices):
    result={}
    for tag in set().union(*(rows[i] for i in indices)):
        data=[rows[i].get(tag,0) for i in indices]
        for q,value in enumerate(walsh(data)):
            if value:
                result[tag|q]=value
    return result


def full_cardinality_stream(T,seed_n):
    # Count all distinct data labels without materializing their four outer
    # frequency copies.  The four private guard sectors and common guard
    # sector are disjoint. Within a sector, private coordinates, coarse/fine
    # coordinates and seed labels are separate, so these loops enumerate each
    # actual label exactly once.
    D0=2**seed_n
    seed=walsh(([-1]+[1]*(D0-1)) if seed_n else [1])
    N=T*T
    m=T.bit_length()-1
    histogram=Counter()
    private_masks=[walsh([-1]+[1]*7),walsh([1,0,0,1]),
                   walsh([1,0,1,1]),walsh([1,1,0,1])]
    for i,masks in enumerate(private_masks):
        t=T//2 if i==0 else T
        block=choice_tree(t,Bits(),seed,seed_n,False)
        scales=Counter(z*c for c in block.values() for z in masks if z)
        for scale,number in scales.items():
            row=[0]*4
            row[i]=scale
            histogram[tuple(row)]+=number
    for y in range(N):
        u=N-2 if y==0 else -2
        v=T*(-1)**(((y%T)&(y//T)).bit_count()%2)
        for q in range(4):
            for z in seed:
                if z:
                    a,b=u*z*(-1)**(q&1),v*z*(-1)**(q>>1)
                    histogram[0,a+b,a,b]+=1
    return {''.join(map(str,idx)):sum(count*sum(c!=0 for c in walsh(row[i] for i in idx))
                                      for row,count in histogram.items()) for idx in
            [(0,1,2,3),*((i,) for i in range(4)),*combinations(range(4),2)]}


def run(T,seed_n,compressed=False):
    rows,seed=branch_rows(T,seed_n,compressed)
    m=T.bit_length()-1
    seed_support=[i for i,c in enumerate(seed) if c]
    s0=len(seed_support)
    R0=rank(2*i+1 for i in seed_support)
    e=R0+1
    k=(0,2,4,4)
    records={}
    pieces=[(0,1,2,3),*((i,) for i in range(4)),*combinations(range(4),2)]
    for idx in pieces:
        support=gather(rows,idx)
        observed_rank=rank(2*i+1 for i in support)
        if len(idx)==4:
            sz=60
            D=(7*T//2)*e+2*m+R0+4
            only_A={tag for tag in support if (tag&~3) in rows[0]}
            inter=rank(only_A)+rank(set(support)-only_A)-rank(support)
            assert inter==2
        elif idx==(0,):
            sz,D=2,(T//2)*e+3
        elif len(idx)==1:
            sz,D=4+k[idx[0]],T*e+2*m+R0+2
        elif idx[0]==0:
            sz,D=2*(6+k[idx[1]]),(3*T//2)*e+2*m+R0+5
        else:
            sz,D=8+2*(k[idx[0]]+k[idx[1]]),2*T*e+2*m+R0+2
        size=sz*T*T*s0
        assert observed_rank==D-1, (idx,D,observed_rank)
        if not compressed:
            assert len(support)==size
        records[''.join(map(str,idx))]={"s":size,"D":D,"sampled_columns":len(support),
                                       "P_squared":str(Fraction(D*D,size))}
    if compressed:
        counts=full_cardinality_stream(T,seed_n)
        assert all(counts[k]==v["s"] for k,v in records.items())
    parent=Fraction(records['0123']['P_squared'])
    best=max(Fraction(row['P_squared']) for name,row in records.items() if name!='0123')
    return {"T":T,"seed_dimension":seed_n,"seed_s":s0,"seed_R":R0,"intersection_dimension":inter,
            "mode":"full streamed Walsh support counts and sampled generating columns" if compressed else "full exact spectra",
            "all_interface_restrictions_strictly_lose":parent>best,"records":records}


if __name__=="__main__":
    cases=[run(8,0),run(16,0),run(32,0),run(8,3),run(256,0,True)]
    assert cases[-1]["all_interface_restrictions_strictly_lose"]
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,"cases":cases},indent=2))
