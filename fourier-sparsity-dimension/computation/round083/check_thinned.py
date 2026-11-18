"""Exact small matching checks for local and multiple-leaf spectrum folding."""
from collections import Counter, defaultdict
from itertools import product
from random import Random
import json
import platform


def eliminate(rows):
    basis={}
    for v,s in rows:
        while v and v.bit_length() in basis:
            a,b=basis[v.bit_length()]
            v ^= a
            s ^= b
        if v:
            basis[v.bit_length()]=(v,s)
    return basis


def reduce(v,basis):
    sig=0
    for i in sorted(basis,reverse=True):
        if (v>>(i-1))&1:
            a,b=basis[i]
            v ^= a
            sig ^= b
    return v,sig


def fold(spec,basis):
    w=defaultdict(int)
    for v,a in spec.items():
        p,s=reduce(v,basis)
        w[p]+=(1-2*s)*a
    return {v:a for v,a in w.items() if a}


def rank(rows):
    return len(eliminate([(v,0) for v in rows]))


def rowspace(H,bits):
    out=[(0,0)]
    for v,s in zip(H,bits):
        out+= [(v^a,s^b) for a,b in out]
    return out


def audit():
    n=3;M=1<<n;N=M//2
    # Three two-bit bent seeds + common dummy on the last edge.
    edges=[(0,1),(2,4),(3,7),(5,6)]
    qcount=n+3*3+1
    aug=1<<qcount
    data_mask=~(M-1)
    masks=[7<<(n+3*i) for i in range(3)]+[1<<(qcount-1)]
    guards=[1<<(n+3*i+2) for i in range(3)]+[masks[-1]]
    height=sum(guards)
    spectra=[]; seeds=[]
    for i,(a,b) in enumerate(edges):
        seed={guards[i]:2} if i==3 else {
            guards[i]|(j<<(n+3*i)):1 if j!=3 else -1 for j in range(4)}
        seeds.append(seed)
        sp={}
        for outer in range(M):
            if ((outer&(a^b)).bit_count()%2)==0:
                for core,c in seed.items():
                    sp[core|outer]=c*(1-2*((outer&a).bit_count()%2))
        spectra.append(sp)
    S={v:a for sp in spectra for v,a in sp.items()}
    r=3;s=4;D=1+rank(aug|v for v in S)
    assert D==14 and len(S)==52
    # Verify the parent representation independently on all 8,192 inputs.
    for x in range(1<<qcount):
        u=x&(M-1)
        i=next(i for i,edge in enumerate(edges) if u in edge)
        e=x>>(n+3*i)
        tg=((x&guards[i]).bit_count()%2)^int(i!=3 and e&3==3)
        assert sum(c*(1-2*((v&x).bit_count()%2)) for v,c in S.items())==8*(1-2*tg)
    effective=[v for v in range(1,aug) if not((v&height).bit_count()%2)]
    rng=Random(8303)
    tests=set()
    # Sample all ways a queried variable may be supported on one leaf, two
    # leaves, a pure address, or generic data. Keep explicit mixed addresses.
    elementary=set(range(1,M))
    for sp in seeds:
        for a in sp:
            for b in sp:
                if a!=b:
                    elementary.update(a^b^i for i in range(M))
    labels=list(set.union(*(set(sp) for sp in seeds)))
    for a in labels:
        for b in labels:
            if a!=b:
                elementary.update(a^b^i for i in range(M))
    el=sorted(elementary)
    tests.update((v,) for v in el)
    for h in (2,3):
        for _ in range(1000):
            pool=el if rng.getrandbits(1) else effective
            H=tuple(sorted(rng.sample(pool,h)))
            if rank(H)==h: tests.add(H)
    tally=Counter()
    for H in sorted(tests):
        h=len(H)
        for ans in product((0,1),repeat=h):
            allrows=rowspace(H,ans)
            pure=[(v,a) for v,a in allrows if not(v&data_mask)]
            Bs=eliminate(list(zip(H,ans)))
            fused=fold(S,Bs)
            newD=rank(aug|v for v in fused)+1
            if len(pure)>1:
                C=[x for x in range(M) if all((x&v).bit_count()%2==a for v,a in pure)]
                f=sum(a in C and b in C for a,b in edges[:3])
                t=sum((a in C)^(b in C) for a,b in edges[:3])
                assert 2*len(fused)>=len(C)*((f+2*t-2*h)*s-8*(1<<h))
                assert newD<=(f+t)*r+n+2
                tally['pure_address_component']+=1
                continue
            local=[]; shifts=[]
            for i,mask in enumerate(masks):
                arr=[(v,a) for v,a in allrows if v&data_mask and not(v&data_mask&~mask)]
                local.extend(arr)
                if i==3: continue
                edge_delta=edges[i][0]^edges[i][1]
                keep=[(v&data_mask,a^((v&edges[i][0]).bit_count()%2)) for v,a in arr
                      if ((v&edge_delta).bit_count()%2)==0]
                cs=fold(seeds[i],eliminate(keep))
                shifts.append((r-rank(aug|v for v in cs),len(cs)))
                assert len(fold(spectra[i],eliminate(arr)))==N*len(cs)
            Bl=eliminate(local)
            localspec=[fold(sp,Bl) for sp in spectra]
            projlabels=[]
            for sp in localspec:
                projlabels.append({reduce(f,Bs)[0]&data_mask for f in sp})
            collision=int(any(projlabels[i]&projlabels[j]
                              for i in range(4) for j in range(i)))
            X=sum(x for x,_ in shifts)
            assert newD<=D-X-collision
            assert len(fused)>=sum(map(len,localspec))-4*N*(1<<h)*collision
            assert all(si*r*r>=s*(r-mi)**2 for mi,si in shifts)
            if not X and not collision:
                assert len(fused)>=len(S) and newD<=D-h
            tally[f'no_pure_X{int(X>0)}_collision{collision}']+=1
    return dict(seed=8303,parents_full_truth_entries=1<<qcount,query_spaces=len(tests),
                answer_counts=dict(sorted(tally.items())))


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
                          finite_check=audit()),indent=2))
