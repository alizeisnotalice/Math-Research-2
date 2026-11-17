#!/usr/bin/env python3
"""Exact R27 checks: D4 twin classes and k=3 pair-block partial Walsh sums."""
from __future__ import annotations

import hashlib, itertools, json, platform, sys, time
from pathlib import Path

OUT = Path(__file__).with_name("twin_closure_and_pairings_result.json")

def par(x: int) -> int: return x.bit_count() & 1

def fwht(a: list[int]) -> list[int]:
    a = a[:]; h = 1
    while h < len(a):
        for s in range(0, len(a), 2*h):
            for j in range(h): a[s+j], a[s+j+h] = a[s+j]+a[s+j+h], a[s+j]-a[s+j+h]
        h *= 2
    return a

def rank(words: list[int]) -> int:
    bas: dict[int,int] = {}
    for word in words:
        while word:
            p = word.bit_length()-1
            if p in bas: word ^= bas[p]
            else: bas[p] = word; break
    return len(bas)

def code_span(gen: list[int]) -> list[int]:
    out=[0]
    for word in gen: out += [x ^ word for x in out]
    return out

def twin_classes(q: int, generators: list[int]) -> list[list[int]]:
    ans: dict[tuple[int,...], list[int]] = {}
    for w in range(q): ans.setdefault(tuple((d>>w)&1 for d in generators), []).append(w)
    return sorted(ans.values(), key=lambda c:(c[0], len(c)))

def rm1_orthogonal(k: int, generators: list[int]) -> bool:
    q=1<<k
    return all(par(d & sum(1<<w for w in range(q) if par(a&w)^c)) == 0
               for d in generators for a in range(q) for c in (0,1))

def connected(q: int, edges: list[int]) -> int:
    p=list(range(q))
    def find(x: int) -> int:
        while p[x]!=x: p[x]=p[p[x]]; x=p[x]
        return x
    for edge in edges:
        vs=[w for w in range(q) if (edge>>w)&1]
        for w in vs[1:]:
            x,y=find(vs[0]),find(w)
            if x!=y:p[y]=x
    return len({find(w) for w in range(q)})

def fixed_pair_family(k: int) -> dict[str,object]:
    q=1<<k; m=q//2
    pairs=[(1<<(2*i))|(1<<(2*i+1)) for i in range(m)]
    gen=[pairs[0]^pairs[i] for i in range(1,m)]
    D=code_span(gen); wt4=[d for d in D if d.bit_count()==4]
    return {"k":k,"q":q,"D_dim":rank(gen),"D_size_enumerated":len(D),
            "D4_rank":rank(wt4),"D4_weight4_count":len(wt4),"H_dim":q-rank(gen),
            "expected_H_dim":m+1,"RM1_in_H":rm1_orthogonal(k,gen),
            "twin_classes":twin_classes(q,gen),"hypergraph_components":connected(q,wt4)}

def prescribed_k4_models() -> list[dict[str,object]]:
    models=[]
    for name, sets in [
        ("six_point_even", ((0,1,2,3),(0,1,4,5))),
        ("odd_twin_obstruction", ((0,1,2,3),(0,4,8,12))),
    ]:
        gen=[sum(1<<w for w in S) for S in sets]; D=code_span(gen); wt4=[d for d in D if d.bit_count()==4]
        models.append({"name":name,"D_dim":rank(gen),"D4_rank":rank(wt4),"D4_words":wt4,
                       "RM1_in_H":rm1_orthogonal(4,gen),"twin_classes":twin_classes(16,gen),
                       "hypergraph_components":connected(16,wt4)})
    return models

def bent4_masks() -> list[int]:
    ans=[]
    for mask in range(1<<16):
        if all(abs(x)==4 for x in fwht([1 if ((mask>>i)&1)==0 else -1 for i in range(16)])): ans.append(mask)
    return ans

def pair_block_test() -> dict[str,object]:
    q=8; choices=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
    masks=bent4_masks(); cases=nonzero=0; hist:dict[int,int]={}; failure=None
    for mask in masks:
        sign=[[0]*q for _ in range(q)]
        for z in range(q):
            for w in range(q):
                bit=((mask>>((z&3)|((w&3)<<2)))&1)^(((z>>2)&1)&((w>>2)&1))
                sign[z][w]=1 if bit==0 else -1
        for c0,c1 in itertools.product(choices,repeat=2):
            cases+=1; cover=[0]*64
            blocks=[tuple(lo|(hi<<2) for lo in pair) for hi,C in enumerate((c0,c1)) for pair in C]
            for pair in blocks:
                for a in range(q):
                    for b in range(q):
                        t=sum(sign[z][w]*(-1 if par(a&z) else 1)*(-1 if par(b&w) else 1) for z in range(q) for w in pair)
                        if t not in (0,-q,q): failure={"mask":mask,"pair":pair,"alpha":a,"beta":b,"value":t}; break
                        if t: cover[a*q+b]+=1; nonzero+=1
                    if failure: break
                if failure: break
            if failure: break
            for x in cover: hist[x]=hist.get(x,0)+1
            if any(x!=1 for x in cover): failure={"mask":mask,"reason":"coverage","coverage":cover}; break
        if failure: break
    return {"construction":"g(z,w)=F(z_low,w_low)+z_high*w_high","bent_F_count":len(masks),
            "pairing_choices_per_component":3,"cases_checked":cases,"blocks_per_case":4,
            "allowed_values":[0,-8,8],"nonzero_block_frequency_count":nonzero,
            "coverage_histogram":hist,"all_passed":failure is None,"first_failure":failure}

def main() -> None:
    t=time.perf_counter()
    ans={"scope":"Exact finite checks only, not a general theorem.","python":sys.version,"platform":platform.platform(),
         "fixed_a_one_pair_codes":[fixed_pair_family(k) for k in (3,4,5)],
         "rm1_k3_twin_classes":twin_classes(8,[d for d in range(256) if all(par(d & sum(1<<w for w in range(8) if par(a&w)^c))==0 for a in range(8) for c in (0,1))]),
         "k4_models":prescribed_k4_models(),"k3_pair_block_walsh":pair_block_test(),"elapsed_seconds":time.perf_counter()-t,
         "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(ans,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"); print(json.dumps(ans,ensure_ascii=False,indent=2))
    if not ans["k3_pair_block_walsh"]["all_passed"]: raise SystemExit(1)
if __name__=="__main__": main()
