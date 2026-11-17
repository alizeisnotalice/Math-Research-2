#!/usr/bin/env python3
"""Exact FWHT and small legal-direction enumeration for R29."""
from __future__ import annotations

import hashlib, itertools, json, platform, sys, time
from pathlib import Path

OUT = Path(__file__).with_name("r29_counterexamples_and_permutations_result.json")

def par(x: int) -> int: return x.bit_count() & 1

def fwht(a: list[int]) -> list[int]:
    a=a[:]; h=1
    while h<len(a):
        for s in range(0,len(a),2*h):
            for j in range(h): a[s+j],a[s+j+h]=a[s+j]+a[s+j+h],a[s+j]-a[s+j+h]
        h*=2
    return a

def rank(vs: list[int]) -> int:
    b={}
    for v in vs:
        while v:
            p=v.bit_length()-1
            if p in b: v^=b[p]
            else: b[p]=v; break
    return len(b)

def vector(k: int, fun) -> int:
    q=1<<k
    return sum(fun(z,w)<<(z*q+w) for z in range(q) for w in range(q))

def all_bent(k: int, gs) -> tuple[bool, int | None]:
    q=1<<k
    for mask in range(1<<len(gs)):
        a=[]
        for z in range(q):
            for w in range(q):
                bit=par(z&w)
                for j,g in enumerate(gs):
                    if (mask>>j)&1: bit^=g(z,w)
                a.append(1 if bit==0 else -1)
        if not all(abs(x)==q for x in fwht(a)): return False,mask
    return True,None

def block_alpha(k: int, g, w: int) -> list[int]:
    q=1<<k; ans=[]
    for alpha in range(q):
        total=sum(1 if par(z&w)^g(z,w)^par(alpha&z)==0 else -1 for z in range(q))
        if total:
            assert abs(total)==q; ans.append(alpha)
    return ans

def two_block(k: int) -> dict[str,object]:
    gs=[lambda z,w:int(w==0),lambda z,w:int(w==1),lambda z,w:(z&1) if w in (0,1) else 0]
    ok,bad=all_bent(k,gs)
    zero=lambda z,w:0
    before=[block_alpha(k,zero,w) for w in (0,1)]
    after=[block_alpha(k,gs[2],w) for w in (0,1)]
    return {"k":k,"q":1<<k,"rank":rank([vector(k,g) for g in gs]),"phases":8,
            "all_bent":ok,"failure":bad,"R_size":2*(1<<k),"local_budget_sum":2,
            "partial_alpha_before":before,"partial_alpha_after_third":after,"supports_change":before!=after}

def r4_case(k: int) -> dict[str,object]:
    E=(0,1,2,3)
    gs=[(lambda t:(lambda z,w:int(w==t)))(t) for t in E]
    gs += [lambda z,w:(z&1) if w in E else 0,
           lambda z,w:((z>>1)&1) if w in (0,2) else 0,
           lambda z,w:((z>>1)&1) if w in (1,3) else 0]
    ok,bad=all_bent(k,gs); images=[]
    for mask in range(8):
        image=[]
        for w in E:
            out=w
            if mask&1: out^=1
            if mask&2 and not(w&1): out^=2
            if mask&4 and (w&1): out^=2
            image.append(out)
        assert sorted(image)==list(E); images.append(image)
    return {"k":k,"q":1<<k,"rank":rank([vector(k,g) for g in gs]),"phases":128,
            "all_bent":ok,"failure":bad,"slope_images":images,"sharp_2r_minus_1":7}

def legal_direction_subspaces(S: tuple[int,...]) -> dict[str,object]:
    k=max(S).bit_length(); r=len(S)
    legal=set()
    for perm in itertools.permutations(S):
        code=sum((x^y)<<(k*i) for i,(x,y) in enumerate(zip(S,perm)))
        legal.add(code)
    # Breadth-first closure: every retained frozenset is a linear subspace contained in legal.
    todo=[frozenset({0})]; seen={frozenset({0})}; dims={0:1}
    while todo:
        sub=todo.pop()
        for v in legal-set(sub):
            ext=frozenset(set(sub)|{x^v for x in sub})
            if ext <= legal and ext not in seen:
                seen.add(ext); todo.append(ext); dims[(len(ext)).bit_length()-1]=dims.get((len(ext)).bit_length()-1,0)+1
    return {"S":list(S),"ambient_coordinate_bits":k,"permutation_count":len(legal),
            "legal_linear_subspace_count":len(seen),"dimension_histogram":{str(d):n for d,n in sorted(dims.items())},
            "max_dimension":max((len(sub)).bit_length()-1 for sub in seen)}

def main() -> None:
    t=time.perf_counter()
    sets=[(0,1,2),(0,1,2,3),(0,1,2,4),(0,1,2,4,8)]
    ans={"scope":"Exact checks of displayed constructions and four finite permutation sets only.",
         "python":sys.version,"platform":platform.platform(),
         "two_block_examples":[two_block(k) for k in range(2,7)],
         "r4_parallel_examples":[r4_case(k) for k in (3,4)],
         "small_permutation_direction_enumeration":[legal_direction_subspaces(S) for S in sets],
         "elapsed_seconds":time.perf_counter()-t,
         "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(ans,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");print(json.dumps(ans,ensure_ascii=False,indent=2))
    assert all(x["all_bent"] and x["rank"]==3 and x["supports_change"] for x in ans["two_block_examples"])
    assert all(x["all_bent"] and x["rank"]==7 for x in ans["r4_parallel_examples"])
    assert [x["max_dimension"] for x in ans["small_permutation_direction_enumeration"]]==[1,3,2,2]
if __name__=='__main__': main()
