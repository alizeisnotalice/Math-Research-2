#!/usr/bin/env python3
"""Full 3x3 encode/decode audit for the r033 parent-path records."""
from collections import Counter
from itertools import permutations
import json
from pathlib import Path


def make_adj(n=3):
    ans = []
    for y in range(n):
        for x in range(n):
            ans.append({yy*n+xx for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))
                        if 0 <= (xx:=x+dx) < n and 0 <= (yy:=y+dy) < n})
    return ans


ADJ = make_adj()


def births(word):
    seen, out = set(), []
    for pos,v in enumerate(word):
        if not (ADJ[v] & seen):
            out.append((pos,v))
        seen.add(v)
    return out


def component(mask, root):
    seen, todo = {root}, [root]
    while todo:
        v = todo.pop()
        for z in ADJ[v] & mask:
            if z not in seen:
                seen.add(z); todo.append(z)
    return seen


def encode(word):
    bs = births(word)
    assert len(bs) == 2
    p, r = word[0], bs[1][1]
    seen, parent = set(), {}
    q = w = None
    for pos,v in enumerate(word):
        earlier = ADJ[v] & seen
        if earlier:
            parent[v] = min(earlier)
        if pos > bs[1][0] and q is None:
            cp = component(seen,p)
            if r in seen and (ADJ[v]&cp) and (ADJ[v]&(seen-cp)):
                q,w = pos,v
                break
        seen.add(v)
    assert q is not None
    pre = set(word[:q]); cp = component(pre,p); cr = pre-cp

    def chain(z,root):
        out=[z]
        while out[-1] != root:
            out.append(parent[out[-1]])
        return list(reversed(out))

    pa = chain(min(ADJ[w]&cp),p)
    rb = chain(min(ADJ[w]&cr),r)
    path = pa + [w] + list(reversed(rb))
    ell, j = len(path), len(pa)
    assert len(path) == len(set(path)) and 1 <= j <= ell-2 and 3 <= ell <= q+1
    pset=set(path)
    output=tuple(path+[v for v in word if v not in pset])

    pin=list(path[1:j]); rorig=list(reversed(path[j+1:]))
    posa=tuple(i for i in range(1,q) if word[i] in pset)
    bita=tuple("P" if word[i] in pin else "R" for i in posa)
    assert len(posa)==ell-2 and bita.count("P")==len(pin) and bita.count("R")==len(rorig)
    rec_a=(q,ell,j,posa,bita)

    rec_b=None
    if word[1] == r:
        posb=tuple(i for i in range(2,q) if word[i] in pset)
        rin=rorig[1:]
        bitb=tuple("P" if word[i] in pin else "R" for i in posb)
        assert len(posb)==ell-3 and bitb.count("P")==len(pin) and bitb.count("R")==len(rin)
        rec_b=(q,ell,j,posb,bitb)
    return q,output,rec_a,rec_b


def decode_a(out,rec):
    q,ell,j,positions,bits=rec
    p,w=out[0],out[j]
    ip,ir=iter(out[1:j]),iter(reversed(out[j+1:ell]))
    marked=dict(zip(positions,bits)); other=iter(out[ell:]); pre=[p]
    for i in range(1,q):
        pre.append(next(ip) if marked.get(i)=="P" else next(ir) if marked.get(i)=="R" else next(other))
    return tuple(pre+[w]+list(other))


def decode_b(out,rec):
    q,ell,j,positions,bits=rec
    p,r,w=out[0],out[ell-1],out[j]
    ip,ir=iter(out[1:j]),iter(reversed(out[j+1:ell-1]))
    marked=dict(zip(positions,bits)); other=iter(out[ell:]); pre=[p,r]
    for i in range(2,q):
        pre.append(next(ip) if marked.get(i)=="P" else next(ir) if marked.get(i)=="R" else next(other))
    return tuple(pre+[w]+list(other))


def main():
    A2=B=0; aq=Counter(); bq=Counter(); ima=set(); imb=set()
    for word in permutations(range(9)):
        if len(births(word)) != 2:
            continue
        A2 += 1
        q,out,ra,rb=encode(word)
        assert len(births(out)) == 1
        assert decode_a(out,ra) == word and (out,ra) not in ima
        ima.add((out,ra)); aq[q]+=1
        if rb is not None:
            B += 1
            assert decode_b(out,rb) == word and (out,rb) not in imb
            imb.add((out,rb)); bq[q]+=1
    assert (A2,B,len(ima),len(imb)) == (143808,81984,143808,81984)
    ans={"status":"passed","A2":A2,"B_top2":B,"A2q":dict(sorted(aq.items())),
         "Bq":dict(sorted(bq.items())),"injective_A2_records":len(ima),
         "injective_B_records":len(imb)}
    Path(__file__).with_name("path_encoder_n3_audit.json").write_text(json.dumps(ans,sort_keys=True,indent=2)+"\n")
    print(json.dumps(ans,sort_keys=True,indent=2))


if __name__ == "__main__":
    main()
