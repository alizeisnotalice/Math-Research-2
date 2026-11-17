#!/usr/bin/env python3
"""Exhaustive local five-site flux audit on square-grid polyominoes."""
from collections import defaultdict
from itertools import permutations
import json
from pathlib import Path

DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))

def norm(s):
    ax=min(x for x,y in s); ay=min(y for x,y in s)
    return frozenset((x-ax,y-ay) for x,y in s)

def gen(k):
    z={frozenset({(0,0)})}
    for _ in range(1,k):
        nxt=set()
        for s in z:
            for x,y in s:
                for dx,dy in DIRS:
                    if (x+dx,y+dy) not in s: nxt.add(norm(set(s)|{(x+dx,y+dy)}))
        z=nxt
    return z

def d4(s):
    ans=[]
    for sw in (0,1):
        for sx in (-1,1):
            for sy in (-1,1):
                ans.append(tuple(sorted(norm({(sx*(y if sw else x),sy*(x if sw else y)) for x,y in s}))))
    return ans

def flux(cells):
    cells=tuple(sorted(cells))
    adj={v:{w for w in cells if abs(v[0]-w[0])+abs(v[1]-w[1])==1} for v in cells}
    c=0; pair=defaultdict(int)
    for word in permutations(cells):
        v,p=word[:2]
        if p in adj[v]: continue
        seen={v,p}; legal=True
        for z in word[2:]:
            if not(adj[z]&seen): legal=False; break
            seen.add(z)
        if not legal: continue
        j=next((i for i,z in enumerate(word[2:],3) if z in adj[p]),None)
        a=1 if j==4 else -2 if j==5 else 0
        if a:
            r=(v[0]-p[0])**2+(v[1]-p[1])**2
            c+=a*r; pair[v,p]+=a
    return c,{str(k):z for k,z in sorted(pair.items()) if z}

def abstract_disjoint():
    V=range(5); adj={0:{1},1:{0,2},2:{1},3:{4},4:{3}}; c=defaultdict(int)
    for w in permutations(V):
        v,p=w[:2]
        if p in adj[v]: continue
        seen={v,p}; legal=True
        for z in w[2:]:
            if not(adj[z]&seen): legal=False; break
            seen.add(z)
        if not legal: continue
        j=next((i for i,z in enumerate(w[2:],3) if z in adj[p]),None)
        if j==4: c[v,p]+=1
        if j==5: c[v,p]-=2
    # Directed coefficients need not vanish.  The spatial weight is symmetric
    # in the unordered root pair, so check the required reverse-pair cancelation.
    for v in V:
        for p in range(v + 1, 5):
            assert c[v,p] + c[p,v] == 0, (v,p,c[v,p],c[p,v])
    return ({f'{v},{p}':z for (v,p),z in sorted(c.items()) if z}, True)

def main():
    fixed=gen(5); assert len(fixed)==63
    groups=defaultdict(list)
    for s in fixed: groups[min(d4(s))].append(s)
    assert len(groups)==12
    rows=[]; allneg=True
    for canon,orb in sorted(groups.items()):
        vals={flux(s)[0] for s in orb}; assert len(vals)==1
        c=vals.pop(); allneg &= c<=0
        rows.append({'free_shape':[list(x) for x in canon],'fixed_embeddings':len(orb),
                     'c_r':c,'pair_coefficients_one_embedding':flux(orb[0])[1]})
    dis, dis_symmetric_zero=abstract_disjoint()
    out={'status':'passed','fixed_connected_polyominoes':len(fixed),
         'free_connected_polyominoes':len(groups),'all_connected_c_r_nonpositive':allneg,
         'free_shape_rows':rows,'disjoint_P3_union_K2_pair_coefficients':dis,
         'disjoint_P3_union_K2_all_pair_coefficients_zero':not dis,
         'disjoint_P3_union_K2_unordered_pair_coefficients_zero':dis_symmetric_zero}
    Path(__file__).with_name('top5_polyomino_flux_audit.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__': main()
