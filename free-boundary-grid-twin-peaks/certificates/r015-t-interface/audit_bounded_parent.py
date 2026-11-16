#!/usr/bin/env python3
"""Integer audits for the bounded-parent gate construction.

The first audit checks the hook-length expression for the gate-and-leaf tree.
The second constructs the proposed family on the 3x3 T partition, enumerates
all constructed words, and independently checks that each has exactly two
births and that no two construction parameters yield the same word.
"""
from itertools import combinations, permutations
from math import ceil, comb, factorial
from pathlib import Path
import json


def extensions(k, delta):
    q=ceil(k/delta); r=k-delta*(q-1)
    # nodes: gates 0,...,q-1 and k named leaves.  Gate j precedes its group;
    # gate j precedes gate j+1.  Brute force all topological growth words.
    nodes=[('g',j) for j in range(q)]
    groups=[];cursor=0
    for j in range(q):
        size=r if j==0 else delta
        group=[('c',cursor+i) for i in range(size)]
        cursor+=size;groups.append(group);nodes+=group
    predecessors={x:set() for x in nodes}
    for j in range(q):
        for x in groups[j]:predecessors[x].add(('g',j))
        if j:predecessors[('g',j)].add(('g',j-1))
    dp={frozenset():1}
    for _ in nodes:
        nxt={}
        for seen,z in dp.items():
            for x in nodes:
                if x not in seen and predecessors[x] <= seen:
                    key=seen|{x};nxt[key]=nxt.get(key,0)+z
        dp=nxt
    got=dp[frozenset(nodes)]
    expected=factorial(k+q-1)//((delta+1)**(q-1)*factorial(q-1))
    assert got==expected
    # q-1 is a binomial mode; this verifies the stated (weaker) inequality.
    lower_num=(delta+1)**max(0,k-delta)
    lower_den=(delta**max(0,k-delta))*(k+q)
    assert got*lower_den >= factorial(k)*lower_num
    return {'k':k,'delta':delta,'q':q,'r':r,'E':got,
            'matching_delta_1': delta==1,
            'stated_E_over_factorial_lower_bound_verified':True}


def grid3():
    n=3;V={(x,y) for x in range(n) for y in range(n)}
    A={(0,0),(0,1)};B={(2,0),(2,1)};C=V-A-B
    adj={v:set() for v in V}
    for x,y in V:
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            q=x+dx,y+dy
            if q in V:adj[x,y].add(q)
    return A,B,C,adj


def onepeak_words(S,adj):
    out=[]
    for word in permutations(S):
        seen=set()
        for v in word:
            if seen and not(adj[v]&seen):break
            seen.add(v)
        else:out.append(word)
    return out


def root_fixed_interleavings(x,y):
    # preserve each word and keep both roots as positions 1,2.
    for first,second in ((x,y),(y,x)):
        for slots in combinations(range(len(x)+len(y)-2),len(x)-1):
            slots=set(slots);ix=iy=1;out=[first[0],second[0]]
            for z in range(len(x)+len(y)-2):
                if z in slots:out.append(x[ix]);ix+=1
                else:out.append(y[iy]);iy+=1
            yield tuple(out)


def tree_extensions(gates,children,parent_rank):
    # Here k=2,delta=4,q=1, so gate is forced then two leaves are free.
    # Keep a generic recursion to audit the actual construction word.
    q=len(gates);k=len(children);delta=4;r=k-delta*(q-1)
    groups=[children[:r]]+[children[r+(j-1)*delta:r+j*delta] for j in range(1,q)]
    pred={x:set() for x in [*gates,*children]}
    for j,g in enumerate(gates):
        if j:pred[g].add(gates[j-1])
        for c in groups[j]:pred[c].add(g)
    nodes=list(pred);result=[]
    def rec(word,seen):
        if len(word)==len(nodes):result.append(tuple(word));return
        for x in nodes:
            if x not in seen and pred[x] <= seen:rec(word+[x],seen|{x})
    rec([],set());return result


def audit_grid3():
    A,B,C,adj=grid3();selected=((1,0),(1,1));parents={(1,0):(0,0),(1,1):(2,1)}
    suffixes=[]
    Z=A|B|set(selected)
    for word in permutations(C-set(selected)):
        seen=set(Z)
        for v in word:
            if not(adj[v]&seen):break
            seen.add(v)
        else:suffixes.append(tuple(word))
    common=[v for v in C if adj[v]&A and adj[v]&B]
    fnum=0
    for word in permutations(C):
        if word[0] not in common:continue
        seen=A|B
        for v in word:
            if not(adj[v]&seen):break
            seen.add(v)
        else:fnum+=1
    assert fnum <= factorial(len(C))//factorial(len(C)-len(selected))*len(suffixes)
    words=set();checked=0
    for aw in onepeak_words(A,adj):
      for bw in onepeak_words(B,adj):
       for ab in root_fixed_interleavings(aw,bw):
        ranks={v:i+1 for i,v in enumerate(ab)}
        child=tuple(sorted(selected,key=lambda x:(ranks[parents[x]],x)))
        # k=2, delta=4 gives q=1 and gate the final AB vertex.
        prefix=ab[:-1]; gates=(ab[-1],)
        for tail in tree_extensions(gates,child,ranks):
            for suffix in suffixes:
                word=prefix+tail+suffix
                births=sum(not(adj[v]&set(word[:i])) for i,v in enumerate(word))
                assert births==2
                assert word not in words
                words.add(word);checked+=1
    # a(A)=a(B)=2, M=4, root fixed AB count=2*C(2,1), E=2.
    base=2*2*2*comb(2,1)*2
    assert checked==base*len(suffixes)
    return {'base_gate_words':base,'h_Z':len(suffixes),'constructed_words':checked,
            'distinct_words':len(words),'Fnum':fnum,'falling_factorial_times_hZ':20*len(suffixes),
            'restriction_bound_verified':True,'all_have_exactly_two_births':True}


def audit_noncommon_grid3():
    # The selected bottom-left point only touches A (and C), not B.
    A,B,C,adj=grid3();selected=((0,2),);parents={(0,2):(0,1)}
    assert not (adj[selected[0]]&B) and adj[selected[0]]&A
    suffixes=[];Z=A|B|set(selected)
    for word in permutations(C-set(selected)):
        seen=set(Z)
        for v in word:
            if not(adj[v]&seen):break
            seen.add(v)
        else:suffixes.append(tuple(word))
    words=set();checked=0
    for aw in onepeak_words(A,adj):
      for bw in onepeak_words(B,adj):
       for ab in root_fixed_interleavings(aw,bw):
        prefix=ab[:-1];gate=(ab[-1],)
        for tail in tree_extensions(gate,selected,{parents[selected[0]]:4}):
            for suffix in suffixes:
                word=prefix+tail+suffix
                births=sum(not(adj[v]&set(word[:i])) for i,v in enumerate(word))
                assert births==2 and word not in words
                words.add(word);checked+=1
    base=2*2*2*comb(2,1)
    assert checked==base*len(suffixes)
    return {'selected_point':selected[0],'is_not_common_bridge':True,'base_gate_words':base,
            'h_Z':len(suffixes),'constructed_words':checked,'distinct_words':len(words),
            'all_have_exactly_two_births':True}


def audit_grid_constant_integer_rounding():
    """Exact finite audit of the discrete implications used in the corollary.

    The proof itself is symbolic for every g.  This checks that the threshold
    k >= ceil(63g/64) and q=ceil(k/4) cause no small-g exception.
    """
    cases=0
    for g in range(1,4097):
        for k in range((63*g+63)//64,g+1):
            s=g-k;q=(k+3)//4
            assert 64*s <= g
            assert q <= k and k+q <= 2*g
            if s:
                assert g >= 64
            cases+=1
    return {'g_range':[1,4096],'integer_cases_checked':cases,
            'threshold_and_q_rounding_verified':True}


rows=[extensions(k,d) for d in range(1,5) for k in range(1,7)]
out={'status':'passed','tree_rows':rows,'grid3_actual_word_audit':audit_grid3(),
     'grid3_noncommon_child_audit':audit_noncommon_grid3(),
     'grid_constant_integer_rounding_audit':audit_grid_constant_integer_rounding()}
Path('bounded-parent-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
