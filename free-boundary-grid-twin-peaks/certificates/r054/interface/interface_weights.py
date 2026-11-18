"""Exact color-bridge kernels and gates for relative-rank seam conditions."""
from argparse import ArgumentParser
from collections import Counter,defaultdict
from decimal import Decimal,localcontext
from fractions import Fraction
from itertools import combinations
import hashlib
import json
from math import comb,factorial,isqrt,prod
from pathlib import Path
import random
import sys


def dec(value):
    with localcontext() as c:
        c.prec=57
        v=Fraction(value)
        return str(Decimal(v.numerator)/Decimal(v.denominator))


def layouts(n):
    assert n%2==0
    vertices=range(n*n)
    adjacency=[{j for j in vertices if abs(j%n-i%n)+abs(j//n-i//n)==1} for i in vertices]
    return adjacency,[{v for v in vertices if (v%n<n//2)==left} for left in (True,False)]


def local_peaks(words,adj):
    ranks=[{v:j for j,v in enumerate(w,1)} for w in words]
    return [{v for v in ranks[i] if all(ranks[i][v]<ranks[i][w] for w in adj[v]&ranks[i].keys())}
            for i in (0,1)]


def bridge_exact(words,adj,bound=None):
    """Count complete bridges by exact global peak bitset, optionally in a corridor."""
    first,second=words
    n,m=len(first),len(second)
    owner=[set(w) for w in words]
    peaks=local_peaks(words,adj)
    ranks=[{v:j for j,v in enumerate(w,1)} for w in words]
    displayed=[]
    for i,w in enumerate(words):
        displayed.append([(v in peaks[i], min((ranks[1-i][u] for u in adj[v]&owner[1-i]),default=len(words[1-i])+1)) for v in w])
    table=[[defaultdict(int) for _ in range(m+1)] for _ in range(n+1)]
    table[0][0][0]=1
    for a in range(n+1):
        for b in range(m+1):
            if a<n and (bound is None or abs(a+1-b)<=bound):
                trend,threshold=displayed[0][a]
                flag=1<<first[a] if trend and b<threshold else 0
                for t,k in table[a][b].items():table[a+1][b][t|flag]+=k
            if b<m and (bound is None or abs(a-b-1)<=bound):
                trend,threshold=displayed[1][b]
                flag=1<<second[b] if trend and a<threshold else 0
                for t,k in table[a][b].items():table[a][b+1][t|flag]+=k
    return dict(table[n][m])


def direct_peak_mask(word,adj):
    seen={v:k for k,v in enumerate(word)}
    return sum(1<<v for v in word if all(seen[v]<seen[u] for u in adj[v]))


def direct_shuffles(words,adj):
    L,R=words
    n,m=len(L),len(R)
    out=Counter()
    for left in combinations(range(n+m),n):
        pos=set(left)
        a,b=iter(L),iter(R)
        word=tuple(next(a) if d in pos else next(b) for d in range(n+m))
        out[direct_peak_mask(word,adj)]+=1
    return dict(out)


def finite_gate(seed=2026111854):
    adj,blocks=layouts(4)
    rng=random.Random(seed)
    cases=[]
    for _ in range(6):
        words=[rng.sample(sorted(H),len(H)) for H in blocks]
        dp=bridge_exact(words,adj)
        assert dp==direct_shuffles(words,adj)
        assert sum(dp.values())==comb(16,8)
        cases.append({'peak_masks':len(dp),'exact_two_peak_count':sum(c for m,c in dp.items() if m.bit_count()==2)})
    return {'seed':seed,'cases':cases,'literal_color_words_checked':len(cases)*comb(16,8)}


def hit_count(M,d):
    return comb(2*M,M+d) if d<=M else 0


def reflection_gate(M):
    seen_higher=[0]*(M+2)
    seen_lower=[0]*(M+2)
    survived=[0]*(M+1)
    for left in combinations(range(2*M),M):
        pos=set(left)
        path=mx=mn=0
        for j in range(2*M):
            path+=1 if j in pos else -1
            mx=max(mx,path);mn=min(mn,path)
        assert path==0
        for d in range(1,M+2):
            seen_higher[d]+=mx>=d
            seen_lower[d]+=-mn>=d
        for width in range(M+1):survived[width]+=mx<=width and -mn<=width
    assert all(seen_higher[d]==seen_lower[d]==hit_count(M,d) for d in range(1,M+2))
    assert all(survived[k]>=comb(2*M,M)-2*hit_count(M,k+1) for k in range(M+1))
    return {'M':M,'bridges':comb(2*M,M),'corridor_counts':survived}


def ordered_fraction(r,M):
    return Fraction(factorial(M)**r*prod(factorial(i) for i in range(r)),prod(factorial(M+i) for i in range(r)))


def direct_multicolor(r,M):
    """Enumerate multicolor words and test both count constraints literally."""
    state=[0]*r
    all_words=0
    count_ordered=0
    gap_counts=[0]*(M+1)
    def visit(t,good,width):
        nonlocal all_words,count_ordered
        if t==r*M:
            all_words+=1
            count_ordered+=good
            for k in range(width,M+1):gap_counts[k]+=1
            return
        for c in range(r):
            if state[c]==M:continue
            state[c]+=1
            visit(t+1,good and all(a>=b for a,b in zip(state,state[1:])),max(width,max(state)-min(state)))
            state[c]-=1
    visit(0,True,0)
    assert all_words==factorial(r*M)//factorial(M)**r
    assert Fraction(count_ordered,all_words)==ordered_fraction(r,M)
    for k,v in enumerate(gap_counts):
        assert Fraction(v,all_words)>=1-2*comb(r,2)*Fraction(hit_count(M,k+1),comb(2*M,M))
    return {'r':r,'M':M,'all_words':all_words,'ordered_count':count_ordered,'corridor_counts':gap_counts}


def grow(initial,region,adj):
    word=list(initial)
    seen=set(word)
    assert seen<=region and len(seen)==len(word)
    while seen!=region:
        front=set().union(*(adj[v] for v in seen))&region-seen
        assert front
        v=min(front)
        word.append(v);seen.add(v)
    return word


def same_mask_obstruction(n):
    h=n//2
    adj,(L,R)=layouts(n)
    u,p,v,q=h-1,0,h,h+1
    assert p not in adj[u] and v in adj[q]
    left=grow([u,p],L,adj)
    early=grow([q,v],R,adj)
    late=grow([q],R-{v},adj)+[v]
    assert local_peaks([left,early],adj)==local_peaks([left,late],adj)==[{u,p},{q}]
    chosen=(1<<p)|(1<<q)
    M=len(L)
    weights=[]
    for right in (early,late):
        exact=bridge_exact([left,right],adj)
        assert sum(exact.values())==comb(2*M,M)
        local_rank=right.index(v)+1
        expected=comb(2*M-local_rank,M)
        assert exact[chosen]==expected and all(mask.bit_count()!=2 or mask==chosen for mask in exact)
        weights.append({'killer_local_rank':local_rank,'accepted_words':str(expected),
                        'acceptance_exact':str(Fraction(expected,comb(2*M,M))),
                        'acceptance_decimal':dec(Fraction(expected,comb(2*M,M)))})
    return {'n':n,'M':M,'fixed_peak_pair':[p,q],'distance':abs(p-q),'same_induced_peaks':[[u,p],[q]],
            'fixed_cross_pair':[u,v],'two_word_pairs':weights,
            'early_over_late_ratio':str(Fraction(int(weights[0]['accepted_words']),int(weights[1]['accepted_words'])))}


def reciprocal_slack_example(n):
    h=n//2
    adj,(L,R)=layouts(n)
    # Left root on top, right root on bottom. Each early seam connector
    # faces the late unwanted peak of the opposite block.
    u=(n-1)*n+h-1
    z=h
    p=h-2
    q=(n-1)*n+h+1
    c=h-1
    d=(n-1)*n+h
    words=[]
    for P,r,active,killed in ((L,p,c,u),(R,q,d,z)):
        safe=P-(adj[killed]|{killed})
        prefix=grow([r,active],safe,adj)
        words.append(grow(prefix+[killed],P,adj))
    assert local_peaks(words,adj)==[{p,u},{q,z}]
    i,j=words[0].index(u)+1,words[1].index(z)+1
    M=len(L)
    D=isqrt(4*M)+(isqrt(4*M)**2<4*M)
    assert min(i,j)>=2+D
    chosen=(1<<p)|(1<<q)
    whole=bridge_exact(words,adj)
    corridor=bridge_exact(words,adj,D)
    assert set(corridor)=={chosen}
    N=comb(2*M,M)
    accepted=whole[chosen]
    assert 0<sum(corridor.values())<=accepted<=N
    assert sum(corridor.values())>=N-2*hit_count(M,D+1)
    return {'n':n,'M':M,'undesired_peak_ranks':[i,j],'suppressor_ranks':[2,2],
            'corridor_width':D,'protection_anchors':[p,q],
            'global_anchor_manhattan':abs(p%n-q%n)+abs(p//n-q//n),
            'accepted_exact':str(Fraction(accepted,N)),
            'accepted_decimal':dec(Fraction(accepted,N)),
            'corridor_exact':str(Fraction(sum(corridor.values()),N)),
            'certified_corridor_lower':str(Fraction(N-2*hit_count(M,D+1),N))}


def main():
    p=ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
    out={'scope':'Exact interleaving counts for fixed local words, no claim on their normalized mass under K=2',
         'python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'original_grid_bruteforce_gate':finite_gate(),'reflection_gates':[reflection_gate(M) for M in range(1,7)],
         'multicolor_gates':[direct_multicolor(r,M) for r,M in ((2,4),(3,3),(4,2),(4,3))],
         'same_mask_examples':[same_mask_obstruction(n) for n in (6,8,10,12)],
         'two_way_slack_example':reciprocal_slack_example(8)}
    with a.output.open('x') as f:
        json.dump(out,f,indent=2);f.write('\n')
    print('all fixed-word gates passed; n12 early/late',out['same_mask_examples'][-1]['early_over_late_ratio'])
    print('two-way example acceptance',out['two_way_slack_example']['accepted_decimal'])

if __name__=='__main__':main()
