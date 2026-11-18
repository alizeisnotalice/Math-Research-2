from fractions import Fraction
from math import factorial
from random import Random


def last_row_factor(a, b, c, safe_only_right=True, exists=None, safe=()):
    W=len(b)
    if exists is None:
        exists=(True,)*W
    safe=set(safe)
    if safe_only_right:
        safe.add(W-1)
    return all(i in safe or not exists[i] or any(
        y>b[i] for y in (a[i], c[i], b[i-1] if i else -1, b[i+1] if i+1<W else -1)
    ) for i in range(W))


def exact_word_counts(points, free):
    ids={v:i for i,v in enumerate(points)}
    N=len(points)
    full=(1<<N)-1
    allowed=[p in free for p in points]
    masks=[sum(1<<ids[w] for w in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if w in ids) for x,y in points]
    dp=[0]*(full+1)
    dp[0]=1
    for s,z in enumerate(dp):
        if z==0:
            continue
        new=full^s
        while new:
            bit=new&-new
            new^=bit
            j=bit.bit_length()-1
            if allowed[j] or masks[j]&s:
                dp[s|bit]+=z
    return dp[full]


def case_graph(W,n,b):
    B={(0,0)}
    if b==2:
        B.add((0,1))
    A={(x,y) for x in range(W,n) for y in range(n)}|{(x,0) for x in range(2,W)}
    K=A|{(1,0)}
    S=K|B
    entire={(x,y) for x in range(n) for y in range(n)}
    H=entire-S
    def boundary(P):
        return {v for v in H if any(q in P for q in ((v[0]-1,v[1]),(v[0]+1,v[1]),(v[0],v[1]-1),(v[0],v[1]+1)))}
    return sorted(H), boundary(K), boundary(S)


def quantitative_c(W):
    c=Fraction(1,factorial(W))
    eps=c*c/(2*W)
    return c*eps**W*Fraction(1,2**(5*W+1)*factorial(4*W))


def exact_finite():
    recorded=[]
    for W, ns in ((2,range(3,11)),(3,range(4,8)),(4,range(5,7))):
        for n in ns:
            for b in (1,2):
                H, free_target,free_source=case_graph(W,n,b)
                same=exact_word_counts(H,free_target)
                enlarged=exact_word_counts(H,free_source)
                test=Fraction(same,enlarged)
                assert 0<test<=1
                if n>=10:
                    assert test>=quantitative_c(W)
                recorded.append((W,n,b,str(test)))
    for k in recorded:
        print("width, n, b, exact quotient:",*k)
    return len(recorded)


def audit_append():
    checked=0
    for w in (2,3,4):
        probs=[]
        for k in range(2, min(20//w, 6)+1):
            V=[(x,j) for j in range(k) for x in range(w)]
            free={(x,0) for x in range(w)}|{(w-1,j) for j in range(k)}
            probs.append(Fraction(exact_word_counts(V,free),factorial(w*k)))
        for a,b in zip(probs,probs[1:]):
            assert b*factorial(w)>=a
            checked+=1
    return checked


def audit_bridges():
    rng=Random(20261118)
    total=0
    for W in range(2,15):
        for _ in range(100):
            a,b,d,e=[[Fraction(t,10**8) for t in rng.sample(range(1,10**8),W)] for _ in range(4)]
            h=max(b+d)
            x=[h+(1-h)*Fraction(t,W+1) for t in range(1,W+1)]
            assert last_row_factor(a,b,x)
            assert last_row_factor(b,x,d)
            assert last_row_factor(x,d,e)
            total+=1
    return total


def audit_row_products():
    checked=0
    rng=Random(20261119)
    for W in (2,3,4,7):
        for n in (10,15):
            for b in (1,2):
                H,free_new,free_old=case_graph(W,n,b)
                Hset=set(H)
                for trials in range(120):
                    # Mix admissible horizontal-row starts, local perturbations and random words.
                    if trials<60:
                        ranked=sorted(H,key=lambda v:v[0])
                        if trials>0:
                            p,q=rng.sample(range(len(ranked)),2)
                            ranked[p],ranked[q]=ranked[q],ranked[p]
                    else:
                        ranked=rng.sample(H,len(H))
                    label={v:i for i,v in enumerate(ranked)}
                    rows=[[label.get((x,y),-1) for x in range(W)] for y in range(n+1)]
                    for free in (free_new,free_old):
                        direct=all(p in free or any(label.get(q,-1)>label[p] for q in ((p[0]-1,p[1]),(p[0]+1,p[1]),(p[0],p[1]-1),(p[0],p[1]+1))) for p in H)
                        prefix=all(last_row_factor(rows[y-1],rows[y],rows[y+1],False,
                                                   [(x,y) in Hset for x in range(W)],
                                                   [x for x in range(W) if (x,y) in free])
                                   for y in range(1,4))
                        tail=all(last_row_factor(rows[y-1],rows[y],rows[y+1])
                                 for y in range(4,n))
                        assert direct==(prefix and tail)
                        checked+=1
    return checked


if __name__=="__main__":
    print("finite deletion count cases:",exact_finite())
    print("one-free-row probabilities checked under appending:",audit_append())
    print("bridging row constructions:",audit_bridges())
    print("full direct no peaks versus row-product equality:",audit_row_products())
