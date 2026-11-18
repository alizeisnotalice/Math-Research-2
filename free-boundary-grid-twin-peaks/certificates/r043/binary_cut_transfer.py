"""Exact four-state count of the binary central-column cut multiplicity.

All labels are represented by a descending word; smaller rank is higher.
The transfer identity is for odd n>=7, outer quarter safe roots.
"""
import itertools,json,random,sys
from pathlib import Path

def geometry(n):
    nb=[]
    for v in range(n*n):
        x,y=v%n,v//n
        nb.append([Y*n+X for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=X<n and 0<=Y<n])
    return nb

def ranks(word):
    pos=[0]*len(word)
    assert sorted(word)==list(range(len(word)))
    for i,v in enumerate(word):pos[v]=i
    return pos

def transfer(n,word):
    assert n>=7 and n%2
    m=n//2;d=(n+3)//4;nb=geometry(n);pos=ranks(word)
    A={v for v in range(n*n) if v%n<m};B={v for v in range(n*n) if v%n>m}
    lp={v for v in A if all(pos[v]<pos[u] for u in nb[v] if u in A)}
    rp={v for v in B if all(pos[v]<pos[u] for u in nb[v] if u in B)}
    li=[v for v in lp if v%n!=m-1];ri=[v for v in rp if v%n!=m+1]
    if len(li)!=1 or len(ri)!=1 or li[0]%n>=d or ri[0]%n<n-d:return 0
    force=[-1]*n
    for y in range(n):
        c=y*n+m
        for p,s,peaks in ((c-1,0,lp),(c+1,1,rp)):
            if p in peaks:
                if pos[c]>=pos[p] or force[y] not in (-1,s):return 0
                force[y]=s
    high=[[pos[y*n+m-1]<pos[y*n+m] for y in range(n)],
          [pos[y*n+m+1]<pos[y*n+m] for y in range(n)]]
    def good(y,prev,s,nxt):
        return (high[s][y] or (y>0 and prev==s and pos[(y-1)*n+m]<pos[y*n+m])
                or (y+1<n and nxt==s and pos[(y+1)*n+m]<pos[y*n+m]))
    dp={(2,s):1 for s in (0,1) if force[0] in (-1,s)}
    for y in range(n-1):
        new={}
        for (prev,s),count in dp.items():
            for t in (0,1):
                if force[y+1] in (-1,t) and good(y,prev,s,t):new[(s,t)]=new.get((s,t),0)+count
        dp=new
    return sum(count for (prev,s),count in dp.items() if good(n-1,prev,s,2))

def brute(n,word):
    nb=geometry(n);pos=ranks(word);m=n//2;d=(n+3)//4;ans=0
    for bits in itertools.product((0,1),repeat=n):
        H={v for v in range(n*n) if v%n<m+bits[v//n]};K=set(range(n*n))-H
        peaks=[]
        for region in (H,K):peaks.append([v for v in region if all(pos[v]<pos[u] for u in nb[v] if u in region)])
        if len(peaks[0])!=1 or len(peaks[1])!=1:continue
        a,b=peaks[0][0],peaks[1][0]
        if a%n>=d or b%n<n-d:continue
        if any(u not in H for u in nb[a]) or any(u not in K for u in nb[b]):continue
        ans+=1
    return ans

def full_peaks(n,word):
    nb=geometry(n);pos=ranks(word)
    return [v for v in range(n*n) if all(pos[v]<pos[u] for u in nb[v])]

def maximum_word(n):
    assert n>=5 and n%2;m=n//2
    return ([y*n+x for y in range(n) for x in range(m)]
            +[y*n+x for y in range(n) for x in range(n-1,m,-1)]
            +[y*n+m for y in range(n)])

def no_hit_word(n):
    assert n>=5 and n%2;m=n//2;nb=geometry(n)
    word=[0,n-1]+list(range(1,m))+[y*n+m for y in range(3)]+[2*n+m-1,2*n+m+1]
    assert len(set(word))==len(word);seen=set(word)
    for v in word[2:]:assert any(u in word[:word.index(v)] for u in nb[v])
    front={u for v in seen for u in nb[v]}-seen
    while front:
        v=min(front);front.remove(v);word.append(v);seen.add(v);front|=set(nb[v])-seen
    assert len(word)==n*n
    return word

def growth_word(region,root,nb,rng):
    S={root};word=[root];front=(set(nb[root])&region)-S
    while front:
        v=rng.choice(sorted(front));front.remove(v);word.append(v);S.add(v);front|=(set(nb[v])&region)-S
    assert S==region
    return word

def constructed(n,rng):
    nb=geometry(n);m=n//2;d=(n+3)//4;bits=[rng.randrange(2) for _ in range(n)]
    L={v for v in range(n*n) if v%n<m+bits[v//n]};R=set(range(n*n))-L
    a=rng.choice([v for v in sorted(L) if v%n<d and all(u in L for u in nb[v])])
    b=rng.choice([v for v in sorted(R) if v%n>=n-d and all(u in R for u in nb[v])])
    lw=growth_word(L,a,nb,rng);rw=growth_word(R,b,nb,rng);word=[];i=j=0
    while i<len(lw) or j<len(rw):
        if rng.randrange(len(lw)+len(rw)-i-j)<len(lw)-i:word.append(lw[i]);i+=1
        else:word.append(rw[j]);j+=1
    return word

def audit(out):
    out.mkdir(exist_ok=False);rng=random.Random(2026111779);rows=[];positive=0
    for n,cases in ((7,40),(9,30),(11,8)):
        observed=[]
        for k in range(cases):
            word=constructed(n,rng);peaks=full_peaks(n,word);assert len(peaks)==2
            x=transfer(n,word);y=brute(n,word);assert x==y and x>0
            positive+=1;observed.append(x)
        rows.append(dict(n=n,constructed_cases=cases,multiplicities=observed))
    # Arbitrary words test the necessary zero conditions as well.
    for n in (7,9):
        word=list(range(n*n))
        for _ in range(10):rng.shuffle(word);assert transfer(n,word)==brute(n,word)
    maxima=[]
    for n in (5,7,9,11,15,31,65):
        word=maximum_word(n);peaks=full_peaks(n,word);assert peaks==[0,n-1]
        value=brute(n,word) if n==5 else transfer(n,word)
        assert value==1<<n
        if n in (7,9,11):assert brute(n,word)==value
        maxima.append(dict(n=n,X=value,distance=n-1))
    zero_examples=[]
    for n in (5,7,9,11,31,65):
        word=no_hit_word(n);assert full_peaks(n,word)==[0,n-1]
        value=brute(n,word) if n==5 else transfer(n,word);assert value==0
        if n in (7,9,11):assert brute(n,word)==0
        pos=ranks(word);nb=geometry(n);m=n//2
        for v,side in ((2*n+m-1,0),(2*n+m+1,1)):
            assert all(pos[v]<pos[u] for u in nb[v] if (u%n<m if side==0 else u%n>m))
        zero_examples.append(dict(n=n,X=0,distance=n-1))
    result=dict(scope='Exact identity checks only; constructed draws are not conditional probability samples',seed=2026111779,positive_cases=positive,arbitrary_words=20,rows=rows,maxima=maxima,zero_examples=zero_examples)
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(positive_cases=positive,arbitrary_words=20,maxima=maxima,zero_examples=zero_examples)))

if __name__=='__main__':audit(Path(sys.argv[1]))
