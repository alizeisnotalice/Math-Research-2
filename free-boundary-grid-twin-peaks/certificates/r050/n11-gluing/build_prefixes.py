"""Plan symmetry-reduced shortest and short connected-prefix constraints."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def adjacency(w,h):
    return [sum(1<<u for u in range(w*h) if abs(u%w-v%w)+abs(u//w-v//w)==1) for v in range(w*h)]


def bits(s):
    while s:
        b=s & -s
        yield b.bit_length()-1
        s^=b


def neighbors(adj,s):
    out=0
    for v in bits(s):
        out|=adj[v]
    return out


def canonical(w,h,T):
    points=[(v%w,v//w) for v in bits(T)]
    if w>h:
        w,h=h,w
        points=[(y,x) for x,y in points]
    families=(points,[(y,x) for x,y in points]) if w==h else (points,)
    val=min(sum(1<<((h-1-y if f_y else y)*w+(w-1-x if f_x else x)) for x,y in pool)
            for pool in families for f_x in (0,1) for f_y in (0,1))
    return w,h,val


@lru_cache(maxsize=None)
def count_independent(w,h,avail):
    last={0:[1]}
    for y in range(h):
        after={}
        row=(avail>>(w*y))&((1<<w)-1)
        s=row
        while True:
            if not s & (s<<1):
                shift=s.bit_count()
                hist=[0]*(max(map(len,last.values()))+shift)
                for t,prev in last.items():
                    if s & t:
                        continue
                    for k,c in enumerate(prev):
                        hist[k+shift]+=c
                after[s]=hist
            if s==0:
                break
            s=(s-1)&row
        last=after
    return sum(map(sum,last.values())),sum(k*c for v in last.values() for k,c in enumerate(v))


def small_words(w,h,u,K,adj):
    bins=defaultdict(Counter)
    current={1<<u:1}
    seam=((1<<w)-1)<<(w*(h-1))
    if (1<<u)&seam:
        bins[u%w,1][1<<u]+=1
    for j in range(2,K+1):
        nxt=defaultdict(int)
        for S,c in current.items():
            front=neighbors(adj,S)&~S
            for v in bits(front):
                T=S|(1<<v)
                if j<K:
                    nxt[T]+=c
                if (1<<v)&seam:
                    bins[v%w,j][T]+=c
        current=nxt
    return bins


def minimal_words(w,h,u,target):
    """Every path decreases distance to target by exactly one at each step."""
    x0,y0=u%w,u//w
    xt,yt=target%w,target//w
    def visit(x,y,T):
        if (x,y)==(xt,yt):
            yield T
        if x!=xt:
            x1=x+(1 if xt>x else -1)
            yield from visit(x1,y,T|(1<<(y*w+x1)))
        if y!=yt:
            y1=y+(1 if yt>y else -1)
            yield from visit(x,y1,T|(1<<(y1*w+x)))
    yield from visit(x0,y0,1<<u)


def make_plan(shapes,cap):
    contexts={}
    data=[]
    skipped=0
    for w,h,K in shapes:
        adj=adjacency(w,h)
        for y in range(h):
            for x in range(w-1):
                u=y*w+x
                bins=small_words(w,h,u,K,adj)
                for c in range(w):
                    min_j=abs(c-x)+(h-1-y)+1
                    if min_j>K:
                        bins[c,min_j].update(minimal_words(w,h,u,(h-1)*w+c))
                for (c,j),proposal in sorted(bins.items()):
                    row=Counter()
                    for T,num in proposal.items():
                        ww,hh,H=canonical(w,h,T)
                        if H.bit_count()==1:
                            # A singleton's completion is the archived rooted total.
                            row[H]+=num
                            continue
                        key=ww,hh,H
                        if key not in contexts:
                            F=((1<<(ww*hh))-1)&~(H|neighbors(adjacency(ww,hh),H))
                            states,links=count_independent(ww,hh,F)
                            contexts[key]=(states,links)
                        if contexts[key][0] <= cap:
                            row[H]+=num
                        else:
                            skipped+=num
                    if row:
                        data.append({'w':w,'h':h,'root':u,'connector_x':c,'local_rank':j,
                                     'prefixes':[[str(b),f] for b,f in sorted(row.items())]})
    tasks=[{'w':w,'h':h,'seen':str(T),'states':z,'predecessor_lookups':l}
           for (w,h,T),(z,l) in sorted(contexts.items()) if z<=cap]
    return {'shapes':shapes,'max_per_context':cap,'unique_observed_contexts':len(contexts),
            'pruned_ordered_prefixes':skipped,'rows':data,'tasks':tasks}


def main():
    ap=ArgumentParser()
    ap.add_argument('output',type=Path)
    ap.add_argument('--rect-cutoff',type=int,default=4)
    ap.add_argument('--square-cutoff',type=int,default=3)
    ap.add_argument('--max-per-context',type=int,default=1_200_000)
    ap.add_argument('--small-gate',action='store_true')
    a=ap.parse_args()
    shapes=[(3,3,4)] if a.small_gate else [(5,6,a.rect_cutoff),(6,5,a.rect_cutoff),(6,6,a.square_cutoff)]
    data=make_plan(shapes,a.max_per_context)
    with a.output.open('x') as f:
        json.dump(data,f,indent=2);f.write('\n')
    for w,h in sorted({(r['w'],r['h']) for r in data['tasks']}):
        T=[r for r in data['tasks'] if (r['w'],r['h'])==(w,h)]
        print(w,h,'contexts',len(T),'states total',sum(r['states'] for r in T),
              'max',max(r['states'] for r in T), 'lookups',sum(r['predecessor_lookups'] for r in T))
    print('excluded',data['pruned_ordered_prefixes'])


if __name__=='__main__':
    main()
