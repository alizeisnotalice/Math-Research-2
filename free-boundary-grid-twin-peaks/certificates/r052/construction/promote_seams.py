"""Reversible four-quadrant promotion map for the all-n far-peak bound."""
from argparse import ArgumentParser
from itertools import combinations
from math import factorial
from pathlib import Path
import json
import random
import sys


class Square:
    def __init__(self,n):
        if n<4:raise ValueError('need both block widths at least two')
        self.n=n
        m=n//2
        self.vertices=tuple(range(n*n))
        self.full=set(self.vertices)
        self.adj=[]
        for v in self.vertices:
            x,y=v%n,v//n
            self.adj.append({Y*n+X for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=X<n and 0<=Y<n})
        self.blocks=[]
        self.seams=[]
        self.roots=[]
        self.block_of=[0]*(n*n)
        self.to_seam={}
        intervals=((0,m),(m,n))
        for j,(y0,y1) in enumerate(intervals):
            for i,(x0,x1) in enumerate(intervals):
                block={y*n+x for x in range(x0,x1) for y in range(y0,y1)}
                assert len(block)>=4
                self.blocks.append(block)
                for v in block:self.block_of[v]=2*j+i
                x_end=x1-1 if i==0 else x0
                y_end=y1-1 if j==0 else y0
                root_x=x0 if i==0 else x1-1
                self.roots.append(y_end*n+root_x)
                seam={v for v in block if self.adj[v]-block}
                assert seam=={v for v in block if v%n==x_end or v//n==y_end}
                self.seams.append(seam)
                for v in block-seam:
                    x,y=v%n,v//n
                    if abs(x-x_end)<=abs(y-y_end):
                        path={n*y+X for X in range(min(x,x_end),max(x,x_end)+1)}
                    else:
                        path={n*Y+x for Y in range(min(y,y_end),max(y,y_end)+1)}
                    self.to_seam[v]=path-seam
        assert len(self.to_seam)==(n-2)**2
        assert sum(map(len,self.seams))==4*n-4

    def support(self,protected_peaks):
        sets=[set(s) for s in self.seams]
        for v in protected_peaks:
            sets[self.block_of[v]].update(self.to_seam[v])
        return sets

    def global_peaks(self,word):
        loc={v:i for i,v in enumerate(word)}
        return tuple(v for v in self.vertices if all(loc[v]<loc[u] for u in self.adj[v]))

    def ordered_support(self,root,S):
        todo=S-{root}
        order=[root]
        frontier=set(self.adj[root])
        while todo:
            v=min(todo&frontier)
            todo.remove(v)
            order.append(v)
            frontier.update(self.adj[v])
        return order

    def promote(self,word):
        peaks=self.global_peaks(word)
        assert len(peaks)==2
        original=[[v for v in word if self.block_of[v]==j] for j in range(4)]
        loc={v:i for i,v in enumerate(word)}
        protected_local=tuple(v for v in self.to_seam if all(loc[v]<loc[u] for u in self.adj[v]&self.blocks[self.block_of[v]]))
        assert set(protected_local)==(set(peaks)&self.to_seam.keys())
        sets=self.support(protected_local)
        assert sum(map(len,sets))<=4*self.n+2*((self.n+1)//2)-6
        words=[]
        slots=[]
        for j,S in enumerate(sets):
            slots.append(tuple(original[j].index(v) for v in sorted(S)))
            words.append(iter(self.ordered_support(self.roots[j],S)+[v for v in original[j] if v not in S]))
        dest=tuple(next(words[self.block_of[v]]) for v in word)
        return dest,(protected_local,slots)

    def undo(self,dest,record):
        protected,slots=record
        sets=self.support(protected)
        full_order=[]
        for j,S in enumerate(sets):
            leftover=iter(v for v in dest if self.block_of[v]==j and v not in S)
            block_order=[None]*len(self.blocks[j])
            for v,k in zip(sorted(S),slots[j],strict=True):
                assert block_order[k] is None
                block_order[k]=v
            full_order.append(iter(tuple(next(leftover) if v is None else v for v in block_order)))
        return tuple(next(full_order[self.block_of[v]]) for v in dest)

    def capacity(self):
        ans=0
        P=tuple(self.to_seam)
        for k in range(3):
            for protected in combinations(P,k):
                factor=1
                for H,S in zip(self.blocks,self.support(protected),strict=True):
                    factor*=factorial(len(H))//factorial(len(H)-len(S))
                ans+=factor
        return ans


def grow_k2(grid,rg,planned_late):
    start=rg.choice(grid.vertices)
    seen={start}
    word=[start]
    front=set(grid.adj[start])
    while len(word)<planned_late and grid.full-(seen|front):
        v=rg.choice(sorted(front))
        seen.add(v);word.append(v);front=(front|grid.adj[v])-seen
        if not grid.full-(seen|front):
            # Undo this last ordinary growth step to keep a possible birth.
            word.pop();seen.remove(v)
            front=set().union(*(grid.adj[w] for w in seen))-seen
            break
    late=len(word)+1
    newborn=rg.choice(sorted(grid.full-(seen|front)))
    seen.add(newborn);word.append(newborn);front=(front|grid.adj[newborn])-seen
    while front:
        v=rg.choice(sorted(front))
        seen.add(v);word.append(v);front=(front|grid.adj[v])-seen
    assert len(word)==len(grid.vertices) and len(grid.global_peaks(word))==2
    return tuple(word),late


def main():
    parser=ArgumentParser()
    parser.add_argument('output',type=Path)
    a=parser.parse_args()
    seed=2026111852
    rng=random.Random(seed)
    entries=[]
    for n in (4,5,6,7,12):
        graph=Square(n)
        res=graph.capacity()
        peak_patterns=set()
        smallest=0
        largest=0
        for tries in range(400):
            word,late=grow_k2(graph,rng,1+tries%n)
            dest,record=graph.promote(word)
            back=graph.undo(dest,record)
            assert back==word and len(set(dest))==n*n
            peaks=graph.global_peaks(dest)
            assert len(peaks)==2 and abs(peaks[0]%n-peaks[1]%n)==n-1
            assert {graph.block_of[v] for v in peaks} in ({0,1},{0,3},{2,1},{2,3})
            smallest=min(smallest or late,late);largest=max(largest,late)
            peak_patterns.add(record[0])
        entries.append({'n':n,'sampled_valid_inputs':400,'earliest_input_second_peak_rank':smallest,
                        'latest_input_second_peak_rank':largest,'protected_peak_patterns_seen':len(peak_patterns),
                        'exact_uniform_fiber_upper_from_slot_capacity':str(res)})
    out={'scope':'Deterministic proof separate; random growth here only tests map/inverse, not conditional probabilities',
         'python':sys.version,'seed':seed,'positive_inverse_checks':sum(r['sampled_valid_inputs'] for r in entries),'entries':entries}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(out['positive_inverse_checks'],'map and inverse gates')

if __name__=='__main__':main()
