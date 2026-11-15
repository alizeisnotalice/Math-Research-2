#!/usr/bin/env python3
"""Exact Q(sqrt(2)) verification of the finite-column mixed-halo example.
Only this finite example is computed; no claim of computing the general hierarchy.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
from itertools import product, combinations
from pathlib import Path
import hashlib, json, os, platform, sys

@total_ordering
@dataclass(frozen=True)
class Q2:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))
    @staticmethod
    def coerce(x):
        return x if isinstance(x,Q2) else Q2(F(x))
    def __add__(self,x):
        x=self.coerce(x); return Q2(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Q2(-self.a,-self.b)
    def __sub__(self,x): return self+-self.coerce(x)
    def __rsub__(self,x): return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x); return Q2(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.coerce(x)
        den=x.a*x.a-2*x.b*x.b
        if den==0: raise ZeroDivisionError
        return self*Q2(x.a/den,-x.b/den)
    def sign(self):
        a,b=self.a,self.b
        if not b: return (a>0)-(a<0)
        if not a: return (b>0)-(b<0)
        if a>0 and b>0: return 1
        if a<0 and b<0: return -1
        s=(a*a>2*b*b)-(a*a<2*b*b)
        return s if a>0 else -s
    def __lt__(self,x): return (self-self.coerce(x)).sign()<0
    def __eq__(self,x):
        x=self.coerce(x); return self.a==x.a and self.b==x.b
    def record(self):
        return {'rational':str(self.a),'sqrt2_coefficient':str(self.b)}

ZERO=Q2(); ONE=Q2(1); HALF=Q2(F(1,2)); G=Q2(3,-2); OPT=Q2(4,-2)
assert G>0

def rectangle_union(rects):
    """Exact slab union; clipped endpoints stay in Q(sqrt(2))."""
    xs=sorted({x for r in rects for x in r[:2]})
    area=ZERO
    for left,right in zip(xs,xs[1:]):
        middle=(left+right)/2
        intervals=sorted((lo,hi) for x0,x1,lo,hi in rects if x0<middle<x1)
        total=ZERO
        if intervals:
            lo,hi=intervals[0]
            for a,b in intervals[1:]:
                if a<=hi: hi=max(hi,b)
                else: total+=hi-lo; lo,hi=a,b
            total+=hi-lo
        area+=(right-left)*total
    return area

def full_local_area(window):
    # Receiver [0,1] x [0,4], cap R=1. All possible vertical lifts
    # that can meet it occur at y=0,4; other lifts are distance >=4.
    pts=[(j,y) for j,c in zip(range(-2,3),window) if c for y in (0,4)]
    rects=[]; examined=0
    # Enumerate every nonempty subset, rejecting span >=2R before roots.
    for count in range(1,len(pts)+1):
        for subset in combinations(pts,count):
            examined+=1
            x0=min(p[0] for p in subset); x1=max(p[0] for p in subset)
            y0=min(p[1] for p in subset); y1=max(p[1] for p in subset)
            if x1-x0>=2 or y1-y0>=2: continue
            assert count in (1,2), (window,subset)
            radius=HALF if count==1 else Q2(0,F(1,2))
            left=max(ZERO,Q2(x1)-radius); right=min(ONE,Q2(x0)+radius)
            bottom=max(ZERO,Q2(y1)-radius); top=min(Q2(4),Q2(y0)+radius)
            if left<right and bottom<top: rects.append((left,right,bottom,top))
    return rectangle_union(rects),examined,len(rects)

def main():
    assert os.path.ismount('/Volumes/data'), 'data volume must be mounted'
    windows=list(product((0,1),repeat=5))
    area_by_window={}; geometry=[]
    for w in windows:
        actual,count,nrect=full_local_area(w)
        expected=Q2(F(w[2]+w[3],2))+G*(w[2]*w[3])
        assert actual==expected,(w,actual,expected)
        area_by_window[w]=actual
        geometry.append({'window':w,'area':actual.record(),'subsets_checked':count,'nonempty_rectangles':nrect})
    vertices=list(product((0,1),repeat=4))
    ids={v:i for i,v in enumerate(vertices)}
    edges={i:[] for i in range(len(vertices))}
    potential=[]
    for w in windows:
        u=ids[w[:-1]]; v=ids[w[1:]]
        a=area_by_window[w]; b=w[2]
        p_u=Q2(F(w[2],2)); p_v=Q2(F(w[3],2))
        slack=p_v-p_u-(a-OPT*b)
        assert slack==G*w[2]*(1-w[3]) and slack>=0
        edges[u].append((v,a,b))
        potential.append({'edge':[u,v],'slack':slack.record()})
    # Enumerate all directed simple cycles once, starting at their minimum ID.
    cycles=[]
    def dfs(start,cur,path,seen,A,M):
        for nxt,a,b in edges[cur]:
            if nxt==start:
                cycles.append((tuple(path),A+a,M+b))
            elif nxt>start and nxt not in seen:
                dfs(start,nxt,path+[nxt],seen|{nxt},A+a,M+b)
    for start in range(len(vertices)):
        dfs(start,start,[start],{start},ZERO,0)
    positive=[]; zero=[]
    for path,A,M in cycles:
        if M:
            ratio=A/M
            assert ratio<=OPT
            positive.append((path,A,M,ratio))
        else:
            assert A==ZERO and path==(ids[(0,0,0,0)],)
            zero.append(path)
    best=max(x[3] for x in positive)
    assert best==OPT
    maximizers=[path for path,A,M,ratio in positive if ratio==OPT]
    assert (ids[(1,1,1,1)],) in maximizers
    # Symbolic root identity at z=OPT and strict positive cost slopes.
    assert max(A-OPT*M for _,A,M,_ in positive)==ZERO
    assert min(M for _,_,M,_ in positive)>0
    source=Path(__file__).resolve()
    result={
      'status':'PASS','scope':'fixed 16-vertex binary one-row finite graph only',
      'python_version':sys.version,'platform':platform.platform(),
      'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'arithmetic':'fractions.Fraction exact Q(sqrt(2)); no floating comparison',
      'seed':None,'tolerance':None,
      'parameters':{'h':'1','vertical_period':'4','radius_cap':'1','k':2,'column_alphabet':[0,1]},
      'vertices':len(vertices),'edges':len(windows),'windows_checked':len(geometry),
      'source_subsets_checked':sum(g['subsets_checked'] for g in geometry),
      'simple_cycles':len(cycles),'positive_mass_simple_cycles':len(positive),'zero_mass_cycles':len(zero),
      'optimum':best.record(),'maximizing_cycles':maximizers,
      'geometry':geometry,'potential_edge_checks':potential,
      'cycles':[{'vertices':path,'area':A.record(),'mass':M,'ratio':ratio.record()} for path,A,M,ratio in positive],
      'limitations':'This run does not compute higher hierarchy levels or prove a general C2 upper bound.'}
    out=source.parent/'results/recursive-boundary/cycle-example.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('status','python_version','vertices','edges','windows_checked','source_subsets_checked','simple_cycles','positive_mass_simple_cycles','zero_mass_cycles','optimum','maximizing_cycles')},indent=2))
    print('output='+str(out))

if __name__=='__main__': main()
