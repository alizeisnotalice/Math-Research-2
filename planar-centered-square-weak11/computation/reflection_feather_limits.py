#!/usr/bin/env python3
"""Exact, bounded verification of reflection degeneration examples.

All arithmetic lies in Q(sqrt(2)); subset boxes with larger radicals are
proved empty by comparing squared box spans with subset mass before roots.
No floating point is used for any assertion. This is not a search.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
from itertools import combinations
from math import isqrt
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
    def of(x):
        return x if isinstance(x, Q2) else Q2(x)
    def __add__(self, x):
        x=Q2.of(x); return Q2(self.a+x.a, self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Q2(-self.a,-self.b)
    def __sub__(self, x): return self+-Q2.of(x)
    def __rsub__(self, x): return Q2.of(x)+-self
    def __mul__(self, x):
        x=Q2.of(x)
        return Q2(self.a*x.a+2*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self, x):
        x=Q2.of(x); d=x.a*x.a-2*x.b*x.b
        assert d
        return self*Q2(x.a/d,-x.b/d)
    def sign(self):
        a,b=self.a,self.b
        if not a: return (b>0)-(b<0)
        if not b: return (a>0)-(a<0)
        if a>0 and b>0: return 1
        if a<0 and b<0: return -1
        d=a*a-2*b*b
        # sqrt(2) is irrational, so this is nonzero for rational a,b != 0.
        assert d
        return ((d>0)-(d<0)) * ((a>0)-(a<0))
    def __eq__(self, x):
        x=Q2.of(x); return self.a==x.a and self.b==x.b
    def __lt__(self,x): return (self-Q2.of(x)).sign()<0
    def __hash__(self): return hash((self.a,self.b))
    def as_json(self): return {'rational':str(self.a),'sqrt2_coefficient':str(self.b)}


def sqrt_fraction(x):
    """Return an exact Q2 root, or fail if proof did not remove a new radical."""
    x=F(x)
    for factor in [1,2]:
        y=x/factor
        n,d=isqrt(y.numerator),isqrt(y.denominator)
        if n*n==y.numerator and d*d==y.denominator:
            return Q2(F(n,d),0) if factor==1 else Q2(0,F(n,d))
    raise AssertionError(('unexpected nonempty radical',x))


def complete_rectangles(atoms):
    out=[]; empty=0
    for mask in range(1,1<<len(atoms)):
        pts=[atoms[i] for i in range(len(atoms)) if mask>>i&1]
        mass=sum((q[2] for q in pts),F(0))
        xmin=min(q[0] for q in pts); xmax=max(q[0] for q in pts)
        ymin=min(q[1] for q in pts); ymax=max(q[1] for q in pts)
        # Strict threshold 1, common radius. Equality gives an empty open side.
        if (xmax-xmin)*(xmax-xmin)>=mass or (ymax-ymin)*(ymax-ymin)>=mass:
            empty+=1; continue
        r=sqrt_fraction(mass/4)
        box=(Q2.of(xmax)-r,Q2.of(xmin)+r,Q2.of(ymax)-r,Q2.of(ymin)+r)
        assert box[0]<box[1] and box[2]<box[3]
        out.append(box)
    return out,empty


def union_area(rectangles):
    xs=sorted(set(x for q in rectangles for x in q[:2]))
    area=Q2()
    for l,r in zip(xs,xs[1:]):
        mid=(l+r)/2
        intervals=sorted((q[2],q[3]) for q in rectangles if q[0]<mid<q[1])
        if not intervals: continue
        lo,hi=intervals[0]; length=Q2()
        for a,b in intervals[1:]:
            if a>hi: length+=hi-lo; lo,hi=a,b
            elif b>hi: hi=b
        length+=hi-lo
        area+=(r-l)*length
    return area


def main():
    assert os.path.ismount('/Volumes/data')
    sq2=Q2(0,1); k=sq2-1; g=Q2(F(3,2),-1)
    assert k*k==Q2(3,-2) and k*k/2==g and g>0
    core=[(F(-1,2),F(0),F(1)),(F(1,2),F(0),F(1))]
    cr,_=complete_rectangles(core); core_area=union_area(cr)
    assert core_area==Q2(5,-2)
    records=[]
    for epsilon in [F(1),F(1,4),F(1,16),F(1,256),F(1,65536)]:
        root=sqrt_fraction(epsilon)
        assert not root.b
        s=root.a; D=F(4); total=2+epsilon
        # A=1/2. This rational squared comparison implies all mixed boxes empty.
        assert (D-F(1,2))**2>2*total
        atoms=core+[(D,F(0),epsilon)]
        c=D+s/2
        reflected=[(2*c-x,y,m) for x,y,m in atoms]
        old_rect,old_empty=complete_rectangles(atoms)
        new_rect,new_empty=complete_rectangles(atoms+reflected)
        old_area=union_area(old_rect); new_area=union_area(new_rect)
        assert max(q[1] for q in old_rect)==c
        assert old_area==core_area+epsilon
        assert new_area==2*core_area+(2+k*k)*epsilon
        excess=new_area-2*old_area
        gain=new_area/(2*total)-old_area/total
        assert excess==k*k*epsilon and gain==g*epsilon/total
        records.append({'epsilon':str(epsilon),'old_atoms':len(atoms),'new_atoms':len(atoms+reflected),
            'old_nonempty_subsets':len(old_rect),'new_nonempty_subsets':len(new_rect),
            'old_empty_subsets':old_empty,'new_empty_subsets':new_empty,
            'seam':str(c),'old_area':old_area.as_json(),'new_area':new_area.as_json(),
            'actual_excess_area':excess.as_json(),'actual_ratio_gain':gain.as_json()})
    # The old delta-based certificate vanishes although a fixed true gain remains.
    delta_records=[]
    for t in [F(1,16),F(1,256),F(1,65536)]:
        eta=k/2-t
        assert eta>0
        a=[(Q2(0),Q2(0),F(1)),(eta,Q2(0),F(1))]
        old,_=complete_rectangles(a)
        c=sq2/2
        ends=sorted(set(q[1] for q in old))
        assert max(ends)==c and min(c-e for e in ends if e<c)==t
        # All old rectangles have |y|<sqrt(2)/2.
        assert max(q[3] for q in old)==sq2/2
        # All four reflected sources have x-span sqrt(2), mass 4, y-span 0.
        both=a+[(2*c-x,y,m) for x,y,m in a]
        assert max(p[0] for p in both)-min(p[0] for p in both)==sq2
        excess_lower=(2-sq2)*(2-sq2)
        assert excess_lower/4==g
        delta_records.append({'t':str(t),'eta':eta.as_json(),'delta':str(t),
                              'true_ratio_gain_lower':g.as_json()})
    output={'status':'PASS','scope':'Exact finite examples and scalar identities; general-parameter proofs are in the paired markdown.',
        'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
        'arithmetic':'fractions.Fraction and exact Q(sqrt(2)) sign comparisons',
        'random_seed':None,'floating_tolerance':None,'lambda':'1','single_process':True,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256':{'proofs/reflection-nonattainment.md':'384968620d3132cddf77252ab864b42b178531b547e2a0abbeb902e78e00c745'},
        'feather_cases':records,'delta_cases':delta_records}
    path=Path(__file__).parent/'results/reflection-feather-limits/exact-check.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','feather_cases':len(records),'delta_cases':len(delta_records),
                      'max_atom_count':6,'max_subsets':63,'output':str(path)},ensure_ascii=False))

if __name__=='__main__': main()
