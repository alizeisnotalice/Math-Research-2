#!/usr/bin/env python3
"""Complete strict-halo unions in Q(sqrt(2),sqrt(3),sqrt(5)); directed rational comparisons."""
from fractions import Fraction as F
from functools import lru_cache,total_ordering
from itertools import combinations
from math import isqrt
from pathlib import Path
import hashlib,json,os,sys
assert os.path.ismount('/Volumes/data')
PRIMES=(2,3,5)
RADS=tuple(__import__('functools').reduce(lambda a,j:a*PRIMES[j],(j for j in range(3) if m>>j&1),1) for m in range(8))
SIGN_BITS=[]
@lru_cache(None)
def sign(v):
    if not any(v):return 0
    for bits in (16,32,64,128,256,512,1024):
        den=1<<bits;lo=F(0);hi=F(0)
        for a,rad in zip(v,RADS):
            if not a:continue
            if rad==1:l=u=F(1)
            else:
                n=isqrt(rad*den*den);l=F(n,den);u=F(n+1,den)
            lo+=a*(l if a>0 else u);hi+=a*(u if a>0 else l)
        if lo>0:SIGN_BITS.append(bits);return 1
        if hi<0:SIGN_BITS.append(bits);return -1
    raise ArithmeticError('comparison could not be certified')
@total_ordering
class R:
    def __init__(self,v=0):
        if isinstance(v,R):self.v=v.v
        elif isinstance(v,(tuple,list)):self.v=tuple(F(x) for x in v)
        else:self.v=(F(v),)+(F(0),)*7
    def __add__(self,o):o=R(o);return R(tuple(a+b for a,b in zip(self.v,o.v)))
    __radd__=__add__
    def __neg__(self):return R(tuple(-a for a in self.v))
    def __sub__(self,o):return self+-R(o)
    def __rsub__(self,o):return R(o)+-self
    def __mul__(self,o):
        o=R(o);v=[F(0)]*8
        for i,a in enumerate(self.v):
            if not a:continue
            for j,b in enumerate(o.v):
                if b:v[i^j]+=a*b*RADS[i&j]
        return R(v)
    __rmul__=__mul__
    def __truediv__(self,q):return R(tuple(a/F(q) for a in self.v))
    def __eq__(self,o):return self.v==R(o).v
    def __lt__(self,o):return sign((self-R(o)).v)<0
    def __hash__(self):return hash(self.v)
    def record(self):return {str(rad):str(a) for rad,a in zip(RADS,self.v) if a}
def sqrt(n):
    if n in (1,4):return R(isqrt(n))
    v=[F(0)]*8;v[RADS.index(n)]=F(1);return R(v)
S=sqrt(2)
def rectangles(points):
    rects=[];counts={};empty=0
    for mask in range(1,1<<len(points)):
        pts=[p for i,p in enumerate(points) if mask>>i&1]
        n=len(pts);ys=[y for x,y in pts];xs=[x for x,y in pts]
        span=max(ys)-min(ys);mixed=min(xs)!=max(xs)
        if n<=span*span or (mixed and n<=2):empty+=1;continue
        assert n<=6
        r=sqrt(n)/2
        x0=R(0) if min(xs)==0 else S;x1=R(0) if max(xs)==0 else S
        rect=(x1-r,x0+r,R(max(ys))-r,R(min(ys))+r)
        assert rect[0]<rect[1] and rect[2]<rect[3]
        rects.append(rect)
        key=f'{"mixed" if mixed else "one_column"}_span{span}_mass{n}'
        counts[key]=counts.get(key,0)+1
    return rects,{'all_nonempty_subsets':(1<<len(points))-1,'nonempty_receivers':len(rects),'empty_receivers':empty,'classes':counts}
def union_area(rects):
    xs=sorted(set(t for r in rects for t in r[:2]));area=R(0)
    for x0,x1 in zip(xs,xs[1:]):
        ints=sorted((b,t) for l,r,b,t in rects if l<x1 and x0<r)
        if not ints:continue
        lo,hi=ints[0];length=R(0)
        for a,b in ints[1:]:
            if a<=hi:hi=max(hi,b)
            else:length+=hi-lo;lo,hi=a,b
        length+=hi-lo;area+=(x1-x0)*length
    return area
rows=[]
for n in range(2,7):
    old=[(0,j) for j in range(n)];new=old+[(1,j) for j in range(n)]
    ro,co=rectangles(old);rn,cn=rectangles(new)
    ao=union_area(ro);an=union_area(rn)
    k=S-1;e=3*S-4
    assert ao==n+(n-1)*k*k
    assert an-2*ao==(n-1)*e
    assert an/(2*n)==1+(1-F(1,n))*(1-S/2)
    if n>=3:assert an/(2*n)>4-2*S
    rows.append({'N':n,'old_enumeration':co,'new_enumeration':cn,'old_area':ao.record(),'new_area':an.record(),'true_added_area':(an-2*ao).record(),'new_ratio':(an/(2*n)).record()})
assert (5*S-7)/3>F(1,50)
assert F(12,1000)<F(1,50)
assert F(353,250)**2<2
here=Path(__file__).resolve();out=here.parent/'results/reflection-column/exact-unions.json';out.parent.mkdir(parents=True,exist_ok=True)
result={'status':'PASS','python':sys.version,'arithmetic':'Exact multiquadratic field Q(sqrt(2),sqrt(3),sqrt(5)); sign comparisons use directed rational isqrt intervals','comparison_max_bits':max(SIGN_BITS),'random_seed':None,'floating_tolerance':None,'threshold':1,'atom_mass':1,'N':list(range(2,7)),'scope':'All nonempty subsets, including unsaturated mass-3 and mass-5 mixed receivers; complete vertical-slab union area. General N is proved separately.','script_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'cases':rows,'L1_eta':'1/1000','L1_ratio_loss_upper':'3/250','N3_gain_over_infinite_old_template_lower':'1/50'}
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','output':str(out),'N':result['N'],'max_subsets':rows[-1]['new_enumeration']['all_nonempty_subsets'],'comparison_max_bits':result['comparison_max_bits']},ensure_ascii=False))
