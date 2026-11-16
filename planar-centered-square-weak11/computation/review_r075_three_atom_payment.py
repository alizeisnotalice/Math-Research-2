"""One independent exact radical horizontal-slab audit of the directed fixture."""
from pathlib import Path
from functools import cmp_to_key
import hashlib, itertools, json, os, platform, sys, time
import sympy as S

start=time.monotonic(); source=Path(sys.argv[1]); raw=source.read_bytes(); data=json.loads(raw)
points=[tuple(S.Rational(x) for x in p) for p in data['positions']]
masses=list(map(S.Rational,data['masses']))
assert points==[(0,0),(-2,0),(0,-2)] and masses==[1,4,4]
def cmp(a,b):
    d=S.simplify(a-b)
    if d==0:return 0
    if bool(d<0):return -1
    assert bool(d>0);return 1
def ordered(xs):return sorted(set(map(S.simplify,xs)),key=cmp_to_key(cmp))
def less(a,b):return cmp(a,b)<0
def leq(a,b):return cmp(a,b)<=0
def minimum(a,b):return a if leq(a,b) else b
def maximum(a,b):return b if less(a,b) else a
def merge(intervals):
    out=[]
    for a,b in sorted(intervals,key=cmp_to_key(lambda i,j:cmp(i[0],j[0]))):
        if not out or less(out[-1][1],a):out.append([a,b])
        else:out[-1][1]=maximum(out[-1][1],b)
    return out
catalog=[];full=[];res=[]
for mask in range(1,8):
    ids=[i for i in range(3) if mask>>i&1]; mass=sum(masses[i] for i in ids)
    rad=S.sqrt(mass)/2
    x0=max(points[i][0] for i in ids)-rad; x1=min(points[i][0] for i in ids)+rad
    y0=max(points[i][1] for i in ids)-rad; y1=min(points[i][1] for i in ids)+rad
    rect=tuple(map(S.simplify,(x0,minimum(x1,0),y0,minimum(y1,0))))
    assert less(rect[0],rect[1]) and less(rect[2],rect[3])
    full.append(rect)
    if not mask&1:res.append(rect)
    catalog.append({'mask':mask,'mass':str(mass),'SW_receiver':list(map(str,rect))})
assert len(full)==7 and len(res)==3
def slabs(fr,rr):
    yy=ordered([v for rect in fr+rr for v in rect[2:]]); areas=[]
    for y0,y1 in zip(yy,yy[1:]):
        ymid=(y0+y1)/2
        f=merge([(a,b) for a,b,c,d in fr if less(c,ymid) and less(ymid,d)])
        r=merge([(a,b) for a,b,c,d in rr if less(c,ymid) and less(ymid,d)])
        length=sum((b-a for a,b in f),S.Integer(0))
        for a,b in f:
            for c,d in r:
                lo=maximum(a,c); hi=minimum(b,d)
                if less(lo,hi):length-=hi-lo
        assert leq(0,length)
        areas.append((y0,y1,S.simplify((y1-y0)*length)))
    return S.simplify(sum((a for y0,y1,a in areas),S.Integer(0))),areas
area,rows=slabs(full,res)
u=S.sqrt(5)/2; v=S.sqrt(2); formula=4*v+u-S.Rational(25,4)
assert S.simplify(area-formula)==0 and bool(area>S.Rational(131,250))
old=S.Rational(1,4); stems=u-1; central=4*v-S.Rational(11,2)
assert S.simplify(old+stems+central-area)==0
# The residual pair rectangle must not be dropped.
wrong,_=slabs(full,[rect for row,rect in zip(catalog,full) if row['mask'] in (2,4)])
assert bool(wrong>area)
def clip(rect,lo,hi):
    a,b,c,d=rect;z=(maximum(a,lo),minimum(b,hi),maximum(c,lo),minimum(d,hi))
    return z if less(z[0],z[1]) and less(z[2],z[3]) else None
def clipped_area(lo,hi):
    ff=[z for r in full if (z:=clip(r,lo,hi)) is not None]
    rr=[z for r in res if (z:=clip(r,lo,hi)) is not None]
    return slabs(ff,rr)[0]
near=clipped_area(-S.Rational(1,2),0)
gamma3=clipped_area(-S.Rational(1,2),-S.Rational(1,9))
assert near==S.Rational(1,4) and gamma3==S.Rational(49,324)
assert S.Rational(707,500)**2<2 and S.Rational(559,500)**2<S.Rational(5,4)
assert 4*S.Rational(707,500)+S.Rational(559,500)-S.Rational(25,4)==S.Rational(131,250)
result={'status':'PASS','algorithm':'exact radical horizontal y slabs and merged x interval difference; no producer imports or numeric comparison',
 'full_subsets':7,'residual_subsets':3,'catalog':catalog,
 'horizontal_slabs':[{'y_lower':str(a),'y_upper':str(b),'deletion_area':str(c)} for a,b,c in rows],
 'whole_SW_area':str(area),'whole_SW_decimal':float(area),'analytic_partition':{'old_singleton':'1/4','two_stems':str(stems),'central_square_after_complete_residual':str(central)},
 'simple_strict_lower':'131/250','strict_excess_over_half_mass':str(S.simplify(area-S.Rational(1,2))),
 'wrong_area_if_residual_pair_omitted':str(wrong),'overcount_without_pair':str(S.simplify(wrong-area)),
 'clipped_minus_half_to_zero_area':str(near),'clipped_minus_half_to_minus_one_ninth_area':str(gamma3),
 'qualification':{'actual_maximum_eligible_mass':'9','original_internal_mass':'1','S_is_maximum':False,'unit_caps_hold':False,'general_caps_hold':False},
 'scope':'Uncapped whole-SW two-effective-outside collapse counterexample only. Clipped Gamma is different. No uniform-grid/full-cap/3.5/C2 counterexample.',
 'source_input':str(source),'source_input_sha256':hashlib.sha256(raw).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'sympy':S.__version__,'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('catalog','horizontal_slabs')}))
