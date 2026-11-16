"""Exact finite fixture: D4 averaging need not preserve a deleted receiver cell."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, itertools, json, os, platform, sys, time
start=time.monotonic(); n=3; w=F(8,9); e=F(1,9)
old=[(F(4*i,9),F(4*j,9)) for i in range(3) for j in range(3)]
q=(-F(1,3),-F(1,4)); b=F(3,50)
orbit=sorted({(x,y) for a in (q,(q[1],q[0])) for x in (a[0],w-a[0]) for y in (a[1],w-a[1])})
assert len(orbit)==8
def audit(points,masses):
    xs=sorted({x for x,y in points}); ys=sorted({y for x,y in points}); closures={}
    for a,d in itertools.combinations_with_replacement(xs,2):
        for f,g in itertools.combinations_with_replacement(ys,2):
            mask=sum(1<<i for i,(x,y) in enumerate(points) if a<=x<=d and f<=y<=g)
            if not mask or mask in closures: continue
            inds=[i for i in range(len(points)) if mask>>i&1]
            width=max(points[i][0] for i in inds)-min(points[i][0] for i in inds)
            height=max(points[i][1] for i in inds)-min(points[i][1] for i in inds)
            mass=sum((masses[i] for i in inds),F(0)); diameter=max(width,height)
            assert mass<=max(F(1),diameter*diameter)
            closures[mask]=(mass,diameter)
    eligible={mask:v for mask,v in closures.items() if v[0]>v[1]*v[1]}
    maximum=max(m for m,d in eligible.values()); ties=[k for k,(m,d) in eligible.items() if m==maximum]
    assert maximum==1 and ties==[(1<<len(old))-1]
    component=set(range(len(old)))
    changed=True
    while changed:
        changed=False
        for mask in eligible:
            inds={i for i in range(len(points)) if mask>>i&1}
            if inds&component and not inds<=component: component|=inds; changed=True
    assert len(component)==len(points)
    diameter=max(max(xs)-min(xs),max(ys)-min(ys))
    gap=diameter**2-F(35,33)*sum(masses)
    assert gap>0
    return {'atoms':len(points),'closures':len(closures),'eligible':len(eligible),
            'unique_max_mask':ties[0],'nonterminal_gap':str(gap),'connected':True}
original=audit(old+[q],[e]*9+[b])
averaged=audit(old+orbit,[e]*9+[b/8]*8)
lo=-F(201,1000); hi=-F(199,1000); cell=(lo,hi,lo,hi)
r=F(101,500)
def far(p): return max(abs(p[0]-lo),abs(p[0]-hi),abs(p[1]-lo),abs(p[1]-hi))
assert max(far((0,0)),far(q))<r and e+b>4*r*r
def delta(points):
    xmin=min(x for x,y in points); xmax=max(x for x,y in points)
    ymin=min(y for x,y in points); ymax=max(y for x,y in points)
    return max((xmax-xmin)/2,(ymax-ymin)/2,lo-xmin,xmax-hi,lo-ymin,ymax-hi,F(0))
residual_margin=4*delta([q])**2-b; assert residual_margin>0
p=F(1,2)-e
local=[(point,mass) for point,mass in zip(old+orbit,[e]*9+[b/8]*8)
       if -1<point[0]<p and -1<point[1]<p]
assert len(local)==3 and local[0][0]==(0,0)
subset_margins=[]
for mask in range(1,1<<len(local)):
    selected=[point for i,(point,mass) in enumerate(local) if mask>>i&1]
    mass=sum((mass for i,(point,mass) in enumerate(local) if mask>>i&1),F(0))
    margin=4*delta(selected)**2-mass; assert margin>0
    subset_margins.append(str(margin))
result={'status':'PASS','n':3,'outside_point':list(map(str,q)),'outside_mass':str(b),
        'd4_orbit':[[str(x),str(y)] for x,y in orbit], 'original':original,'d4_average':averaged,
        'receiver_cell':list(map(str,cell)),'receiver_cell_area':str((hi-lo)**2),
        'original_common_open_radius':str(r),'original_full_margin':str(e+b-4*r*r),
        'original_residual_exclusion_margin':str(residual_margin),
        'averaged_relevant_sources':3,'averaged_all_nonempty_local_subset_exclusion_margins':subset_margins,
        'conclusion':'Original deletion contains this positive-area cell; D4-average full halo misses the entire cell. Both inputs have exact caps, unique S3 maximum, connectivity and nonterminal strictness.',
        'scope':'Failure of pointwise/receiver-cell preservation only; no claim about total corner area comparison, mass-half bound, or main corner bound.',
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result))
