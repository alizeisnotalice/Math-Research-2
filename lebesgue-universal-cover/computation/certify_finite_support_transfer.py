"""Rational Fourier samples with positive and negative local transfer controls."""
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys

from flint import arb, ctx

from certify_five_harmonic import old_common_polygon
from certify_large_support import curvature
from certify_nine import PI_UPPER, S, ar, append_point, arc_body, hull

here = Path(__file__).resolve().parent
ctx.prec = 160
pi = arb.pi()
assert pi < ar(PI_UPPER)
source = here / "sparse-eight-escape-input.json"
fit = here / "competing-fits-input.json"
box = here / "five-body-contact-basin.json"
reference = here / "five-body-contact-stationary.json"
spec = json.loads(source.read_bytes())
terms = [(n,F(a)) for n,a in json.loads(fit.read_bytes())["terms"]]
assert curvature(terms,pi) == F(1,100)
# Width one makes rho(theta+pi)=1-rho(theta), hence rho<=0.99.
assert F(1,2)-sum(abs(a) for n,a in terms) > F("0.46")
denominator = S+10


def rotate(p,c,s):
    x,y = p
    return x*c-y*s,x*s+y*c


def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def sample(count,phase):
    points=[]
    for j in range(count):
        theta=2*pi*ar((j+phase)/count)
        h=ar(F(1,2))+sum(ar(a)*(n*theta).cos() for n,a in terms)
        hp=-sum(ar(n*a)*(n*theta).sin() for n,a in terms)
        c,s=theta.cos(),theta.sin()
        append_point(points,h*c-hp*s,h*s+hp*c)
    squared=max((x-u)**2+(y-v)**2 for x,y in points for u,v in points)
    assert squared < denominator**2
    # Grid points are within 2/S of S; since S contains the radius
    # 0.46 disk the 10/S rescaling puts the rational samples in S.
    return points,squared


def fit_entire_host_box(points,placement):
    box_receipt=json.loads(box.read_bytes())
    center=list(map(F,box_receipt["raw_center"]))
    half=list(map(F,box_receipt["raw_halfwidth"]))
    polygon,shrink=old_common_polygon(center,half,pi,k=384)
    ch=placement["reflection_sign"]
    assert ch in (1,-1)
    a=ar(placement["angle"])
    dx,dy=map(ar,placement["translation"])
    values=[]
    for x,y in points:
        u,v=rotate((ar(F(x,denominator)),ar(F(ch*y,denominator))),a.cos(),a.sin())
        values.append((u+dx,v+dy))
    margin=arb(1)
    for (x,y),(u,v) in zip(polygon,polygon[1:]+polygon[:1],strict=True):
        nx,ny=v-y,x-u
        length=arb(nx*nx+ny*ny).sqrt()
        offset=ar(F(nx*x+ny*y,S))
        for a,b in values:
            slack=(offset-nx*a-ny*b)/length
            assert slack>ar(F("0.00188"))
            margin=min(margin,slack.lower())
    return margin,len(polygon)


def transfer(old_lower,count):
    distance=F(99,100)*PI_UPPER**2/(2*count*count)+F(12,S)
    lower=(old_lower-PI_UPPER*distance**2)/(1+4*distance)
    return distance,lower


def first_count(old_lower,target):
    lo,hi=1,2
    while transfer(old_lower,hi)[1] <= target:
        hi*=2
    while hi-lo>1:
        m=(lo+hi)//2
        if transfer(old_lower,m)[1] > target:
            hi=m
        else:
            lo=m
    return hi


def area_upper(params, points, pi):
    params=list(map(ar,params))
    grid=2048
    material=[]
    for j in range(8*grid):
        a=2*pi*j/(8*grid)
        append_point(material,a.cos()/2,a.sin()/2)
    arc_body(material,grid,pi,3,arb(0),*params[:2])
    for n,k in ((5,2),(7,5),(9,8)):
        arc_body(material,grid,pi,n,*params[k:k+3])
    a,x,y=params[11:]
    for p in points:
        u,v=rotate(tuple(map(ar,p)),a.cos(),a.sin())
        append_point(material,u+x,v+y)
    poly=hull(material)
    edges=list(zip(poly,poly[1:]+poly[:1],strict=True))
    area=F(sum(cross(a,b) for a,b in edges),2*S*S)
    per=0
    for (x,y),(u,v) in edges:
        square=(x-u)**2+(y-v)**2
        r=isqrt(square)
        per+=r+(r*r<square)
    epsilon=PI_UPPER**2/F(72*grid**2)+F(2,S)
    upper=area+F(per,S)*epsilon+PI_UPPER*epsilon**2
    assert upper<F("0.83729965")
    return upper,grid


assert spec["n_points"] == 8
eight,diameter2 = sample(8,F(spec["seed_arc_phase"]))
margin,vertices=fit_entire_host_box(eight,spec["common_inner_placement"])
upper,arc_count=area_upper(spec["params"],[tuple(F(z,denominator) for z in p)
                                         for p in eight],pi)
stationary=json.loads(reference.read_bytes())
assert arb(stationary["area_ball"]) < ar(F("0.837295006222"))
raw_box=json.loads(box.read_bytes())
for value,c,w in zip(stationary["placement"],raw_box["raw_center"],
                     raw_box["raw_halfwidth"],strict=True):
    c,w=F(c),F(w)
    assert ar(c-w)<arb(value)<ar(c+w)
target=F("0.83730")
narrow=here/"basin-083945-parent.json"
wide=here/"wide-083730-parent.json"
narrow_area=F("0.83945005239")
assert F(json.loads(narrow.read_bytes())["lower"]) > narrow_area
count=first_count(narrow_area,target)
assert count==88
dense,dense_diameter2=sample(count,F(0))
delta,local_floor=transfer(narrow_area,count)
assert local_floor>F("0.8373341036")
wide_area=F("0.83730023870")
assert F(json.loads(wide.read_bytes())["lower"]) > wide_area
wide_count=first_count(wide_area,target)
assert wide_count==8284
print(json.dumps({
    "python":sys.version, "arb_precision_bits":ctx.prec,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in (source,fit,box,reference,narrow,wide)},
    "point_denominator":denominator,
    "eight_point_coordinate_numerators":eight,
    "eight_point_squared_diameter_numerator_exact":str(diameter2),
    "eight_points_fit_all_R102_hosts_same_pose_margin_lower":str(margin),
    "common_R102_inner_polygon_vertex_count":vertices,
    "eight_point_free_full_five_area_upper_exact":str(upper),
    "eight_point_free_full_five_area_upper_decimal":"0.83729965",
    "eight_points_at_R99_stationary_host_area_upper":"0.837295006222",
    "area_upper_cells_per_original_arc":arc_count,
    "dense_point_coordinate_numerators":dense,
    "dense_point_squared_diameter_numerator_exact":str(dense_diameter2),
    "dense_support_loss_upper_exact":str(delta),
    "dense_full_five_area_lower_on_R102_box_exact":str(local_floor),
    "dense_full_five_area_lower_on_R102_box_decimal":"0.8373341036",
    "minimum_N_for_this_uniform_bound_at_R102":"88",
    "minimum_N_for_this_uniform_bound_at_R108_wide_box":"8284",
    "passed":True,
},indent=2))
