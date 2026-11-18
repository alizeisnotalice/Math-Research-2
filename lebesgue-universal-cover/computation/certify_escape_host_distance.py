"""Rational support transfer from the R108 escape host to all R95 orientations."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

from certify_b9_padding_four_minimizer import Hull, balance
from certify_nine import S, ar, append_point, arc_body, hull

here = Path(__file__).resolve().parent
ctx.prec = 160
pi = arb.pi()
input_path = here / "escape-host-distance-input.json"
pose_path = here / "competing-fits-input.json"
spec = json.loads(input_path.read_bytes())
pose = json.loads(pose_path.read_bytes())["zero_gain_five"]["params"]
eight = here / "four-body-contact-stationary.json"


def rotate(point, c, s):
    x, y = point
    return x*c-y*s, x*s+y*c


def dot(x, y):
    return sum(a*b for a,b in zip(x,y,strict=True))


def det(x, y):
    return x[0]*y[1]-x[1]*y[0]


host = Hull()
assert spec["normal_seed_denominator"] == 10**8
normals=[host.support((sign,q*10**4)) for sign,q in spec["normal_seeds"]]
samples=[]
for j in range(256):
    angle=2*pi*j/256
    append_point(samples,angle.cos()/2,angle.sin()/2)
arc_body(samples,128,pi,3,arb(0),*map(ar,pose[:2]))
for n,index in ((5,2),(7,5),(9,8)):
    arc_body(samples,128,pi,n,*map(ar,pose[index:index+3]))
denominator=S+4
integer_hull=hull(samples)
material=[tuple(F(a,denominator) for a in point) for point in integer_hull]
assert len(material)>=3
# Rounded points are within 2/S of the escape hull, which contains the
# radius-1/2 disk. Thus dividing their coordinates by (1+4/S)
# places each point inside the actual escape hull.
assert all(x*x+y*y<1 for x,y in material)
approximate=np.array(material,dtype=float)

total_cells=spec["angular_cells"]
used=set()
worst=arb(1)
accepted=0
for sign in (-1,1):
    rotated_material=approximate*[1,sign]
    data=spec["cell_support_indices"][str(sign)]
    assert len(data)==total_cells

    def evaluate(left,right,indices):
        assert len(indices)==len(set(indices)) in (2,3)
        theta=ar(left+right)*pi
        width=ar(right-left)*pi
        c,s=theta.cos(),theta.sin()
        cc,ss=float(c),float(s)
        floated=rotated_material@np.array([[cc,ss],[-ss,cc]])
        directions=[normals[j][0] for j in indices]
        if len(directions)==2:
            assert all(x==-y for x,y in zip(*directions,strict=True))
            weights=[F(1,2),F(1,2)]
        else:
            weights=balance(directions)
        A=B=C=arb(0)
        selected=[]
        for j,w in zip(indices,weights,strict=True):
            d,u,_,h=normals[j]
            p=int(np.argmax(floated@np.array(d,dtype=float)))
            selected.append(p)
            x,y=map(ar,material[p])
            turned=rotate((x,sign*y),c,s)
            A+=ar(w)*dot(turned,u)
            B+=ar(w)*det(turned,u)
            C+=ar(w)*h
        assert A>0
        slack=A*width.cos()-abs(B)*width.sin()-C
        return slack,selected

    def verify_cell(left,right,choices,depth=0):
        global accepted,worst
        selected=None
        for indices in choices:
            slack,selected=evaluate(left,right,indices)
            if slack>ar(F("0.012")):
                worst=min(worst,slack.lower())
                used.update(selected)
                accepted+=1
                return
        assert depth<4,(sign,left,right,slack)
        middle=(left+right)/2
        verify_cell(left,middle,choices,depth+1)
        verify_cell(middle,right,choices,depth+1)

    for i,indices in enumerate(data):
        choices=[indices]+[data[(i+j)%total_cells] for j in (-1,1,-2,2)]
        verify_cell(F(i,total_cells),F(i+1,total_cells),choices)

point_numbers=(( -477094447069,173648177667),(-253856652971,439692620786),
               (88163490354,500000000000),(388930956715,326351822333),(507713305943,0),
               (388930956715,-326351822333),(88163490354,-500000000000),
               (-253856652971,-439692620786),(-477094447069,-173648177667))
point_denominator=10**12+10
assert all((ax-bx)**2+(ay-by)**2<point_denominator**2
           for ax,ay in point_numbers for bx,by in point_numbers)
angle,tx,ty=map(ar,spec["inside_nine_pose_exact"])
cc,ss=angle.cos(),angle.sin()
nine=[]
for x,y in point_numbers:
    u,v=rotate((ar(F(x,point_denominator)),ar(F(y,point_denominator))),cc,ss)
    nine.append((u+tx,v+ty))
least=arb(1)
for (x,y),(u,v) in zip(integer_hull,integer_hull[1:]+integer_hull[:1],strict=True):
    a,b=v-y,x-u
    length=arb(a*a+b*b).sqrt()
    offset=ar(F(a*x+b*y,denominator))
    for px,py in nine:
        slack=(offset-a*px-b*py)/length
        assert slack>ar(F("0.000035"))
        least=min(least,slack.lower())

escape_receipt = here / "competing-fits-replay.json"
escape_results = json.loads(escape_receipt.read_bytes())["five_with_s_strictly_inside"]
assert F(escape_results["minimum_support_slack"])>0
assert F(escape_results["area_upper"]) < F("0.837657369")

print(json.dumps({
    "python":sys.version, "numpy":np.__version__, "precision_bits":ctx.prec,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),
    "escape_pose_sha256":hashlib.sha256(pose_path.read_bytes()).hexdigest(),
    "escape_replay_sha256":hashlib.sha256(escape_receipt.read_bytes()).hexdigest(),
    "escape_area_upper":"0.837657369",
    "four_body_minimizer_receipt_sha256":hashlib.sha256(eight.read_bytes()).hexdigest(),
    "normal_count":len(normals),"initial_angle_cells_by_chirality":total_cells,"total_accepted_angle_cells":accepted,
    "escape_rational_inner_hull_vertices":len(material),
    "rational_witness_material_used":len(used),
    "directed_padding_strictly_exceeds":"0.012",
    "directed_padding_certified_ball_lower":str(worst),
    "P9_strict_inside_escape_host_margin_lower":str(least),
    "passed":True,
},indent=2))
