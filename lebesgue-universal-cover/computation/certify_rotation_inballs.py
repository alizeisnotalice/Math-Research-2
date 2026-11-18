"""Rebuild rational disks inside all eight bodies and exact pair-distance exclusions."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import prod
from pathlib import Path
import platform

from rational_trig import PI, scale, trig

HERE = Path(__file__).resolve().parent
TAU = F(419, 500)
RADIUS_SCALE = 10**9
DISTANCE_SCALE = 10**4


def fourier_bounds():
    # h=1/2+p(x)/1200 with x=cos(theta). p' is a quadratic in x^2.
    def derivative(x):
        return -1520*x**4 + 1320*x*x - 140

    roots = []
    values = [F(0), F(4)]  # |p(0)| and |p(1)|
    for lo, hi in [(F(3, 10), F(4, 10)), (F(8, 10), F(9, 10))]:
        for _ in range(50):
            mid = (lo + hi) / 2
            if derivative(lo)*derivative(mid) < 0:
                hi = mid
            else:
                lo = mid
            assert derivative(lo)*derivative(hi) < 0
        # Bound |p| by positive-coefficient interval arithmetic.
        pmin = -304*hi**5 + 440*lo**3 - 140*hi
        pmax = -304*lo**5 + 440*hi**3 - 140*lo
        values.append(max(abs(pmin), abs(pmax)))
        roots.append((lo, hi, pmin, pmax))
    assert values[2] > values[3] and roots[0][3] < -values[3]
    assert roots[0][3] < -4
    return F(1, 2)-max(values)/1200, F(1, 2)+roots[0][3]/1200


def inradii():
    low, high = fourier_bounds()
    assert F(47355617284265,10**14) < low <= high < F(47355617284267,10**14)
    radii = {'D': F(1,2)}
    for n in [3,5,7,9,11,19]:
        c_lo, c_hi = trig(scale(PI,F(1,2*n)),'cos')
        if n < 11:
            r = 1 - 1/(2*c_lo)
            label = f'B{n}'
        else:
            r = trig(scale(PI,F(1,n)), 'cos')[0]/(2*c_hi)
            label = f'P{n}'
        radii[label] = F((r*RADIUS_SCALE).__floor__(),RADIUS_SCALE)
    radii['F'] = F((low*RADIUS_SCALE).__floor__(),RADIUS_SCALE)
    return radii


def hull(points):
    def lower(seq):
        chain = []
        for x,y in seq:
            while len(chain) >= 2:
                a,b = chain[-2]; c,d = chain[-1]
                if (c-a)*(y-b)-(d-b)*(x-a) > 0:
                    break
                chain.pop()
            chain.append((x,y))
        return chain[:-1]
    p = sorted(set(points))
    return lower(p)+lower(reversed(p))


def run():
    radii = inradii()
    base_file = HERE/'outer_domain_certificate.json'
    data = json.loads(base_file.read_bytes())
    den = int(data['disk_polygon_denominator'])
    circle = [tuple(map(int,p)) for p in data['disk_polygon_points']]
    assert len(circle)==512 and all(4*(x*x+y*y)<=den*den for x,y in circle)
    unit = den*RADIUS_SCALE
    polygons = {}
    for name,r in radii.items():
        f = 2*r*RADIUS_SCALE
        assert f.denominator == 1
        polygons[name] = [(x*f.numerator,y*f.numerator) for x,y in circle]

    pairs = {}
    distances = {}
    for a,b in combinations(radii,2):
        def area(digit):
            xshift = digit*unit//DISTANCE_SCALE
            assert xshift*DISTANCE_SCALE == digit*unit
            points = hull(polygons[a]+[(x+xshift,y) for x,y in polygons[b]])
            cross_sum = sum(x*v-y*u for (x,y),(u,v) in zip(points,points[1:]+points[:1]))
            return F(cross_sum,2*unit*unit)
        lo,hi=0,DISTANCE_SCALE
        while hi-lo>1:
            mid=(hi+lo)//2
            if area(mid)>TAU:
                hi=mid
            else:
                lo=mid
        d = F(hi,DISTANCE_SCALE)
        area_bound = area(hi)
        assert area_bound > TAU
        distances[a,b] = d
        pairs[a+'-'+b] = {'distance':str(d), 'inner_area':str(area_bound)}
    relative_box = prod((distances['D',b] / (F(6,25) if b=='F' else F(21,100)))**2 for b in radii if b!='D')
    t3lo,t3hi=F(1732050807568877,10**15),F(1732050807568878,10**15)
    assert t3lo**2<3<t3hi**2
    pair_area=[PI[0]/8+t3lo/4,PI[1]/8+t3hi/4]
    assert pair_area[0]>F(8257117835909,10**13)
    return {'python':platform.python_version(),'threshold':str(TAU),'distance_mesh':'1/10000',
            'radius_mesh':'1/10^9','disk_data_sha256':sha256(base_file.read_bytes()).hexdigest(),
            'radii':{k:str(v) for k,v in radii.items()},'pairs':pairs,
            'fourier_centered_inradius_enclosure':['0.47355617284265','0.47355617284267'],
            'translation_box_volume_fraction':str(relative_box),
            'translation_box_volume_display':float(relative_box),
            'disk_B3_optimum_enclosure':[str(z) for z in pair_area],
            'scope':'28 necessary radial constraints for area <=419/500 and exact two-body optimum; no eight-body minimum certified'}

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
