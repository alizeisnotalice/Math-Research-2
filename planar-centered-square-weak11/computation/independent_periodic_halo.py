#!/usr/bin/env python3
"""Independent exact enclosure of a periodic centered-square atomic halo.

Masses are accumulated anew for each left x endpoint by adding complete x
columns to y rows. There is no 2-D prefix table. Rectangle area uses an x-slab
sweep and explicit merging of a sorted multiset of active y intervals; there
is no segment tree. All arithmetic used for the enclosure is integer/rational.

For a closed coordinate box B of mass M, centers whose distance to every
point of B is less than R=sqrt(M/(4*lambda)) form the open rectangle
[xmax-R,xmin+R] x [ymax-R,ymin+R]. Capping R by Rmax suffices after an exact
periodic tail bound proves no radius >=Rmax can witness the strict halo.
The bounds below use a strict lower dyadic radius and a strict upper radius
before capping; open/closed endpoint choices do not alter area.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import bisect
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import sys
import time


def frac(value):
    return F(str(value))


def floor(q):
    return q.numerator // q.denominator


def ceil(q):
    return -((-q.numerator) // q.denominator)


def lcm_all(values):
    return math.lcm(*values)


def hash_rectangles(rectangles):
    h = hashlib.sha256()
    for rectangle in sorted(rectangles):
        h.update((','.join(map(str, rectangle)) + '\n').encode())
    return h.hexdigest()


def slab_union_area(rectangles):
    """Exact union area; merge all active y intervals in each x slab."""
    events = defaultdict(list)
    for x0, x1, y0, y1 in rectangles:
        assert x0 < x1 and y0 < y1
        events[x0].append((1, (y0, y1)))
        events[x1].append((-1, (y0, y1)))
    if not events:
        return 0, {'x_events': 0, 'peak_active_intervals': 0}
    multiplicity = Counter()
    ordered = []
    area = 0
    length = 0
    previous = min(events)
    peak = 0
    for x, changes in sorted(events.items()):
        area += (x - previous) * length
        for sign, interval in changes:
            old = multiplicity[interval]
            new = old + sign
            assert new >= 0
            if old == 0:
                bisect.insort(ordered, interval)
            if new == 0:
                index = bisect.bisect_left(ordered, interval)
                assert ordered[index] == interval
                ordered.pop(index)
                del multiplicity[interval]
            else:
                multiplicity[interval] = new
        peak = max(peak, len(ordered))
        length = 0
        if ordered:
            low, high = ordered[0]
            for nextlow, nexthigh in ordered[1:]:
                if nextlow > high:
                    length += high - low
                    low, high = nextlow, nexthigh
                elif nexthigh > high:
                    high = nexthigh
            length += high - low
        previous = x
    assert not ordered and not multiplicity and length == 0
    return area, {'x_events': len(events), 'peak_active_intervals': peak}


def cell_area(rectangles):
    """Slower independent coordinate-cell union used only in tiny tests."""
    if not rectangles:
        return 0
    xs = sorted({r[i] for r in rectangles for i in (0, 1)})
    ys = sorted({r[i] for r in rectangles for i in (2, 3)})
    area = 0
    for a, b in zip(xs, xs[1:]):
        for c, d in zip(ys, ys[1:]):
            if any(x0 <= a and b <= x1 and y0 <= c and d <= y1
                   for x0, x1, y0, y1 in rectangles):
                area += (b - a) * (d - c)
    return area


def setup_grid(points, periods, bits, radius_cap, lam):
    D = 1 << bits
    S = lcm_all([D, radius_cap.denominator, *(q.denominator for q in periods),
                 *(q.denominator for p in points for q in p[:2])])
    W = lcm_all([p[2].denominator for p in points])
    integer_points = [(int(x*S), int(y*S), int(m*W)) for x,y,m in points]
    return {'D': D, 'S': S, 'W': W, 'lam': lam, 'cap': int(radius_cap*S),
            'cell': (0, int(periods[0]*S), 0, int(periods[1]*S)),
            'points': integer_points}


def rounded_radii(mass, grid):
    """isqrt of a floored rational equals floor of its positive square root."""
    D, S, W, lam = (grid[k] for k in ('D', 'S', 'W', 'lam'))
    numerator = mass * D * D * lam.denominator
    denominator = 4 * W * lam.numerator
    root = math.isqrt(numerator // denominator)
    # Explicitly verify the floor-root inequality; no floating point involved.
    assert root * root * denominator <= numerator
    assert (root + 1) * (root + 1) * denominator > numerator
    quantum = S // D
    return (min(grid['cap'], max(0, root - 1) * quantum),
            min(grid['cap'], (root + 1) * quantum))


def rectangle_for_bounds(bounds, radius, cell):
    xmin, xmax, ymin, ymax = bounds
    cx0, cx1, cy0, cy1 = cell
    result = (max(cx0, xmax-radius), min(cx1, xmin+radius),
              max(cy0, ymax-radius), min(cy1, ymin+radius))
    return result if result[0] < result[1] and result[2] < result[3] else None


def enumerate_boxes(grid, progress=False):
    """All coordinate boxes with spans <=2*cap, via x-column accumulation."""
    xs = sorted({x for x,y,m in grid['points']})
    ys = sorted({y for x,y,m in grid['points']})
    yindex = {y:i for i,y in enumerate(ys)}
    columns = defaultdict(list)
    for x,y,m in grid['points']:
        columns[x].append((yindex[y],m))
    ymax = []
    for i,y in enumerate(ys):
        ymax.append(bisect.bisect_right(ys,y+2*grid['cap']))
    cache = {}
    inner, outer = set(), set()
    counts = Counter({'positive_outer_with_empty_inner': 0})
    for li, left in enumerate(xs):
        row_mass = [0]*len(ys)
        for right in xs[li:]:
            if right-left > 2*grid['cap']:
                break
            for yi,m in columns[right]:
                row_mass[yi] += m
            counts['x_intervals'] += 1
            for low, bottom in enumerate(ys):
                mass = 0
                for hi in range(low, ymax[low]):
                    counts['all_span_admissible_boxes'] += 1
                    mass += row_mass[hi]
                    if not mass:
                        counts['empty_boxes'] += 1
                        continue
                    counts['positive_mass_boxes'] += 1
                    radii = cache.get(mass)
                    if radii is None:
                        radii = rounded_radii(mass,grid)
                        cache[mass] = radii
                    bounds = left, right, bottom, ys[hi]
                    inside = rectangle_for_bounds(bounds,radii[0],grid['cell'])
                    outside = rectangle_for_bounds(bounds,radii[1],grid['cell'])
                    if inside is not None:
                        assert outside is not None
                        assert outside[0] <= inside[0] < inside[1] <= outside[1]
                        assert outside[2] <= inside[2] < inside[3] <= outside[3]
                        inner.add(inside)
                        counts['positive_inner_boxes'] += 1
                    if outside is not None:
                        outer.add(outside)
                        counts['positive_outer_boxes'] += 1
                        if inside is None:
                            counts['positive_outer_with_empty_inner'] += 1
        if progress and li%20 == 0:
            print(json.dumps({'x_left_completed':li+1,'x_coordinates':len(xs),
                              'inner_unique':len(inner),'outer_unique':len(outer)}),flush=True)
    counts.update({'x_coordinates':len(xs),'y_coordinates':len(ys),
                   'unique_masses':len(cache),'unique_inner_rectangles':len(inner),
                   'unique_outer_rectangles':len(outer)})
    # Rectangle inclusion is checked for every box, not by equality of tuple sets.
    return inner, outer, dict(counts)


def subset_rectangles(grid):
    result = [set(),set()]
    points = grid['points']
    for mask in range(1,1<<len(points)):
        selected = [p for i,p in enumerate(points) if mask>>i&1]
        bounds = (min(p[0] for p in selected),max(p[0] for p in selected),
                  min(p[1] for p in selected),max(p[1] for p in selected))
        if bounds[1]-bounds[0]>2*grid['cap'] or bounds[3]-bounds[2]>2*grid['cap']:
            continue
        for index,radius in enumerate(rounded_radii(sum(p[2] for p in selected),grid)):
            rectangle = rectangle_for_bounds(bounds,radius,grid['cell'])
            if rectangle is not None:
                result[index].add(rectangle)
    return result


def run_tests():
    records = []
    examples = [
        [(F(0),F(0),F(4))],
        [(F(0),F(0),F(1)),(F(1),F(0),F(1))],
        [(F(0),F(0),F(1)),(F(0),F(0),F(2)),(F(1),F(1),F(3))],
        [(F(-1),F(0),F(1)),(F(1),F(0),F(3)),(F(0),F(1),F(1))],
        [(F(0),F(0),F(100)),(F(10),F(0),F(1)),(F(0),F(10),F(1))],
        [(F(0),F(0),F(2)),(F(2),F(0),F(2))],
    ]
    rng = random.Random(20261114)
    for _ in range(12):
        examples.append([(F(rng.randrange(-2,5),2),F(rng.randrange(-2,5),2),
                          F(rng.randrange(1,8),3)) for __ in range(5)])
    for index,points in enumerate(examples):
        grid = setup_grid(points,(F(3),F(3)),6,F(5),F(1))
        inner,outer,counts = enumerate_boxes(grid)
        subsets = subset_rectangles(grid)
        if index == 5:
            assert counts['positive_outer_with_empty_inner'] > 0
        areas=[]
        for boxes,sub in zip((inner,outer),subsets):
            a,_=slab_union_area(boxes)
            b,_=slab_union_area(sub)
            assert a==b==cell_area(boxes)==cell_area(sub)
            areas.append(str(F(a,grid['S']**2)))
        records.append({'case':index,'atoms':len(points),'areas':areas,
                        'all_subsets_vs_boxes':'equal','slab_vs_cells':'equal',
                        'positive_outer_with_empty_inner':counts['positive_outer_with_empty_inner']})
    # Exact square-root thresholds, including the cap at 5.
    grid=setup_grid([(F(0),F(0),F(1))],(F(3),F(3)),6,F(5),F(1))
    expected={4:(63,65),100:(319,320),101:(320,320)}
    for mass, radius_units in expected.items():
        got=rounded_radii(mass,grid)
        assert got==radius_units,(mass,got,radius_units)
    # A rectangle crossing and touching the cell boundary is clipped exactly.
    assert rectangle_for_bounds((0,0,0,0),64,(0,192,0,192))==(0,64,0,64)
    assert rectangle_for_bounds((0,128,0,0),64,(0,192,0,192)) is None
    # Every admitted coordinate pair includes equality at a span of 2*cap.
    assert records[4]['all_subsets_vs_boxes']=='equal'
    return {'seed':20261114,'tiny_examples':records,
            'threshold_tests':{'R=1':'passed','R=5':'passed','R>5':'passed',
                               'clipping_and_zero_width':'passed'}}


def rational_record(q):
    return {'fraction':str(q),'decimal_approximation':float(q)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--bits',type=int,default=40)
    parser.add_argument('--radius-cap',type=F,default=F(5))
    parser.add_argument('--compare-primary',type=Path)
    parser.add_argument('--compare-rectangles',type=Path)
    args=parser.parse_args()
    start=time.perf_counter()
    started_at=datetime.now(timezone.utc).isoformat()
    raw=args.input.read_bytes()
    original=json.loads(raw)
    px,py=frac(original['Px']),frac(original['Py'])
    lam=frac(original['lam'])
    points=[tuple(map(frac,p)) for p in original['pts']]
    assert px>0 and py>0 and lam>0 and args.radius_cap>0 and args.bits>=2
    assert all(0<=x<px and 0<=y<py and m>0 for x,y,m in points)
    base_mass=sum((p[2] for p in points),F(0))
    density=base_mass/(px*py)
    cap=args.radius_cap
    tail=density*(1+px/(2*cap))*(1+py/(2*cap))
    assert tail<lam, 'The proposed radius cap does not exclude all large radii.'
    lifted=[]
    for x,y,m in points:
        for k in range(ceil((-cap-x)/px),floor((px+cap-x)/px)+1):
            for ell in range(ceil((-cap-y)/py),floor((py+cap-y)/py)+1):
                xx,yy=x+k*px,y+ell*py
                assert -cap<=xx<=px+cap and -cap<=yy<=py+cap
                lifted.append((xx,yy,m))
    tests=run_tests()
    grid=setup_grid(lifted,(px,py),args.bits,cap,lam)
    inner,outer,counts=enumerate_boxes(grid,True)
    enumeration_seconds=time.perf_counter()-start
    a,inner_sweep=slab_union_area(inner)
    b,outer_sweep=slab_union_area(outer)
    assert 0<=a<=b<=grid['cell'][1]*grid['cell'][3]
    lower,upper=F(a,grid['S']**2),F(b,grid['S']**2)
    ratios=(lam*lower/base_mass,lam*upper/base_mass)
    output={
        'algorithm':'independent x-column accumulation and active interval merging',
        'input_path':str(args.input.resolve()),'input_sha256':hashlib.sha256(raw).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':sys.version,'platform':platform.platform(),'bits':args.bits,
        'started_at_utc':started_at,
        'radius_cap':str(cap),'base_atoms':len(points),'lifted_atoms':len(lifted),
        'periods':[str(px),str(py)],'lambda':str(lam),'base_mass':str(base_mass),
        'density':rational_record(density),'tail_bound_at_cap':rational_record(tail),
        'strict_tail_margin':rational_record(lam-tail),
        'coordinate_grid_denominator':str(grid['S']),'mass_grid_denominator':str(grid['W']),
        'counts':counts,'inner_rectangle_sha256':hash_rectangles(inner),
        'outer_rectangle_sha256':hash_rectangles(outer),
        'inner_sweep':inner_sweep,'outer_sweep':outer_sweep,
        'area_lower':rational_record(lower),'area_upper':rational_record(upper),
        'area_gap':rational_record(upper-lower),
        'weak_ratio_lower':rational_record(ratios[0]),
        'weak_ratio_upper':rational_record(ratios[1]),
        'weak_ratio_gap':rational_record(ratios[1]-ratios[0]),
        'tests':tests,'enumeration_seconds':enumeration_seconds,
        'total_seconds':time.perf_counter()-start,
        'command':[sys.executable,*sys.argv],
        'evidence_status':'Exact rational finite enumeration; periodic completeness uses the documented tail and coordinate-box reduction.',
    }
    if args.compare_primary:
        primary_raw=args.compare_primary.read_bytes()
        primary=json.loads(primary_raw)
        checks={
            'input_sha256':primary['input_sha256']==output['input_sha256'],
            'base_mass':F(primary['mass'])==base_mass,
            'tail_bound':F(primary['tail_bound'])==tail,
            'area_lower':F(primary['inner_area'])==lower,
            'area_upper':F(primary['outer_area'])==upper,
            'ratio_lower':F(primary['inner_ratio'])==ratios[0],
            'ratio_upper':F(primary['outer_ratio'])==ratios[1],
            'all_boxes':primary['counts']['boxes']==counts['all_span_admissible_boxes'],
            'lifted_atoms':primary['window_atoms_distinct']==len(lifted),
            'unique_inner':primary['inner_rectangles_unique']==len(inner),
            'unique_outer':primary['outer_rectangles_unique']==len(outer),
        }
        assert all(checks.values()),checks
        output['primary_comparison']={
            'path':str(args.compare_primary.resolve()),
            'sha256':hashlib.sha256(primary_raw).hexdigest(),'checks':checks}
    if args.compare_rectangles:
        rectangle_raw=args.compare_rectangles.read_bytes()
        primary_rectangles=json.loads(rectangle_raw)
        primary_scale=int(primary_rectangles['scale'])
        checks={}
        for key,ours in [('inner',inner),('outer',outer)]:
            theirs={tuple(F(int(x),primary_scale) for x in r)
                    for r in primary_rectangles[key]}
            our_rationals={tuple(F(x,grid['S']) for x in r) for r in ours}
            checks[key+'_complete_rectangle_set']=theirs==our_rationals
        assert all(checks.values()),checks
        output['primary_rectangle_comparison']={
            'path':str(args.compare_rectangles.resolve()),
            'sha256':hashlib.sha256(rectangle_raw).hexdigest(),'checks':checks}
    output['total_seconds']=time.perf_counter()-start
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ('counts','area_lower','area_upper',
                    'weak_ratio_lower','weak_ratio_upper','total_seconds')},indent=2))

if __name__=='__main__':
    main()
