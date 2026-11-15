#!/usr/bin/env python3
"""Independent R61 audit: coordinate boxes, local overlays, direct cut sums.

No imports from, or execution of, the author's code.  All geometry uses scaled
integer endpoints.  Fraction is used for final ratios.  Zero masses/collisions
are normalized explicitly in the additional edge tests, not silently accepted
as the author's saved-input domain.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

B = F(1994727599733432151936139, 604462909807314587353088)
ROOT = Path('/Volumes/data/CodexProjects/T2/.local/workspaces')
TOPIC = 'planar-centered-square-weak11'
AUTHOR = ROOT / 'pcw-astra-recursive' / TOPIC
MANIFEST = AUTHOR / 'computation/results/round061-integration-manifest.json'
MANIFEST_SHA = 'dfe8b86abd2114051633b70c12e791981a31fc5b5cbaa7ea1d7fd52c8569100b'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def floor_root(q):
    """Independent integer root by bisection, with squared postcondition."""
    lo, hi = 0, 1 << ((q.bit_length() + 1) // 2)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid * mid <= q:
            lo = mid
        else:
            hi = mid
    assert lo * lo <= q < (lo + 1) * (lo + 1)
    return lo


def normalize(points, masses):
    d = defaultdict(int)
    for p, m in zip(points, masses):
        assert m >= 0
        d[tuple(p)] += m
    rows = sorted((p, m) for p, m in d.items() if m)
    return [p for p, _ in rows], [m for _, m in rows]


def boxes(points, masses, bits):
    """Enumerate coordinate-bound boxes, retain only tight nonempty boxes.

    Author enumerates all source subsets, then closes them.  Here no source
    subsets are enumerated: every coordinate interval pair is visited once.
    """
    assert bits >= 4 and len(points) == len(masses) and masses
    assert all(type(m) is int and m > 0 for m in masses)
    assert len(set(map(tuple, points))) == len(points)
    D = 1 << bits
    xs = sorted({p[0] for p in points})
    ys = sorted({p[1] for p in points})
    closed, rects, gamma = set(), [], 0
    for a in xs:
        for b in xs:
            if b < a:
                continue
            for c in ys:
                for d in ys:
                    if d < c:
                        continue
                    ids = [i for i, (x, y) in enumerate(points)
                           if a <= x <= b and c <= y <= d]
                    if not ids:
                        continue
                    if (min(points[i][0] for i in ids), max(points[i][0] for i in ids),
                        min(points[i][1] for i in ids), max(points[i][1] for i in ids)) != (a, b, c, d):
                        continue
                    mask = sum(1 << i for i in ids)
                    assert mask not in closed
                    closed.add(mask)
                    q = sum(masses[i] for i in ids)
                    if max(b - a, d - c) ** 2 >= q:
                        continue  # Exact true-empty or zero-width receiver test.
                    k = floor_root(q * D * D // 64)
                    areas = []
                    for kind, radius in ((0, k - 1), (1, k + 1)):
                        assert radius >= 0
                        # kind0 inner, kind1 outer. These are rational bounds.
                        assert (64 * radius * radius < q * D * D) if kind == 0 else (64 * radius * radius > q * D * D)
                        x0, x1 = b * D // 4 - radius, a * D // 4 + radius
                        y0, y1 = d * D // 4 - radius, c * D // 4 + radius
                        area = max(0, x1 - x0) * max(0, y1 - y0)
                        areas.append(area)
                        if area:
                            rects.append((x0, x1, y0, y1, mask, kind))
                    assert areas[1] >= areas[0]
                    gamma += areas[1] - areas[0]
    return rects, len(closed), gamma


def overlay(rects):
    """Horizontal strips, rebuilt local x-event overlays using active sets.

    No global y arrangement, no per-source overlap counters, and no author
    geometry helper.  Pair labels retain outer incidences ONLY for a lower
    bound on D(S), never as legal allocation edges.
    """
    ycoords = sorted({r[j] for r in rects for j in (2, 3)})
    weights = defaultdict(int)
    pairs = defaultdict(int)
    cells = 0
    for low, high in zip(ycoords, ycoords[1:]):
        active = [r for r in rects if r[2] <= low and high <= r[3]]
        events = defaultdict(list)
        for rid, r in enumerate(active):
            events[r[0]].append((rid, 1))
            events[r[1]].append((rid, -1))
        xcoords = sorted(events)
        live = set()
        for left, right in zip(xcoords, xcoords[1:]):
            for rid, direction in events[left]:
                if direction == 1:
                    live.add(rid)
                else:
                    live.remove(rid)
            inside, outside = 0, 0
            for rid in live:
                mask, kind = active[rid][4:6]
                if kind == 0:
                    inside |= mask
                else:
                    outside |= mask
            assert inside & ~outside == 0
            if outside:
                area = (right - left) * (high - low)
                assert area > 0
                weights[inside] += area
                pairs[inside, outside] += area
                cells += 1
    return dict(weights), dict(pairs), cells


def cut_audit(points, masses, bits, weights, pairs, gamma, saved=None):
    n = len(masses)
    raw, geom, low_cost = F(0), F(0), F(0)
    rm, gm = 0, 0
    compact, noncompact, improved = 0, 0, 0
    z_at_raw, m_at_raw = 0, 0
    D2 = 1 << (2 * bits)
    for s in range(1, 1 << n):
        ids = [j for j in range(n) if s & (1 << j)]
        m = sum(masses[j] for j in ids)
        # Direct sums over labels, intentionally not a zeta transform.
        z = sum(w for t, w in weights.items() if t & ~s == 0)
        lower = sum(w for (inner, outer), w in pairs.items()
                    if inner and outer & ~s == 0)
        assert 0 <= z - lower <= gamma
        ratio, lratio = F(16 * z, D2 * m), F(16 * lower, D2 * m)
        low_cost = max(low_cost, lratio)
        span = max(max(points[j][d] for j in ids) - min(points[j][d] for j in ids)
                   for d in (0, 1))
        compact += span * span <= m
        noncompact += span * span > m
        terminal = B * max(F(1), F(span * span, m))
        improved += terminal < ratio
        bound = min(ratio, terminal)
        if ratio > raw:
            raw, rm, z_at_raw, m_at_raw = ratio, s, z, m
        if bound > geom:
            geom, gm = bound, s
    assert raw - low_cost <= F(16 * gamma, D2 * min(masses))
    if saved is not None:
        assert F(saved['raw_bound']) == raw
        assert F(saved['geometric_bound']) == geom
        assert (saved['raw_critical_mask'], saved['geometric_critical_mask']) == (rm, gm)
        assert (saved['compact_cuts'], saved['noncompact_cuts'], saved['terminal_improved_cuts']) == (compact, noncompact, improved)
        assert saved['all_cuts'] == (1 << n) - 1
        assert saved['certifies_B'] == (geom <= B)
    return dict(raw=str(raw), geometric=str(geom), raw_mask=rm, geometric_mask=gm,
                lower_cost=str(low_cost), compact=compact, noncompact=noncompact,
                improved=improved, z_at_raw=z_at_raw, m_at_raw=m_at_raw,
                gamma=str(F(gamma, D2)), error=str(F(16 * gamma, D2 * min(masses))))


def flow_audit(flow, weights, masses, cut):
    scale = flow['integer_scale']
    cap = int(flow['source_capacity_multiplier'])
    assert scale == cut['m_at_raw'] and cap == cut['z_at_raw']
    assert flow['area_denominator'] == '2^(2*bits)'
    loads = [0] * len(masses)
    demand = defaultdict(int)
    seen = set()
    for i, t, value in flow['arcs']:
        value = int(value)
        assert 0 <= i < len(masses) and t in weights and t > 0
        assert t & (1 << i) and value > 0 and (i, t) not in seen
        seen.add((i, t))
        loads[i] += value
        demand[t] += value
    assert dict(demand) == {t: w * scale for t, w in weights.items() if t}
    assert loads == list(map(int, flow['source_load_numerators']))
    assert all(v <= cap * m for v, m in zip(loads, masses))
    assert sum(loads) == int(flow['total_routed_numerator'])
    # Units: arc area = value/(scale*D²), source budget=raw*m_i/16.
    return len(seen)


def edge_tests():
    tests = []
    # Perfect-square masses, singleton/collinear/2D boxes, strict equality,
    # near-empty receiver rectangles, repeated source labels, and zero faces.
    raw = [
        ('singleton', [(0, 0)], [16]),
        ('strict_equal_bbox', [(0, 0), (4, 0)], [8, 8]),
        ('just_nonempty_bbox', [(0, 0), (4, 0)], [8, 9]),
        ('degenerate_collinear', [(0, 0), (1, 0), (4, 0)], [1, 15, 1]),
        ('cross_sources', [(-2, 0), (2, 0), (0, -2), (0, 2)], [16] * 4),
        ('collision_zero_normalized', [(0, 0), (0, 0), (1, 1), (4, 0)], [7, 9, 0, 16]),
    ]
    # All mass assignments {1,16} on the fixed four-point rectangle.
    for mask in range(16):
        raw.append((f'rectangle_masses_{mask}', [(0, 0), (4, 0), (0, 3), (4, 3)],
                    [16 if mask & (1 << i) else 1 for i in range(4)]))
    for name, ps, ms in raw:
        points, masses = normalize(ps, ms)
        by_bits = []
        for bits in (4, 8, 24):
            rects, closed, gamma = boxes(points, masses, bits)
            w, pairs, cells = overlay(rects)
            row = cut_audit(points, masses, bits, w, pairs, gamma)
            row.update(bits=bits, cells=cells, boxes=closed)
            by_bits.append(row)
        # Dyadic enclosure refinement decreases gamma and contains all costs.
        assert all(F(a['gamma']) >= F(b['gamma']) for a, b in zip(by_bits, by_bits[1:]))
        tests.append(dict(name=name, normalized_points=points, normalized_mass=masses, runs=by_bits))
    return tests


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert os.path.ismount('/Volumes/data')
    started = time.perf_counter()
    assert sha(MANIFEST) == MANIFEST_SHA
    manifest = json.loads(MANIFEST.read_text())
    pins = []
    for row in manifest['files'] + manifest['existing_mathematical_sources']:
        p = Path(row['source'])
        assert sha(p) == row['sha256']
        if 'bytes' in row:
            assert p.stat().st_size == row['bytes']
        pins.append(dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size))
    pinmap = {Path(r['source']).name: Path(r['source']) for r in manifest['files']}
    inputs = json.loads(pinmap['r061-bounded-inputs.json'].read_text())
    summary = json.loads(pinmap['finite-source-cut-summary.json'].read_text())
    srows = {(r['batch'], r['case']): r for r in summary['rows']}
    receipt = json.loads(pinmap['run-receipt.json'].read_text())
    assert sorted(r['batch'] for r in receipt) == list(range(12))
    assert all(r['exit_status'] == 0 for r in receipt)
    rows, identities = [], Counter()
    total_arcs = total_cells = 0
    for batch in range(12):
        data = json.loads(pinmap[f'batch-{batch}.json'].read_text())
        assert data['status'] == 'PASS' and data['batch'] == batch and F(data['B']) == B
        assert data['script_sha256'] == sha(pinmap['finite_source_cut_certificate.py'])
        assert data['input_sha256'] == sha(pinmap['r061-bounded-inputs.json'])
        assert len(data['cases']) == len(inputs['batches'][batch]['inputs']) == 64
        for index, saved in enumerate(data['cases']):
            inp = saved['input']
            assert inp == inputs['batches'][batch]['inputs'][index]
            assert (inp['coordinate_unit'], inp['mass_unit']) == ('1/4', '1/16')
            points, masses = inp['points_quarters'], inp['mass_sixteenths']
            rects, closed, gamma = boxes(points, masses, saved['bits'])
            assert saved['closed_source_masks'] == closed
            assert saved['inner_rectangles'] == sum(r[5] == 0 for r in rects)
            assert saved['outer_rectangles'] == sum(r[5] == 1 for r in rects)
            w, pairs, cells = overlay(rects)
            assert w == {t: int(v) for t, v in saved['area_weight_numerators']}
            assert sum(w.values()) == int(saved['outer_area_numerator'])
            assert w.get(0, 0) == int(saved['uncertain_fringe_numerator'])
            cut = cut_audit(points, masses, saved['bits'], w, pairs, gamma, saved)
            arcs = flow_audit(saved['flow'], w, masses, cut)
            sr = srows[batch, index]
            assert sr['label_uncertainty_area'] == cut['gamma']
            assert sr['transport_cost_error_bound'] == cut['error']
            s = cut['geometric_mask']
            ids = [i for i in range(len(masses)) if s & (1 << i)]
            span = max(max(points[i][d] for i in ids) - min(points[i][d] for i in ids) for d in (0, 1))
            mass = sum(masses[i] for i in ids)
            proper, noncompact = s != (1 << len(masses)) - 1, span * span > mass
            assert sr['critical_span2_over_mass'] == str(F(span * span, mass))
            assert sr['proper_critical'] == proper and sr['critical_compact'] == (not noncompact)
            identity = tuple(sorted(tuple(p) + (m,) for p, m in zip(points, masses)))
            identities[identity] += 1
            rows.append(dict(batch=batch, case=index, n=len(masses), cells=cells, flow_arcs=arcs,
                             proper_critical=proper, noncompact_critical=noncompact, **cut))
            total_arcs += arcs
            total_cells += cells
        assert data['max_bound'] == str(max(F(r['geometric_bound']) for r in data['cases']))
        assert data['total_cuts'] == sum(r['all_cuts'] for r in data['cases'])
        assert data['all_certify_B'] == all(r['certifies_B'] for r in data['cases'])
        print(json.dumps(dict(batch=batch, audited=64, elapsed=time.perf_counter() - started)), flush=True)
    stats = dict(cases=len(rows), distinct_inputs=len(identities), repeated_controls=len(rows)-len(identities),
                 cuts=sum((1 << r['n']) - 1 for r in rows), compact=sum(r['compact'] for r in rows),
                 noncompact=sum(r['noncompact'] for r in rows), improved=sum(r['improved'] for r in rows),
                 proper_critical=sum(r['proper_critical'] for r in rows),
                 noncompact_critical=sum(r['noncompact_critical'] for r in rows),
                 reconstructed_cells=total_cells, verified_flow_arcs=total_arcs)
    assert (stats['cases'], stats['distinct_inputs'], stats['cuts'], stats['compact'], stats['noncompact'],
            stats['proper_critical'], stats['noncompact_critical'], stats['improved']) == (768, 725, 1217024, 397434, 819590, 430, 165, 0)
    best = max(rows, key=lambda r: F(r['geometric']))
    assert (best['batch'], best['case'], best['geometric'], best['geometric_mask']) == (7, 46, '611497158639875/525566558076928', 1791)
    tests = edge_tests()
    payload = dict(status='PASS', scope='All 768 saved R61 geometries/cuts/flows plus 66 exact edge-test runs; no general C2 assertion',
                   arithmetic='integer/Fraction; no floating geometric decisions; no tolerance; no randomization',
                   python=sys.version, executable=sys.executable, platform=platform.platform(),
                   script_sha256=sha(Path(__file__)), manifest_path=str(MANIFEST), manifest_sha256=sha(MANIFEST),
                   source_pins=pins, statistics=stats, maximum_case=best, edge_tests=tests,
                   wall_seconds=time.perf_counter()-started, rows=rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + '\n')
    print(json.dumps(dict(status='PASS', **stats, wall_seconds=payload['wall_seconds'])), flush=True)


if __name__ == '__main__':
    main()
