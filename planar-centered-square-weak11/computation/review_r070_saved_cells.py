#!/usr/bin/env python3
"""Independent R70 saved-leaf/source/area verifier. No producer imports."""
from array import array
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import time


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def subsets(points, masses):
    """Literal subset recurrence, not coordinate-slab generation."""
    rows = [None] * (1 << len(points))
    for mask in range(1, len(rows)):
        bit = mask & -mask
        i = bit.bit_length() - 1
        x, y = points[i]
        prev = rows[mask ^ bit]
        if prev is None:
            bbox, mass = (x, x, y, y), masses[i]
        else:
            (a, b, c, d), mass = prev
            bbox, mass = (min(a, x), max(b, x), min(c, y), max(d, y)), mass + masses[i]
        rows[mask] = bbox, mass
    return rows


def closed(rows, points):
    result = {}
    for bbox, _ in rows[1:]:
        if bbox not in result:
            a, b, c, d = bbox
            result[bbox] = sum(1 << i for i, (x, y) in enumerate(points) if a <= x <= b and c <= y <= d)
    return {mask: (bbox, rows[mask][1]) for bbox, mask in result.items()}


def bounds2(cell, bbox):
    """Project bbox midpoints onto cell, then check corners for maxima."""
    near, far = [], []
    for lo, hi, a, b in ((cell[0], cell[1], bbox[0], bbox[1]), (cell[2], cell[3], bbox[2], bbox[3])):
        projected2 = min(2 * hi, max(2 * lo, a + b))
        near.append(max(abs(projected2 - 2 * a), abs(projected2 - 2 * b)))
        far.append(2 * max(abs(lo - a), abs(lo - b), abs(hi - a), abs(hi - b)))
    return max(near), max(far)


def passes(c2, r2, a, L):
    """Piecewise radical elimination, independently separated by active minimum."""
    if r2 * r2 > L:
        return False
    # Unconstrained critical mass. Here f_min^2=3(L-a)/4.
    if L - a >= 3 * a and L - a >= 3 * r2 * r2:
        return c2 * c2 > 3 * (L - a)
    # Capture-radius endpoint. Here f= sqrt(L+4r^2-a)-r.
    if r2 * r2 >= a:
        return (c2 + r2) ** 2 > 4 * (L + r2 * r2 - a)
    # Internal-mass endpoint. Here f=sqrt(L)-sqrt(a)/2.
    residual = 4 * L - a - c2 * c2
    return residual < 0 or 4 * c2 * c2 * a > residual * residual


def quadtree(leaves, domain):
    trie = {}
    states = {-1: 0, 0: 0, 1: 0}
    areas = {-1: 0, 0: 0, 1: 0}
    for row in leaves:
        target, cell, node = tuple(row['cell']), tuple(domain), trie
        for _ in range(row['depth']):
            assert 'leaf' not in node
            a, b, c, d = cell
            xm, ym = (a + b) // 2, (c + d) // 2
            children = [(a, xm, c, ym), (xm, b, c, ym), (a, xm, ym, d), (xm, b, ym, d)]
            matches = [j for j, cc in enumerate(children) if cc[0] <= target[0] < target[1] <= cc[1] and cc[2] <= target[2] < target[3] <= cc[3]]
            assert len(matches) == 1
            j = matches[0]
            cell, node = children[j], node.setdefault(j, {})
        assert cell == target and not node
        node['leaf'] = True
        state = row['state']
        states[state] += 1
        areas[state] += (cell[1] - cell[0]) * (cell[3] - cell[2])
    def complete(node):
        if 'leaf' in node:
            assert len(node) == 1
        else:
            assert set(node) == {0, 1, 2, 3}
            for child in node.values():
                complete(child)
    complete(trie)
    assert sum(areas.values()) == (domain[1] - domain[0]) * (domain[3] - domain[2])
    return states, areas


def halo_rectangles(rows, D):
    inner, outer = [], []
    # Every literal subset is considered; deduplicate coincident boxes afterward.
    for bbox, mass in rows[1:]:
        a, b, c, d = bbox
        if mass <= max(b - a, d - c) ** 2:
            continue
        aa, bb, cc, dd = [int(v * D) for v in bbox]
        assert all(v * D == int(v * D) for v in bbox)
        k = isqrt(mass.numerator * D * D // (4 * mass.denominator))
        for dst, radius in ((inner, max(0, k - 1)), (outer, k + 1)):
            rec = (bb - radius, aa + radius, dd - radius, cc + radius)
            if rec[0] < rec[1] and rec[2] < rec[3]:
                dst.append(rec)
    return sorted(set(inner)), sorted(set(outer))


def area_grid(layers):
    """Coordinate-compressed 2D difference arrays; no interval-strip merge."""
    xs = sorted({r[j] for layer in layers for r in layer for j in (0, 1)})
    ys = sorted({r[j] for layer in layers for r in layer for j in (2, 3)})
    nx, ny = len(xs), len(ys)
    xi, yi = {x: i for i, x in enumerate(xs)}, {y: i for i, y in enumerate(ys)}
    flags = bytearray((nx - 1) * (ny - 1))
    for bit, rects in enumerate(layers):
        diff = array('i', [0]) * (nx * ny)
        for a, b, c, d in rects:
            i, j, k, l = xi[a], xi[b], yi[c], yi[d]
            diff[k * nx + i] += 1
            diff[k * nx + j] -= 1
            diff[l * nx + i] -= 1
            diff[l * nx + j] += 1
        for y in range(ny - 1):
            offset, prior, foff, running = y * nx, (y - 1) * nx, y * (nx - 1), 0
            for x in range(nx - 1):
                running += diff[offset + x]
                value = running + (diff[prior + x] if y else 0)
                diff[offset + x] = value
                assert value >= 0
                if value:
                    flags[foff + x] |= 1 << bit
    totals = [0] * 10
    widths = [b - a for a, b in zip(xs, xs[1:])]
    for y in range(ny - 1):
        height = ys[y + 1] - ys[y]
        for x, width in enumerate(widths):
            v = flags[y * (nx - 1) + x]
            mi, mo, ii, io, ri, ro, fi, fo = [bool(v & (1 << j)) for j in range(8)]
            conditions = (mi, mo, ii, io, int(ii) + int(mi and not io), int(io) + int(mo and not ii), mi and not (io or ro), mo and not (ii or ri), fi and not ro, fo and not ri)
            area = width * height
            for j, keep in enumerate(conditions):
                if keep:
                    totals[j] += area * int(keep)
    return totals, (nx - 1) * (ny - 1)


def verify(root, job, spec):
    ip = root / job['input']
    assert sha(ip) == job['sha256']
    data = json.loads(ip.read_text())
    base = root / 'planar-centered-square-weak11/computation/results/round070-envelope'
    op = base / (job['id'] + '.json')
    saved = json.loads(op.read_text())
    cp = base / saved['certificate_file']
    assert sha(cp) == saved['certificate_sha256']
    assert saved['input_sha256'] == sha(ip)
    cert = json.loads(gzip.decompress(cp.read_bytes()))
    D, L = cert['scale'], cert['L_units']
    pts = [tuple(map(F, p)) for p in data['points']]
    masses = list(map(F, data['initial_masses']))
    allrows = subsets(pts, masses)
    eligible = [(i, mass) for i, rec in enumerate(allrows) if i for bbox, mass in [rec] if mass > max(bbox[1] - bbox[0], bbox[3] - bbox[2]) ** 2]
    maximum = max(mass for _, mass in eligible)
    assert [mask for mask, mass in eligible if mass == maximum] == saved['all_maximum_masks']
    assert maximum == F(saved['L']) == F(L, D * D)
    S = data['S_mask']
    indices = [i for i in range(len(pts)) if S >> i & 1]
    outside = [i for i in range(len(pts)) if not S >> i & 1]
    pp, mm = [pts[i] for i in indices], [masses[i] for i in indices]
    inside = subsets(pp, mm)
    closures = closed(inside, pp)
    expected = {}
    for mask, (bbox, mass) in closures.items():
        if mass < maximum:
            expected[mask] = (tuple(int(v * D) for v in bbox), int(mass * D * D))
    assert {lab['mask']: (tuple(lab['bbox']), lab['mass']) for lab in cert['labels']} == expected
    assert len(closures) == saved['internal_all_closed_bboxes']
    states, ar = quadtree(cert['leaves'], cert['domain'])
    assert {str(k): v for k, v in states.items()} == saved['cell_counts']
    checks = 0
    for row in cert['leaves']:
        if not row['state']:
            continue
        cmin, cmax = bounds2(row['cell'], cert['S_bbox'])
        answers = []
        for bbox, a in expected.values():
            rmin, rmax = bounds2(row['cell'], bbox)
            checks += 1
            if row['state'] == 1:
                if rmax * rmax < L and passes(cmin, rmax, a, L):
                    answers.append(True)
                    break
            else:
                assert rmin * rmin > L or not passes(cmax, rmin, a, L)
        if row['state'] == 1:
            assert answers
    ii, io = halo_rectangles(inside, D)
    ri, ro = halo_rectangles(subsets([pts[i] for i in outside], [masses[i] for i in outside]), D)
    fi, fo = ([], []) if job['id'] == 'historical16' else halo_rectangles(allrows, D)
    mi = [tuple(r['cell']) for r in cert['leaves'] if r['state'] == 1]
    mo = [tuple(r['cell']) for r in cert['leaves'] if r['state'] != -1]
    vals, cells = area_grid([mi, mo, ii, io, ri, ro, fi, fo])
    names = ('M_S_envelope', 'I', 'I_union_M_S_only_budget', 'M_minus_I_union_R_diagnostic', 'true_Delta')
    for j, name in enumerate(names):
        if name == 'true_Delta' and (job['id'] == 'historical16' or data.get('mode') == 'S_only'):
            continue
        assert [F(saved['regions'][name][q]) for q in ('lower', 'upper')] == [F(v, D * D) for v in vals[2 * j:2 * j + 2]], (job['id'], name)
    return {'id': job['id'], 'literal_subsets': len(allrows) - 1, 'proper_internal_closures': len(expected), 'saved_leaves': len(cert['leaves']), 'leaf_label_tests': checks, 'area_grid_cells': cells, 'exact_envelope_and_internal_union_area_match': True, 'true_Delta_independently_rebuilt': job['id'] == 'balanced8', 'historical_P_Delta_reused_only': job['id'] == 'historical16'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    start = time.perf_counter()
    root = args.archive_root.resolve()
    manifest_path = root / 'planar-centered-square-weak11/computation/results/round070-compute-integration-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    pins = {}
    for row in manifest['files'] + manifest['existing_dependencies']:
        p = root / row['destination']
        assert sha(p) == row['sha256'] and p.stat().st_size == row['bytes']
        pins[row['destination']] = row['sha256']
    specpath = root / 'planar-centered-square-weak11/computation/round070/three-pilot-spec.json'
    spec = json.loads(specpath.read_text())
    cases = [verify(root, job, spec) for job in spec['jobs']]
    output = {'status': 'PASS_INDEPENDENT_SAVED_R70_CELLS', 'source_sha256': sha(Path(__file__)), 'manifest_sha256': sha(manifest_path), 'pins': pins, 'python': platform.python_version(), 'cases': cases, 'elapsed_seconds': time.perf_counter() - start, 'scope': 'Literal subset source audit, complete quadtree coverage, all classified leaves via independent three-branch radical elimination, 2D difference-array areas. Undecided leaves remain outward; no producer imports. Historical true P/Delta only byte-bound, not recomputed.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'status': output['status'], 'seconds': output['elapsed_seconds'], 'cases': cases}))


if __name__ == '__main__':
    main()
