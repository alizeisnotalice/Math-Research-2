#!/usr/bin/env python3
"""Independent exact event-radius audit for R22; standard library, one process."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import time


def atom_mass(atoms, x, r):
    return sum((m for p, m in atoms if max(abs(x[0] - p), abs(x[1])) <= r), F(0))


def atom_max(atoms, x):
    radii = sorted({max(abs(x[0] - p), abs(x[1])) for p, m in atoms})
    assert radii and radii[0] > 0
    return max(atom_mass(atoms, x, r) / (4 * r * r) for r in radii)


def overlap(center, radius, lo, hi):
    return max(F(0), min(center + radius, hi) - max(center - radius, lo))


def square_mass(atoms, eta, x, r):
    return sum((m * overlap(x[0], r, p - eta, p + eta)
                * overlap(x[1], r, -eta, eta) / (4 * eta * eta)
                for p, m in atoms), F(0))


def square_max(atoms, eta, x):
    """Exact all-r maximum for sums of constant square densities.

    Between edge-distance events the mass is A*r^2+B*r+C.
    Test interval endpoints and the only possible interior stationary radius
    r=-2*C/B, plus the r->0 and r->infinity limits.
    """
    events = {F(0)}
    for p, m in atoms:
        events.update(abs(x[0] - a) for a in (p - eta, p + eta))
        events.update(abs(x[1] - a) for a in (-eta, eta))
    knots = sorted(events)
    values = [F(0)]
    for r in knots[1:]:
        values.append(square_mass(atoms, eta, x, r) / (4 * r * r))
    for j, left in enumerate(knots):
        right = knots[j + 1] if j + 1 < len(knots) else None
        step = (right - left) / 3 if right is not None else F(1)
        a, b = left + step, left + 2 * step
        A = B = C = F(0)
        for p, m in atoms:
            lx = [overlap(x[0], r, p - eta, p + eta) for r in (a, b)]
            ly = [overlap(x[1], r, -eta, eta) for r in (a, b)]
            sx, sy = (lx[1] - lx[0]) / (b - a), (ly[1] - ly[0]) / (b - a)
            cx, cy = lx[0] - sx * a, ly[0] - sy * a
            density = m / (4 * eta * eta)
            A += density * sx * sy
            B += density * (sx * cy + sy * cx)
            C += density * cx * cy
        assert A * a * a + B * a + C == square_mass(atoms, eta, x, a)
        if left == 0:
            assert B == C == 0
            values.append(A / 4)
        if B:
            r = -2 * C / B
            if r > left and (right is None or r < right):
                values.append((A * r * r + B * r + C) / (4 * r * r))
    return max(values)


def construction(n, epsilon):
    points = [F(2**i) for i in range(n)]
    masses = [F(4)] + [4 * (points[i]**2 - points[i - 1]**2) for i in range(1, n)]
    return list(zip(points, masses)), epsilon / (144 * points[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source_topic', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    source_paths = [args.source_topic / p for p in (
        'proofs/cross-scale-insertion-obstruction.md',
        'computation/cross_scale_insertion_obstruction.py',
        'computation/results/cross-scale-budget/exact-insertion.json',
        'rounds/022.md')]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    source = json.loads(source_paths[2].read_text())
    source_rows = 0
    for case in source['cases']:
        n, epsilon = case['n'], F(case['epsilon'])
        atoms, eta = construction(n, epsilon)
        assert case['base_points'] == [[str(p), '0'] for p, m in atoms]
        assert list(map(F, case['base_masses'])) == [m for p, m in atoms]
        assert F(case['base_total']) == sum(m for p, m in atoms) == 4**n
        assert F(case['mollifier_square_halfside']) == eta
        assert F(case['disjoint_new_area_lower']) == n * epsilon / 9
        assert F(case['l1_new_area_lower']) == n * epsilon / 36
        for i, row in enumerate(case['rows']):
            p = atoms[i][0]
            w, r = epsilon / (9 * p), p + epsilon / (9 * p)
            mass = sum(m for q, m in atoms[:i + 1])
            assert mass == 4 * p * p == F(row['prefix_mass'])
            assert mass + epsilon - 4 * r * r == F(row['strict_slack']) >= 5 * epsilon / 81
            assert F(row['radius']) == r and F(row['width']) == w
            assert eta <= w / 16 and p + 3*w/4 + eta < r and 7*p/8 + eta < r
            assert 4 * (5*p/8 - eta)**2 > epsilon
            source_rows += 1

    atom_centers = l1_centers = baseline_boundaries = 0
    worst = {'base_atom': F(0), 'base_l1': F(0), 'new_atom': F(0), 'new_l1': F(0)}
    cases = []
    for n in (1, 2, 4, 6):
        for epsilon in (F(1), F(1, 2), F(1, 1024)):
            atoms, eta = construction(n, epsilon)
            new = [(F(0), epsilon)]
            for p, mass in atoms:
                w, r = epsilon/(9*p), p + epsilon/(9*p)
                prefix = sum(m for q, m in atoms if q <= p)
                # Includes the closed boundary of the proposed open atomic strips.
                for u in (-w, -w/2, F(0)):
                    for v in (p/2, 3*p/4, p, -p/2, -3*p/4, -p):
                        x = (u, v)
                        b, g = atom_max(atoms, x), atom_max(new, x)
                        assert b <= 1 and g <= 1
                        assert atom_mass(atoms + new, x, r) == prefix + epsilon > 4*r*r
                        worst['base_atom'] = max(worst['base_atom'], b)
                        worst['new_atom'] = max(worst['new_atom'], g)
                        atom_centers += 1
                # Four corners and one center in each closed contracted strip.
                probes = [(u, v) for u in (-3*w/4, -w/4) for v in (5*p/8, 7*p/8)]
                probes.append((-w/2, 3*p/4))
                for u, positive_v in probes:
                    for sign in (-1, 1):
                        x = (u, sign * positive_v)
                        b, g = square_max(atoms, eta, x), square_max(new, eta, x)
                        assert b <= 1 and g < 1
                        assert square_mass(atoms + new, eta, x, r) >= prefix + epsilon > 4*r*r
                        worst['base_l1'] = max(worst['base_l1'], b)
                        worst['new_l1'] = max(worst['new_l1'], g)
                        l1_centers += 1
                # Test the exact convolution half-plane boundary u=-eta.
                for v in (F(0), p/2, p, -p):
                    assert square_max(atoms, eta, (-eta, v)) <= 1
                    baseline_boundaries += 1
            cases.append({'n': n, 'epsilon': str(epsilon), 'eta': str(eta)})
    result = {
        'status': 'passed', 'python': platform.python_version(),
        'method': 'Fraction arithmetic; complete event-radius maxima for atomic and square densities',
        'source_sha256': hashes, 'source_rows_rechecked': source_rows,
        'small_cases': cases, 'atom_strip_centers': atom_centers,
        'l1_strip_centers': l1_centers, 'l1_halfplane_boundary_centers': baseline_boundaries,
        'largest_individual_halo_maxima_on_strips': {k: str(v) for k, v in worst.items()},
        'random_seed': None, 'tolerance': None,
        'elapsed_seconds': time.perf_counter() - started,
        'limits': 'Finite exact checks complement the all-N analytic review; no complete halo area or C2 computed.'}
    assert os.path.ismount('/Volumes/data'), 'data volume is not mounted'
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'python', 'source_rows_rechecked',
        'atom_strip_centers', 'l1_strip_centers', 'l1_halfplane_boundary_centers', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
