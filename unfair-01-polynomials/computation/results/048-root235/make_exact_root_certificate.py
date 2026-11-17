#!/usr/bin/env python3
"""Make rational root-disk and interval-negative certificates for frozen masks.

All certificate arithmetic after Newton seed refinement is fractions.Fraction.
For a disk |z-c|<=r, the exact Rouché test is
  max(|Re C'(c)|,|Im C'(c)|) r > |Re C(c)|+|Im C(c)|+T(c,r),
where T bounds the Taylor terms of degree >=2 using |c|<=|Re c|+|Im c|.
Thus each disk contains exactly one C-root; pairwise L-infinity separation
greater than 2r makes the disks disjoint.  Coefficient balls are then
propagated exactly through monic group factors and their constant normalization.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def pack(q: F) -> dict[str, str]:
    return {'n': str(q.numerator), 'd': str(q.denominator)}


def unpack(x: dict[str, str]) -> F:
    return F(int(x['n']), int(x['d']))


def cadd(x, y): return (x[0]+y[0], x[1]+y[1])
def cmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def cscale(x, q): return (x[0]*q, x[1]*q)
def cabs1(x): return abs(x[0])+abs(x[1])
def clinf(x): return max(abs(x[0]), abs(x[1]))


def horner(c, z):
    out = (F(0), F(0))
    for a in reversed(c):
        out = cadd(cmul(out, z), (F(a), F(0)))
    return out


def deriv_horner(c, z):
    out = (F(0), F(0))
    for k in range(len(c)-1, 0, -1):
        out = cadd(cmul(out, z), (F(k*c[k]), F(0)))
    return out


def tail_bound(c, z, r):
    """Exact bound for |C(z+w)-C(z)-C'(z)w| on |w|=r."""
    rho = cabs1(z)
    total = F(0)
    for i, ai in enumerate(c):
        if i >= 2 and ai:
            total += abs(ai) * ((rho+r)**i-rho**i-i*r*rho**(i-1))
    return total


def rouche_certificate(c, z, r):
    value, deriv = horner(c, z), deriv_horner(c, z)
    left = clinf(deriv)*r
    right = cabs1(value)+tail_bound(c, z, r)
    return {'left': left, 'right': right, 'margin': left-right, 'value': value, 'derivative': deriv}


def mp_poly(c, z):
    out = mp.mpc(0)
    for a in reversed(c): out = out*z+a
    return out


def mp_deriv(c, z):
    out = mp.mpc(0)
    for k in range(len(c)-1, 0, -1): out = out*z+k*c[k]
    return out


def refine(seed, c, tol, steps):
    z = mp.mpc(mp.mpf(str(seed[0])), mp.mpf(str(seed[1])))
    for count in range(1, steps+1):
        v, d = mp_poly(c, z), mp_deriv(c, z)
        if not mp.isfinite(v.real) or not mp.isfinite(v.imag) or not d:
            raise RuntimeError('Newton hit nonfinite/critical derivative')
        dz = v/d
        z -= dz
        if abs(dz) <= tol*max(1, abs(z)) and abs(mp_poly(c, z)) <= tol:
            return z, count
    raise RuntimeError('Newton iteration limit')


def rounded_fraction(x, places):
    scale = mp.mpf(10) ** places
    return F(int(mp.nint(x*scale)), 10**places)


class Ball:
    __slots__ = ('c', 'r')
    def __init__(self, c, r=F(0)):
        if r < 0: raise ValueError('negative radius')
        self.c, self.r = F(c), F(r)
    def add(self, other): return Ball(self.c+other.c, self.r+other.r)
    def mul(self, other):
        return Ball(self.c*other.c, abs(self.c)*other.r+abs(other.c)*self.r+self.r*other.r)
    def reciprocal(self):
        lo, hi = self.c-self.r, self.c+self.r
        if lo <= 0 <= hi: raise ValueError('normalizing constant interval contains zero')
        a, b = F(1, lo), F(1, hi)
        lo2, hi2 = min(a, b), max(a, b)
        return Ball((lo2+hi2)/2, (hi2-lo2)/2)
    def upper(self): return self.c+self.r


def poly_mul(left, right):
    out = [Ball(0) for _ in range(len(left)+len(right)-1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right): out[i+j] = out[i+j].add(a.mul(b))
    return out


def group_ball_factor(group, r):
    if len(group) == 1:
        a, b = group[0]
        if b != 0: raise ValueError('real group has nonreal rational centre')
        return [Ball(-a, r), Ball(1)]
    if len(group) == 2:
        a, b = group[0]
        a2, b2 = group[1]
        if a2 != a or b2 != -b: raise ValueError('pair centres are not exact conjugates')
        l1 = abs(a)+abs(b)
        return [Ball(a*a+b*b, 2*l1*r+r*r), Ball(-2*a, 2*r), Ball(1)]
    raise ValueError('group cardinality not one or two')


def normalised_product(group_factors, selected):
    raw = [Ball(1)]
    for bit, factor in zip(selected, group_factors):
        if bit: raw = poly_mul(raw, factor)
    inv = raw[0].reciprocal()
    out = [Ball(1)]
    out.extend(x.mul(inv) for x in raw[1:])
    return out, raw[0]


def bits(groups, mask):
    if not 0 <= mask < (1 << (groups-1)): raise ValueError('mask range')
    return [0]+[(mask >> (groups-1-g)) & 1 for g in range(1, groups)]


def group_centres(task, places, newton_tol, steps):
    c = [int(v) for v in task['C']]
    groups, count_log = [], []
    for raw in task['roots']:
        if len(raw) == 1:
            z, count = refine(raw[0], c, newton_tol, steps)
            a = rounded_fraction(z.real, places)
            groups.append([(a, F(0))]); count_log.append(count)
        elif len(raw) == 2:
            z, count = refine(raw[0], c, newton_tol, steps)
            a, b = rounded_fraction(z.real, places), rounded_fraction(z.imag, places)
            if b == 0: raise RuntimeError('complex root rounded to real centre')
            groups.append([(a, b), (a, -b)]); count_log.extend([count, count])
        else: raise RuntimeError('bad root group cardinality')
    return groups, count_log


def serialize_groups(groups):
    return [[{'re': pack(a), 'im': pack(b)} for a, b in group] for group in groups]


def certify_task(task, masks, r, places, newton_tol, steps):
    c = [int(v) for v in task['C']]
    centres, newton_steps = group_centres(task, places, newton_tol, steps)
    flat = [z for g in centres for z in g]
    root_rows = []
    for z in flat:
        row = rouche_certificate(c, z, r)
        if row['margin'] <= 0: raise RuntimeError('Rouche margin nonpositive')
        root_rows.append({'centre': {'re': pack(z[0]), 'im': pack(z[1])},
                          'left': pack(row['left']), 'right': pack(row['right']), 'margin': pack(row['margin'])})
    sep = min(max(abs(x[0]-y[0]), abs(x[1]-y[1])) for i, x in enumerate(flat) for y in flat[:i])
    if sep <= 2*r: raise RuntimeError('root disks not pairwise disjoint in L-infinity')
    factors = [group_ball_factor(g, r) for g in centres]
    masks_out = []
    for mask in masks:
        b = bits(task['groups'], mask)
        a_poly, a_const = normalised_product(factors, b)
        b_poly, b_const = normalised_product(factors, [1-x for x in b])
        candidates = [('A', i, q) for i, q in enumerate(a_poly) if q.upper() < 0]
        candidates += [('B', i, q) for i, q in enumerate(b_poly) if q.upper() < 0]
        if not candidates: raise RuntimeError(f'mask {mask} has no certified negative interval')
        side, index, ball = min(candidates, key=lambda z: z[2].upper())
        masks_out.append({'mask': mask, 'bits_by_group': b, 'negative_factor': side,
                          'negative_index': index, 'negative_upper': pack(ball.upper()),
                          'negative_center': pack(ball.c), 'negative_radius': pack(ball.r),
                          'raw_constant_A': {'center': pack(a_const.c), 'radius': pack(a_const.r)},
                          'raw_constant_B': {'center': pack(b_const.c), 'radius': pack(b_const.r)}})
    return {'task': task['id'], 'degree': len(c)-1, 'groups': task['groups'],
            'centres_by_group': serialize_groups(centres), 'newton_steps': newton_steps,
            'root_disks': root_rows, 'min_linf_center_separation': pack(sep),
            'mask_certificates': masks_out}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--diagnostics', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--tasks', default='18,19,20,21,22,23')
    p.add_argument('--limit-per-task', type=int)
    p.add_argument('--dps', type=int, default=100)
    p.add_argument('--places', type=int, default=80)
    p.add_argument('--radius-exp', type=int, default=60)
    p.add_argument('--newton-tol-exp', type=int, default=85)
    p.add_argument('--max-newton-steps', type=int, default=80)
    args = p.parse_args()
    if not (args.dps > args.newton_tol_exp+10 and args.places > args.radius_exp):
        p.error('insufficient precision/centre places for requested radius')
    mp.mp.dps = args.dps
    ids = {int(x) for x in args.tasks.split(',')}
    inp, diag = json.loads(args.input.read_text()), json.loads(args.diagnostics.read_text())
    tasks = {t['id']: t for t in inp['tasks'] if t['id'] in ids}
    chosen = {i: [] for i in ids}
    for row in diag:
        if row.get('task') in ids and row.get('classification') == 'nearunknown':
            chosen[row['task']].append(int(row['mask']))
    if args.limit_per_task is not None: chosen = {i: v[:args.limit_per_task] for i, v in chosen.items()}
    if set(tasks) != ids or any(not chosen[i] for i in ids): raise RuntimeError('selection mismatch')
    r, tol = F(1, 10**args.radius_exp), mp.mpf(10)**(-args.newton_tol_exp)
    certs = [certify_task(tasks[i], chosen[i], r, args.places, tol, args.max_newton_steps) for i in sorted(ids)]
    report = {'status': 'completed_exact_rational_interval_certificate',
              'proof_scope': 'only listed frozen masks; not a global exclusion',
              'input': str(args.input), 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'diagnostics': str(args.diagnostics), 'diagnostics_sha256': hashlib.sha256(args.diagnostics.read_bytes()).hexdigest(),
              'tasks': sorted(ids), 'selected_counts': {str(i): len(chosen[i]) for i in sorted(ids)},
              'radius': pack(r), 'decimal_centre_places': args.places,
              'newton_dps': args.dps, 'newton_tolerance_exponent': args.newton_tol_exp,
              'certificates': certs, 'mpmath': mp.__version__, 'python': sys.version,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'tasks': report['tasks'],
                      'masks': sum(report['selected_counts'].values()), 'out': str(args.out)}, ensure_ascii=False))


if __name__ == '__main__': main()
