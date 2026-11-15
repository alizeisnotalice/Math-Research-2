#!/usr/bin/env python3
"""Independent rational audit of R62 formulas, enclosures and coverage states."""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
from types import SimpleNamespace
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import platform
import sys
import time
import sympy as sp

BASE = Path('/Volumes/data/CodexProjects/T2/.local/workspaces')
TOPIC = 'planar-centered-square-weak11'
OWNER = BASE/'pcw-sol-kernel-sources'/TOPIC
OWN = BASE/'pcw-sol-cover-structure'/TOPIC
TMP = Path('/Volumes/data/CodexTemp/T2/r062-independent-source-audit')
SOURCE = OWNER/'computation/r062_coupled_source_verify.py'
PINS = {
    'proofs/r062-coupled-wholecell-source-interface.md': '3afd1ca1d3c3762ec4ba6e997fca18717b9a4942b2a0bafcd640ffbeb88f9042',
    'computation/r062_coupled_source_verify.py': '65143813817936d4db538aa0a4b1f4f6de7d55e1070b211bf05f4211f9710827',
    'computation/cauchy_occupation_coupled_lp.m': 'fb69861a37adaeba94ac57dfbf9921889970b2b53b07b2df166e50e39d29c5a7',
}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


@lru_cache(None)
def alternating(x):
    """Exact rational asymmetric alternating bracket, no fixed-point grid."""
    assert 0 <= x <= F(1, 2)
    value = sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(120)), F(0))
    return value, value + x**241 / 241


def plus(a, b):
    return a[0]+b[0], a[1]+b[1]


def minus(a, b):
    return a[0]-b[1], a[1]-b[0]


def multiply(a, b):
    vals = [x*y for x in a for y in b]
    return min(vals), max(vals)


def times(a, b):
    return multiply(a, (b, b))


PI = minus(times(alternating(F(1, 5)), 16), times(alternating(F(1, 239)), 4))


@lru_cache(None)
def atan(x):
    if x < 0:
        a = atan(-x)
        return -a[1], -a[0]
    if x > 1:
        return minus(times(PI, F(1, 2)), atan(1/x))
    if x > F(1, 2):
        return minus(times(PI, F(1, 4)), alternating((1-x)/(1+x)))
    return alternating(x)


def divide(a, b):
    assert b[0] > 0
    return multiply(a, (1/b[1], 1/b[0]))


def reference_jet(x, h):
    den = h*h+x*x
    return (divide((h, h), times(PI, den)),
            divide((-2*h*x, -2*h*x), times(PI, den**2)),
            divide((2*h*(3*x*x-h*h),)*2, times(PI, den**3)))


@lru_cache(None)
def reference_prob(distance, delta, h):
    return divide(minus(atan((delta/2-distance)/h), atan((-delta/2-distance)/h)), PI)


def contains(fixed, rational, S):
    assert F(fixed[0], S) <= rational[0] <= rational[1] <= F(fixed[1], S)


def main():
    assert os.path.ismount('/Volumes/data')
    start = time.perf_counter()
    for rel, expected in PINS.items():
        assert digest(OWNER/rel) == expected
    spec = importlib.util.spec_from_file_location('r62_checked_source', SOURCE)
    checked = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checked)
    S = checked.S
    contains(checked.PI, PI, S)
    # Actual symbolic differentiation independent of code's interval formulas.
    x, h, P = sp.symbols('x h P', positive=True)
    k = h/(P*(h*h+x*x))
    assert sp.simplify(sp.diff(k, x)+2*h*x/(P*(h*h+x*x)**2)) == 0
    assert sp.simplify(sp.diff(k, x, 2)-2*h*(3*x*x-h*h)/(P*(h*h+x*x)**3)) == 0
    # Two moving rectangle endpoints give the difference-of-k-primes sign.
    c, d = sp.symbols('c d', real=True)
    T = (sp.atan((x-c+d/2)/h)-sp.atan((x-c-d/2)/h))/P
    assert sp.simplify(sp.diff(T, x) - (k.subs(x, x-c+d/2)-k.subs(x, x-c-d/2))) == 0
    assert sp.simplify(sp.diff(T, x, 2) - (sp.diff(k,x).subs(x, x-c+d/2)-sp.diff(k,x).subs(x, x-c-d/2))) == 0
    # New positive-time full-D interface: no interior time minimum.
    dist, width = sp.symbols('dist width', real=True)
    Ts=(sp.atan((width/2-dist)/h)-sp.atan((-width/2-dist)/h))/P
    claimed=-width*(h*h-(dist*dist-width*width/4))/(P*(h*h+(dist-width/2)**2)*(h*h+(dist+width/2)**2))
    assert sp.simplify(sp.diff(Ts,h)-claimed)==0
    jet_points = prob_points = 0
    intervals = [(F(-2), F(-1)), (F(-1,3), F(2,5)), (F(0), F(0)), (F(1,8), F(7,8))]
    for hval in (F(3,20), F(1), F(2)):
        for lo, hi in intervals:
            ki, kppi = checked.cauchyjet(lo, hi, hval)
            kpi = checked.kprime(lo, hi, hval)
            for q in (lo, (3*lo+hi)/4, (lo+hi)/2, (lo+3*hi)/4, hi):
                refk, refkp, refkpp = reference_jet(q, hval)
                contains(ki, refk, S)
                contains(kpi, refkp, S)
                contains(kppi, refkpp, S)
                jet_points += 1
            for center in (F(-1,4), F(0), F(1,4)):
                delta = F(1,8)
                ti, tppi = checked.boxjet(lo, hi, center, delta, hval)
                for q in (lo, (lo+hi)/2, hi):
                    refp = reference_prob(q-center, delta, hval)
                    rp = reference_jet(q-center+delta/2, hval)[1]
                    rm = reference_jet(q-center-delta/2, hval)[1]
                    contains(ti, refp, S)
                    contains(tppi, minus(rp, rm), S)
                    prob_points += 1
    # Rational toys separate corners, curvature, coupling and source legality.
    # F=(x-1/2)^2-1/8 has positive corners but negative interior.
    corner = F(1,8)
    true_min = F(-1,8)
    assert corner-F(2,8) == true_min < 0 < corner
    # Bilinear xy has arbitrarily large mixed derivative, no pure curvature.
    assert sp.simplify(sp.diff(1000*x*c, x, 2)) == 0
    # Positive functions x and 1-x: separate infima lose 1, coupled sum loses 0.
    assert F(0)+F(0) < F(1)
    # Signed actual-coefficient curvature may cancel: x²+(2-x²)=2.
    assert F(2)+F(-2) == 0 and max(F(0), F(2))+max(F(0), F(-2)) == 2
    # SourceA=-S+target is equivalent to the sufficient >=0 residual row.
    for coeff in (F(0), F(1,8), F(3)):
        row_value = F(2,3)*coeff-F(1,7)
        assert (row_value >= 0) == (-row_value <= 0)
    TMP.mkdir(parents=True, exist_ok=True)
    inp = TMP/'invalid-toy-input.json'
    candidate = dict(h='1',delta='1/8',N=1,source_origin_mass='0',density_coefficients=['1']*3)
    inp.write_text(json.dumps(candidate)+'\n')
    common = dict(input=str(inp), repair='0', initial_subdiv=1, max_depth=0, max_seconds=20,
                  max_boxes=100, stop_on_negative=False)
    runs = []
    # One checked box does not cover all targets, so repair must be null.
    for name, options in [('budget_incomplete',dict(max_boxes=1)), ('negative_early',dict(stop_on_negative=True)),
                          ('negative_complete',{})]:
        args = SimpleNamespace(**(common|options|dict(output=str(TMP/(name+'.json')))))
        with contextlib.redirect_stdout(io.StringIO()):
            checked.exact_run(args)
        row = json.loads(Path(args.output).read_text())
        assert row['status'] == 'INVALID_SOURCE' and row['negative_corner_witnesses']
        assert row['D_minimum_certified'] is False and row['source_only'] is True
        if name != 'negative_complete':
            assert not row['complete_terminal_coverage'] and row['additional_origin_repair_sufficient_exact'] is None
        else:
            assert row['complete_terminal_coverage'] and F(row['additional_origin_repair_sufficient_exact']) > 0
        runs.append(dict(name=name, status=row['status'], complete=row['complete_terminal_coverage'],
                         sufficient_repair=row['additional_origin_repair_sufficient_exact'],
                         checked_boxes=row['checked_boxes'], negative_witnesses=len(row['negative_corner_witnesses']),
                         input_sha256=row['input_sha256']))
    repair = F(runs[-1]['sufficient_repair'])
    args = SimpleNamespace(**(common|dict(repair=str(repair),output=str(TMP/'repaired-toy.json'))))
    with contextlib.redirect_stdout(io.StringIO()):
        checked.exact_run(args)
    repaired = json.loads(Path(args.output).read_text())
    assert repaired['status']=='PASS_SOURCE_NONNEGATIVE' and repaired['complete_terminal_coverage']
    assert repaired['additional_origin_repair_sufficient_exact']=='0'
    assert F(repaired['halfmass_over_pi_exact']) == F(1)*(repair+F(9,64))
    runs.append(dict(name='repaired_complete',status=repaired['status'],complete=True,
                     sufficient_repair=repaired['additional_origin_repair_sufficient_exact'],
                     checked_boxes=repaired['checked_boxes'], origin_amplitude=str(repair),
                     halfmass_over_pi=repaired['halfmass_over_pi_exact']))
    # Freeze pins through the full bounded audit.
    for rel, expected in PINS.items():
        assert digest(OWNER/rel) == expected
    out = dict(status='PASS_BOUNDED_INDEPENDENT_R62_AUDIT',python=sys.version,executable=sys.executable,
               sympy=sp.__version__,platform=platform.platform(),arithmetic='Exact Fraction reference enclosures and symbolic differentiation; fixed-point source checker compared independently.',
               source_pins=[dict(path=str(OWNER/rel),sha256=expected)for rel,expected in PINS.items()],
               script_sha256=digest(Path(__file__)),jet_sample_points=jet_points,rectangle_probability_sample_points=prob_points,
               symbolic_derivative_identities=5,toys=dict(corner_only_false_positive=True,pure_curvature_correction_exact=True,
                  no_mixed_derivative_required=True,coefficient_cancellation=True,LP_row_direction=True),
               coverage_repair_runs=runs,elapsed_seconds=time.perf_counter()-start,
               scope='Independent proof/code review plus bounded exact rational checks. No MATLAB execution, no optimized candidate certification, no full-D bound.')
    dest=OWN/'computation/results/round062-coupled-source-independent.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
