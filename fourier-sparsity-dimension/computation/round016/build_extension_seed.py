#!/usr/bin/env python3
"""Exhaustive four-variable bent classification and a full-selector regression.

All Boolean phases are represented by complete truth masks. Walsh tests use
exact Hamming counts against all characters, independently of the GPU FWHT.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time


def basis(vs):
    p = {}
    for x in vs:
        for j in sorted(p, reverse=True):
            if x >> j & 1:
                x ^= p[j]
        if x:
            j = x.bit_length()-1
            for k in p:
                if p[k] >> j & 1:
                    p[k] ^= x
            p[j] = x
    return tuple(p[j] for j in sorted(p, reverse=True))


def span(b):
    out = [0]
    for g in b:
        out += [x ^ g for x in out]
    return out


def capacity(b):
    code = span(b); dp = {0: 0}
    for s in sorted(code, key=int.bit_count):
        if s:
            dp[s] = max(sum(g & ~a == 0 for g in code).bit_length() + dp[s ^ a]
                        for a in code if a & (s & -s) and a & ~s == 0)
    return dp[65535]


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--out', type=Path, required=True)
    a.add_argument('--report', type=Path, required=True)
    args = a.parse_args(); began = time.monotonic()
    chars = [sum(((u & x).bit_count() & 1) << x for x in range(16)) for u in range(16)]
    bent = {f for f in range(1 << 16) if all(abs(16-2*(f ^ c).bit_count()) == 4 for c in chars)}
    assert len(bent) == 896
    f0 = sum(((((x >> 0) & 1)*((x >> 2) & 1)) ^ (((x >> 1) & 1)*((x >> 3) & 1))) << x for x in range(16))
    affine = basis([65535, *chars])
    assert len(affine) == 5
    acode = span(affine)
    allowed = {f ^ f0 for f in bent}
    assert all((g ^ c) in allowed for g in allowed for c in acode)

    def reduce(x):
        for v in affine:
            if x >> (v.bit_length()-1) & 1:
                x ^= v
        return x

    quotient = {reduce(h) for h in allowed}
    assert len(quotient) == 28
    levels = [{(0,)}]
    while True:
        nxt = set()
        for space in levels[-1]:
            for g in quotient - set(space):
                enlarged = tuple(sorted(set(space) | {v ^ g for v in space}))
                if set(enlarged) <= quotient:
                    nxt.add(enlarged)
        if not nxt:
            break
        levels.append(nxt)
    maxbases = [list(basis([*affine, *space])) for space in sorted(levels[-1])]
    maxdim = 5 + len(levels)-1
    assert all(len(b) == maxdim and all((f0 ^ h) in bent for h in span(b)) for b in maxbases)
    maximal_counts = {}
    for d, spaces in enumerate(levels):
        nxt = levels[d+1] if d+1 < len(levels) else set()
        maximal_counts[str(5+d)] = sum(not any(set(space) < set(upper) for upper in nxt) for space in spaces)
    assert maximal_counts == {'5': 0, '6': 0, '7': 45}
    max_capacities = [capacity(b) for b in maxbases]
    assert max_capacities == [8]*45

    # Independently enumerate all 65536 possible extra directions of the
    # four-variable selector built from the maximal two-variable affine code.
    old_c = [15, 10, 12]
    def embed(g, w):
        return sum(((g >> x) & 1) << (x+4*z+8*w) for x in range(4) for z in range(2))
    zmask = sum(1 << p for p in range(16) if p >> 2 & 1)
    sb = basis([*(embed(g, w) for w in range(2) for g in old_c), zmask])
    scode = set(span(sb)); assert len(sb) == 7
    base = sum(((((x >> 0) & 1)*((x >> 1) & 1)) ^ (((x >> 2) & 1)*((x >> 3) & 1))) << x for x in range(16))
    actual = {g for g in range(1 << 16) if all((base ^ g ^ h) in bent for h in scode)}
    assert actual == scode
    # Boundary: restoring global affine phases alone is insufficient.
    xz = sum((((x >> 0) & 1)*((x >> 2) & 1)) << x for x in range(16))
    xw = sum((((x >> 0) & 1)*((x >> 3) & 1)) << x for x in range(16))
    coupled = affine
    extension = basis([*coupled, xz, xw])
    assert xz not in scode and len(extension) == 7
    assert all((base ^ h) in bent for h in span(extension))
    yw = sum((((x >> 1) & 1)*((x >> 3) & 1)) << x for x in range(16))
    assert yw in scode and (base ^ xz ^ yw) not in bent
    yz = sum((((x >> 1) & 1)*((x >> 2) & 1)) << x for x in range(16))
    rank2_extension = basis([*affine, xz, yz])
    char3 = [sum(((u & x).bit_count() & 1) << x for x in range(8)) for u in range(8)]
    branch_checks = 0
    for v in span(rank2_extension):
        assert (base ^ v) in bent
        supports = []
        for w in range(2):
            f = sum((((base ^ v) >> (x+8*w)) & 1) << x for x in range(8))
            walsh = [8-2*(f ^ c).bit_count() for c in char3]
            support = {u for u, value in enumerate(walsh) if value}
            assert len(support) == 4 and all(abs(walsh[u]) == 4 for u in support)
            anf = [(f >> x) & 1 for x in range(8)]
            for j in range(3):
                for x in range(8):
                    if x >> j & 1:
                        anf[x] ^= anf[x ^ (1 << j)]
            assert anf[7] == 0
            matrix = [sum(anf[(1 << i) | (1 << j)] << j for j in range(3) if j != i) for i in range(3)]
            assert len(basis(matrix)) == 2
            supports.append(support); branch_checks += 1
        assert not supports[0] & supports[1] and supports[0] | supports[1] == set(range(8))
    plane_rows = []
    for a in range(2):
        for b in range(2):
            v = (xz if a else 0) ^ (yz if b else 0)
            for w in range(2):
                f = sum((((base ^ v) >> (x+8*w)) & 1) << x for x in range(8))
                support = [u for u, c in enumerate(char3) if 8-2*(f ^ c).bit_count()]
                predicted = [u for u in range(8) if (((u >> 2) & 1) ^ (b & (u & 1)) ^ (a & ((u >> 1) & 1)) ^ (a & b)) == w]
                assert support == predicted
                plane_rows.append(dict(a=a, b=b, w=w, frequency_support=support))

    seed = dict(schema='r016-four-variable-affine-seed-v1', old_domain_bits='x0,x1,y0,y1',
                old_base_phase='x0*y0+x1*y1', old_base_phase_mask=f0,
                affine_basis=list(affine), relative_bent_phase_masks=sorted(allowed),
                seed_count=len(allowed), product_family_count=2*len(allowed)**2,
                product_domain_bits='x0,x1,y0,y1,z,w')
    args.out.write_text(json.dumps(seed, indent=2)+'\n')
    out = dict(status='passed', truth_tables_checked=65536, bent_truth_tables=896,
               affine_dimension=5, bent_quotient_classes=28,
               valid_quotient_subspaces_by_dimension=[len(s) for s in levels],
               maximum_base_extension_dimension=maxdim, maximum_base_extension_bases=maxbases,
               inclusion_maximal_base_extensions_by_dimension=maximal_counts,
               maximum_base_extension_capacities=max_capacities,
               base_completion_capacity=8,
               selector_small_test=dict(base_dimension=3, selector_dimension=7,
                                        arbitrary_extra_directions_checked=65536,
                                        common_legal_directions=len(actual), exactly_selector=True),
               coupled_boundary_example=dict(base_phase_mask=base, coupled_basis=list(coupled),
                                             external_direction=xz, its_z_derivative='old x, nonconstant',
                                             obstructing_full_selector_direction=yw,
                                             obstruction_phase_identity='xy+zw+xz+yw=(x+w)(y+z), Walsh magnitude 0 or 8 rather than 4',
                                             seven_dimensional_extension_basis=list(extension),
                                             full_extension_phases_checked=len(span(extension))),
               constant_only_rank2_reduction=dict(extension_basis=list(rank2_extension),
                                                  branch_spectra_checked=branch_checks,
                                                  all_four_point_flat=True, all_quadratic_rank_two=True,
                                                  support_equation_rows=plane_rows,
                                                  joint_support_nonaffine_witness='At w=0, (a,b,beta)=0,(1,0,0),(0,1,0) belong, but (1,1,0) does not.'),
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               seed_sha256=hashlib.sha256(args.out.read_bytes()).hexdigest(),
               python=sys.version, elapsed_seconds=time.monotonic()-began,
               scope='Exhaustive on four-variable truth tables and base quotients; no enumeration of all six-variable truth tables. Large family completeness relies on the separate rigidity proof.')
    args.report.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','bent_truth_tables','valid_quotient_subspaces_by_dimension','maximum_base_extension_dimension','elapsed_seconds']}))


if __name__ == '__main__':
    main()
