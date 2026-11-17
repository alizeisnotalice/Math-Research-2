#!/usr/bin/env python3
"""Complete exact F2 enumeration of affine constant-rank-two 4x4 alternating forms.

The six coordinates are (01,02,03,12,13,23), packed into a number 0..63.
All vector-space and matrix operations are integer/F2 operations; no floating
point arithmetic or random sampling is used.
"""
import datetime as dt
import hashlib
import json
import platform
from pathlib import Path
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def popcount(x):
    return bin(x).count('1')


def xor_span(vectors):
    values = {0}
    for vector in vectors:
        values |= {x ^ vector for x in list(values)}
    return values


def rank_vectors(vectors):
    pivots = {}
    for value in vectors:
        x = value
        while x:
            bit = x.bit_length() - 1
            if bit in pivots:
                x ^= pivots[bit]
            else:
                pivots[bit] = x
                break
    return len(pivots)


def canonical_basis(subspace):
    pivots = {}
    for value in sorted(subspace, reverse=True):
        x = value
        while x:
            bit = x.bit_length() - 1
            if bit in pivots:
                x ^= pivots[bit]
            else:
                pivots[bit] = x
                break
    return tuple(pivots[k] for k in sorted(pivots, reverse=True))


def all_subspaces_six():
    levels = [{1}]  # 64-bit membership set: initial subspace {0}.
    for dimension in range(6):
        next_level = set()
        for membership in levels[-1]:
            elements = [x for x in range(64) if membership >> x & 1]
            for vector in range(1, 64):
                if membership >> vector & 1:
                    continue
                extended = membership
                for x in elements:
                    extended |= 1 << (x ^ vector)
                next_level.add(extended)
        levels.append(next_level)
    expected = [1, 63, 651, 1395, 651, 63, 1]
    if [len(level) for level in levels] != expected:
        raise AssertionError('subspace enumeration count mismatch')
    return levels


def alternating_rows(form):
    c = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    rows = [0, 0, 0, 0]
    for bit, (i, j) in enumerate(c):
        if form >> bit & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rows


def form_rank(form):
    return rank_vectors(alternating_rows(form))


def image_membership(form):
    columns = alternating_rows(form)  # symmetric, so rows are columns
    return sum(1 << x for x in xor_span(columns))


def membership_elements(membership, ambient=16):
    return [x for x in range(ambient) if membership >> x & 1]


RANK = [form_rank(x) for x in range(64)]
IMAGE = [image_membership(x) for x in range(64)]


def affine_cosets(L_membership):
    L = membership_elements(L_membership, 64)
    assigned = 0
    for anchor in range(64):
        if assigned >> anchor & 1:
            continue
        coset = [anchor ^ x for x in L]
        for x in coset:
            assigned |= 1 << x
        yield anchor, coset
    if assigned != (1 << 64) - 1:
        raise AssertionError('coset coverage failure')


def classify_coset(forms):
    common = (1 << 16) - 1
    image_vectors = []
    for form in forms:
        common &= IMAGE[form]
        image_vectors.extend(membership_elements(IMAGE[form]))
    common_dimension = popcount(common).bit_length() - 1
    core_dimension = rank_vectors(image_vectors)
    return common_dimension, core_dimension


def ambient_rank(n, edge_bits, form):
    rows = [0] * n
    for bit, (i, j) in enumerate(edge_bits):
        if form >> bit & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rank_vectors(rows)


def higher_dim_samples():
    # In F2^8: {e0^e7 + sum_{r=1}^6 t_r e0^er}; affine dimension six.
    n = 8
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    edge_index = {edge: k for k, edge in enumerate(edges)}
    def edge(i, j):
        return 1 << edge_index[(min(i, j), max(i, j))]
    star_A = edge(0, 7)
    star_L = [edge(0, r) for r in range(1, 7)]
    star_forms = [star_A ^ x for x in xor_span(star_L)]
    if len(star_forms) != 64 or any(ambient_rank(n, edges, x) != 2 for x in star_forms):
        raise AssertionError('high-dimensional star construction failed')
    # A rank-four difference line embedded in the first four coordinates.
    core_A = edge(0, 1)
    core_difference = edge(0, 1) ^ edge(2, 3)
    core_forms = [core_A, core_A ^ core_difference]
    if any(ambient_rank(n, edges, x) != 2 for x in core_forms):
        raise AssertionError('embedded four-coordinate sample failed')
    if ambient_rank(n, edges, core_difference) != 4:
        raise AssertionError('claimed rank-four difference is wrong')
    return {
        'ambient_dimension': n,
        'star_family': {'affine_dimension': 6, 'forms': len(star_forms),
                        'all_ranks': sorted({ambient_rank(n, edges, x) for x in star_forms}),
                        'fixed_common_factor_vector': 1},
        'embedded_four_coordinate_rank4_difference_line': {
            'affine_dimension': 1, 'forms': len(core_forms),
            'all_ranks': [ambient_rank(n, edges, x) for x in core_forms],
            'difference_rank': ambient_rank(n, edges, core_difference),
            'coordinate_support_bound': 4,
        },
    }


def main():
    root = Path(__file__).resolve().parent
    started = time.monotonic()
    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    levels = all_subspaces_six()
    counts = {str(d): 0 for d in range(7)}
    all_affine = {str(d): 0 for d in range(7)}
    rank4_difference = {str(d): 0 for d in range(7)}
    star = core = both = 0
    nonrank4_family_count = nonrank4_star = nonrank4_core = nonrank4_both = 0
    nonrank4_neither = 0
    nonrank4_neither_witness = None
    max_rank4_difference_dimension = -1
    rank4_witness = None
    strict_witness = None
    for d, level in enumerate(levels):
        for L_membership in level:
            L = membership_elements(L_membership, 64)
            has_rank4_difference = any(RANK[x] == 4 for x in L)
            for anchor, coset in affine_cosets(L_membership):
                all_affine[str(d)] += 1
                if not all(RANK[x] == 2 for x in coset):
                    continue
                counts[str(d)] += 1
                if has_rank4_difference:
                    rank4_difference[str(d)] += 1
                    max_rank4_difference_dimension = max(max_rank4_difference_dimension, d)
                    if rank4_witness is None:
                        rank4_witness = {'dimension': d, 'anchor': anchor,
                                         'L_basis': list(canonical_basis(L)), 'forms': coset}
                common_dimension, core_dimension = classify_coset(coset)
                is_star = common_dimension >= 1
                is_core = core_dimension <= 3
                star += is_star
                core += is_core
                both += is_star and is_core
                if not has_rank4_difference:
                    nonrank4_family_count += 1
                    nonrank4_star += is_star
                    nonrank4_core += is_core
                    nonrank4_both += is_star and is_core
                    if not is_star and not is_core:
                        nonrank4_neither += 1
                        if nonrank4_neither_witness is None:
                            nonrank4_neither_witness = {
                                'dimension': d, 'anchor': anchor,
                                'L_basis': list(canonical_basis(L)), 'forms': coset,
                                'common_support_dimension': common_dimension,
                                'sum_support_dimension': core_dimension,
                            }
                if strict_witness is None and not has_rank4_difference and not is_star and is_core:
                    strict_witness = {'dimension': d, 'anchor': anchor,
                                      'L_basis': list(canonical_basis(L)), 'forms': coset,
                                      'common_support_dimension': common_dimension,
                                      'sum_support_dimension': core_dimension}
    candidate_classification_ok = (nonrank4_neither == 0 and
                                   max_rank4_difference_dimension <= 3)
    expected_affine = {'0': 64, '1': 2016, '2': 10416, '3': 11160,
                       '4': 2604, '5': 126, '6': 1}
    if all_affine != expected_affine:
        raise AssertionError('affine-subspace coverage count mismatch')
    result = {
        'schema': 'a28-r017-affine-alternating-4x4-v1',
        'status': 'passed' if candidate_classification_ok else 'classification_counterexample_found',
        'ambient_dimension': 4,
        'coordinate_order': ['01', '02', '03', '12', '13', '23'],
        'enumerated_all_linear_subspaces_by_dimension': [len(x) for x in levels],
        'enumerated_all_affine_subspaces_by_dimension': all_affine,
        'constant_rank_2_affine_families_by_dimension': counts,
        'rank4_difference_constant_rank_2_families_by_dimension': rank4_difference,
        'max_dimension_with_rank4_difference': max_rank4_difference_dimension,
        'classification_checks': {
            'common_star_factor_count': star,
            'fixed_core_dimension_at_most_3_count': core,
            'both_count': both,
            'nonrank4_difference_family_count': nonrank4_family_count,
            'nonrank4_difference_common_star_count': nonrank4_star,
            'nonrank4_difference_fixed_core_dimension_at_most_3_count': nonrank4_core,
            'nonrank4_difference_both_count': nonrank4_both,
            'nonrank4_difference_neither_count': nonrank4_neither,
            'rank4_difference_dimension_at_most_3': max_rank4_difference_dimension <= 3,
        },
        'rank4_difference_witness': rank4_witness,
        'nonstar_fixed_core_witness': strict_witness,
        'nonrank4_difference_neither_witness': nonrank4_neither_witness,
        'higher_dim_exact_samples': higher_dim_samples(),
        'rank_histogram_all_4x4_alternating_forms': {str(r): RANK.count(r) for r in (0, 2, 4)},
        'source_sha256': sha(__file__), 'python': sys.version, 'platform': platform.platform(),
        'started_utc': started_utc, 'elapsed_seconds': time.monotonic() - started,
        'scope': 'Complete F2 4x4 enumeration plus two exact higher-dimensional constructions; finite evidence only, not a proof of the candidate classification in arbitrary dimension.',
    }
    out = root / 'enumeration-result.json'
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'constant_rank2': counts,
                      'max_rank4_difference_dimension': max_rank4_difference_dimension,
                      'nonrank4_neither': nonrank4_neither,
                      'elapsed_seconds': result['elapsed_seconds']}, sort_keys=True))
    return 0 if candidate_classification_ok else 10


if __name__ == '__main__':
    sys.exit(main())
