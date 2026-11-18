#!/usr/bin/env python3
"""Exact greedy-prefix, density, signed-convolution, and common-query audits."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
import math
import platform
import random
import time


def dot(x, y):
    return bin(x & y).count('1') % 2


def character(x, y):
    return 1 - 2 * dot(x, y)


def basis(vectors):
    rows = []
    for vector in vectors:
        vector = reduce_by(vector, rows)
        if vector:
            rows.append(vector)
    return tuple(rows)


def reduce_by(vector, rows):
    for row in rows:
        if vector & (1 << (row.bit_length() - 1)):
            vector ^= row
    return vector


def extend_space(space, vector):
    return space.union(x ^ vector for x in space)


def rank_space(space):
    return len(space).bit_length() - 1


@lru_cache(maxsize=None)
def all_subspaces(n):
    levels = [(frozenset({0}),)]
    for dimension in range(n):
        next_level = set()
        for space in levels[-1]:
            for representative in quotient_representatives(space, n):
                if representative:
                    next_level.add(extend_space(space, representative))
        levels.append(tuple(sorted(next_level, key=space_key)))
    return levels


def space_key(space):
    return tuple(sorted(space))


def quotient_representatives(space, n):
    rows = basis(space_key(space))
    pivots = sum(1 << (row.bit_length() - 1) for row in rows)
    return [x for x in range(1 << n) if not (x & pivots)]


def quotient_bins(coefficients, core):
    rows = basis(space_key(core))
    bins = defaultdict(list)
    for x, coefficient in coefficients.items():
        bins[reduce_by(x, rows)].append((x, coefficient))
    return dict(bins)


def evaluate_bins(bins, physical):
    return {q: value for q, entries in bins.items()
            if (value := sum(a * character(x, physical) for x, a in entries))}


def audit_convolution(coefficients, denominator):
    energy = sum(a * a for a in coefficients.values())
    convolutions = Counter()
    for (u, a), (v, b) in combinations(sorted(coefficients.items()), 2):
        convolutions[u ^ v] += a * b
    if energy != denominator ** 2 or any(convolutions.values()):
        raise AssertionError(('not a Boolean spectrum', energy, denominator, convolutions))
    return dict(squared_numerator_sum=energy, denominator=denominator,
                nonzero_xor_equations=len(convolutions), integer_coefficients=sorted(coefficients.items()))


def transform(values):
    output = list(values)
    stride = 1
    while stride < len(output):
        for start in range(0, len(output), 2 * stride):
            for x in range(start, start + stride):
                a, b = output[x], output[x + stride]
                output[x], output[x + stride] = a + b, a - b
        stride *= 2
    return output


def truth_spectrum(truth):
    return {x: a for x, a in enumerate(transform(truth)) if a}, len(truth)


def indicator_bands(cells, core_dim):
    denominator = len(cells)
    coefficients, weights = {}, []
    if sorted(x for cell in cells for x in cell) != list(range(1 << core_dim)):
        raise AssertionError('address cells do not partition their cube')
    for i, cell in enumerate(cells):
        band = {}
        for q in range(1 << core_dim):
            numerator = sum(character(q, x) for x in cell)
            if numerator:
                numerator, residual = divmod(numerator * denominator, 1 << core_dim)
                if residual:
                    raise AssertionError('invalid address indicator numerator')
                band[q] = numerator
                coefficients[(1 << (core_dim + i)) ^ q] = numerator
        weights.append(frozenset(band))
    return coefficients, denominator, weights


def address_cells(bits):
    return [[x] for x in range(1 << bits)]


def r49_cells(height):
    pairs = ((0, 1), (2, 4), (3, 7), (5, 6))
    internal_per_tree = (1 << (height - 1)) - 1
    core_dim = 3 + 4 * internal_per_tree
    cells = [[] for _ in range(1 << (height + 1))]
    leaves_per_tree = 1 << (height - 1)
    for x in range(1 << core_dim):
        tree = next(i for i, pair in enumerate(pairs) if (x & 7) in pair)
        node = 0
        while node < internal_per_tree:
            bit = (x >> (3 + tree * internal_per_tree + node)) & 1
            node = 2 * node + 1 + bit
        cells[tree * leaves_per_tree + node - internal_per_tree].append(x)
    return core_dim, cells


def weak_folding_cells(h):
    # Tarannikov Lemmas 3–4: use the affine slice of spread graphs at an x-bit of 1.
    q = 1 << (h + 1)
    irreducible = {1: 0b111, 2: 0b1011}[h]

    def multiply(x, y):
        answer = 0
        while y:
            if y & 1:
                answer ^= x
            y >>= 1
            x <<= 1
            if x & q:
                x ^= irreducible
        return answer

    cells = [[((multiply(2 * x + 1, alpha) << h) ^ x) for x in range(1 << h)]
             for alpha in range(q)]
    return 2 * h + 1, cells


def fraction_json(value):
    return [value.numerator, value.denominator]


def density_general(n, support):
    support = set(support) - {0}
    best, winner = Fraction(0), frozenset({0})
    for level in all_subspaces(n)[1:]:
        for subspace in level:
            ratio = Fraction(len(support & subspace), rank_space(subspace))
            if ratio > best:
                best, winner = ratio, subspace
    return density_result(best, winner, support, 'all ambient subspaces')


def density_bands(core_dim, bands, support):
    best, winner = Fraction(0), frozenset({0})
    for level in all_subspaces(core_dim):
        for core in level:
            weights = sorted(((len(band & core), index) for index, band in enumerate(bands)), reverse=True)
            cumulative = 0
            for j, (weight, _) in enumerate(weights, start=1):
                cumulative += weight
                ratio = Fraction(cumulative, rank_space(core) + j)
                if ratio > best:
                    best, winner = ratio, core
                    for _, index in weights[:j]:
                        winner = extend_space(winner, 1 << (core_dim + index))
    return density_result(best, winner, set(support), 'all core subspaces and optimal active-band subsets')


def density_result(ratio, winner, support, route):
    points = sorted(winner & support)
    return dict(density=fraction_json(ratio), c=math.ceil(ratio), checked_method=route,
                witness_basis=basis(space_key(winner)), witness_support_points=points)


def active_leaves(n, coefficients, denominator, core):
    bins = quotient_bins(coefficients, core)
    annihilator = quotient_representatives(core, n)
    # Representatives for physical points modulo H^perp use the pivot positions of H.
    pivot_mask = reduce_mask(n, annihilator)
    shifts = [x for x in range(pivot_mask + 1) if not (x & ~pivot_mask)]
    leaves = [(shift, evaluate_bins(bins, shift)) for shift in shifts]
    if set().union(*(set(leaf) for _, leaf in leaves)) != set(bins):
        raise AssertionError('incorrect quotient leaf ensemble')
    for _, leaf in leaves:
        audit_convolution(leaf, denominator)
    return bins, leaves


def reduce_mask(n, quotient):
    return ((1 << n) - 1) ^ basis_union(basis(quotient))


def basis_union(vectors):
    mask = 0
    for x in vectors:
        mask |= x
    return mask


def audit_basis(n, coefficients, denominator, core, exhaustive_basis):
    bins, leaves = active_leaves(n, coefficients, denominator, core)
    weights = Counter(q for _, leaf in leaves for q in leaf)
    labels = sorted(set(bins) - {0}, key=lambda q: (-weights[q], q))
    chosen = basis_label_subset(labels)
    chosen_set = set(chosen)
    r = len(chosen)
    if exhaustive_basis:
        all_weights = [sum(weights[x] for x in candidate) for candidate in combinations(labels, r)
                       if len(basis(candidate)) == r]
        if len(labels) and sum(weights[x] for x in chosen) != max(all_weights):
            raise AssertionError('maximum-weight basis optimizer disagrees with enumeration')
    activity = Counter()
    for _, leaf in leaves:
        k, t = len(set(leaf) & chosen_set), len(set(leaf) - chosen_set)
        if k * (k - 1) > 6 * t + t * (t - 1):
            raise AssertionError(('pointwise pair budget', k, t))
        activity[k, t] += 1
    for q in chosen:
        if len(bins[q]) * weights[q] < len(leaves):
            raise AssertionError('uncertainty check failed')
    expected_t = Fraction(sum(t * count for (k, t), count in activity.items()), len(leaves))
    if r * r > len(coefficients) * (2 + expected_t) or (not expected_t and r * r > len(coefficients)):
        raise AssertionError('mean pair budget')
    return dict(L=len(bins), quotient_rank=r, maximum_weight_basis=chosen,
                mean_t=fraction_json(expected_t), activity_histogram=[[k,t,count] for (k,t),count in sorted(activity.items())],
                bin_size_histogram=sorted(Counter(map(len, bins.values())).items()))


def basis_label_subset(labels):
    rows, selected = [], []
    for label in labels:
        reduced = reduce_by(label, rows)
        if reduced:
            selected.append(label)
            rows.append(reduced)
    return selected


def all_tied_frontier(n, coefficients, denominator, depth, exhaustive_basis):
    paths = {frozenset({0}): 1}
    levels = []
    for step in range(depth + 1):
        next_paths = Counter()
        rows = []
        for core, multiplicity in sorted(paths.items(), key=lambda item: space_key(item[0])):
            result = audit_basis(n, coefficients, denominator, core, exhaustive_basis)
            support = tuple(quotient_bins(coefficients, core))
            folds = Counter(a ^ b for a, b in combinations(support, 2))
            max_drop = max(folds.values(), default=0)
            ext = sorted(v for v, count in folds.items() if count == max_drop)
            if not max_drop:
                ext = [v for v in quotient_representatives(core, n) if v]
            if step < depth:
                for direction in ext:
                    new_core = extend_space(core, direction)
                    if len(quotient_bins(coefficients, new_core)) != len(support) - max_drop:
                        raise AssertionError('pair-count quotient drop mismatch')
                    next_paths[new_core] += multiplicity
            result.update(core_basis=basis(space_key(core)), paths=multiplicity, next_tied_directions=ext,
                          maximum_drop=max_drop)
            rows.append(result)
        costs = [Fraction(*row['mean_t']) for row in rows]
        levels.append(dict(depth=step, path_count=sum(paths.values()), cores=len(rows),
                           min_mean_t=fraction_json(min(costs)), max_mean_t=fraction_json(max(costs)), rows=rows))
        paths = next_paths
    return levels


def minimum_certificate(n, leaf):
    pivots = basis(sorted(leaf))
    compact_n = len(pivots)
    compact_mapping = {0: 0}
    for i, row in enumerate(pivots):
        compact_mapping.update({v ^ row: w ^ (1 << i) for v, w in list(compact_mapping.items())})
    compact_leaf = {compact_mapping[v]: a for v,a in leaf.items()}
    for k, level in enumerate(all_subspaces(compact_n)):
        options = []
        for certificate in level:
            bins = quotient_bins(compact_leaf, certificate)
            mask = reduce_mask(compact_n, quotient_representatives(certificate, compact_n))
            for shift in range(mask + 1):
                if not shift & ~mask and set(evaluate_bins(bins, shift)) == {0}:
                    inv = {w:v for v,w in compact_mapping.items()}
                    options.append(tuple(inv[h] for h in space_key(certificate)))
                    break
        if options:
            return k, options
    raise AssertionError(('no certificate', n, leaf))


def scale_trace(n, coefficients, denominator):
    original_s = len(coefficients)
    core = frozenset({0})
    trace = []
    while True:
        bins, leaves = active_leaves(n, coefficients, denominator, core)
        L = len(bins)
        m = max(len(leaf) for _, leaf in leaves)
        if m * original_s < L ** 2:
            raise AssertionError('max-leaf uncertainty bound')
        if L ** 2 <= original_s:
            added = basis_label_subset(sorted(set(bins) - {0}))
            for direction in added:
                core = extend_space(core, direction)
            trace.append(dict(L=L, m=m, kappa=None, actual_added_rank=len(added), delta_L=L-1,
                              adding=added, terminal=True))
            break
        # Tie rule: first minimize the certificate rank among all maximally sparse leaves;
        # then maximize global quotient drop; finally use fixed lexicographic order.
        candidates = []
        for shift, leaf in leaves:
            if len(leaf) != m:
                continue
            kappa, options = minimum_certificate(n, leaf)
            for certificate in options:
                new_core = frozenset(v ^ w for v in core for w in certificate)
                after = quotient_bins(coefficients, new_core)
                candidates.append((kappa, len(after), shift, tuple(sorted(certificate)), new_core))
        kappa, next_L, shift, certificate, next_core = min(candidates, key=lambda row: row[:4])
        actual_rank = rank_space(next_core) - rank_space(core)
        if 2 * (L - next_L) < m or actual_rank != kappa:
            raise AssertionError(('certificate scale budget', L, next_L, m, kappa, actual_rank))
        band = 0
        while (1 << (2 * (band + 1))) * original_s <= L ** 2:
            band += 1
        trace.append(dict(L=L, m=m, kappa=kappa, actual_added_rank=actual_rank, delta_L=L-next_L,
                          scale_band=band, selected_shift=shift, adding=basis(certificate),
                          inherited_bin_histogram=sorted(Counter(map(len,bins.values())).items()),
                          terminal=False))
        core = next_core
    if rank_space(core) != len(basis(sorted(coefficients))):
        raise AssertionError('common queries did not generate the spectral space')
    return trace


def structured_scale_trace(core_dim, cells, coefficients, denominator):
    """Exact certificate ranks for selectors on disjoint address cells.

    Intersect a candidate monochromatic affine flat with a fixed address fibre.
    Its homogeneous data-fibre direction must kill each hit, unfixed data label.
    Conversely, separately fixing those labels pays that same data codimension.
    """
    n = core_dim + len(cells)
    original_s = len(coefficients)
    address_limit = 1 << core_dim
    address_of = {x: i for i, cell in enumerate(cells) for x in cell}
    queried = frozenset({0})
    address_core = frozenset({0})
    known_data = set()
    trace = []
    geometry_cache = {}
    while True:
        bins, leaves = active_leaves(n, coefficients, denominator, queried)
        L = len(bins)
        m = max(len(leaf) for _, leaf in leaves)
        if m * original_s < L ** 2:
            raise AssertionError('max-leaf uncertainty bound')
        if L ** 2 <= original_s:
            added = basis_label_subset(sorted(set(bins) - {0}))
            for direction in added:
                queried = extend_space(queried, direction)
            trace.append(dict(L=L,m=m,kappa=None,actual_added_rank=len(added),delta_L=L-1,
                              adding=added,terminal=True))
            break
        current_address_rank = rank_space(address_core)
        current_address_basis = basis(space_key(address_core))
        current_address_dim = core_dim - current_address_rank
        compatible_subspaces = [v for level in all_subspaces(core_dim)[0:current_address_dim + 1]
                                for v in level
                                if all(not dot(x, row) for x in basis(space_key(v))
                                       for row in current_address_basis)]
        candidates = []
        for shift, leaf in leaves:
            if len(leaf) != m:
                continue
            address_shift = shift & (address_limit - 1)
            key = (address_core, address_shift)
            if key not in geometry_cache:
                domain = frozenset(x for x in range(address_limit)
                                   if all(not dot(x ^ address_shift, row) for row in current_address_basis))
                affine_flats = []
                for subspace in compatible_subspaces:
                    left = set(domain)
                    while left:
                        x = min(left)
                        coset = {x ^ v for v in subspace}
                        left.difference_update(coset)
                        hit = frozenset(address_of[v] for v in coset)
                        affine_flats.append((current_address_dim - rank_space(subspace), hit, x, subspace))
                geometry_cache[key] = affine_flats
            best = n + 1
            best_flats = []
            for address_cost, hit, x, subspace in geometry_cache[key]:
                fixed_signs = {((shift >> (core_dim + i)) & 1) for i in hit & known_data}
                if len(fixed_signs) > 1:
                    continue
                new_data = sorted(hit - known_data)
                cost = address_cost + len(new_data)
                if cost > best:
                    continue
                if cost < best:
                    best, best_flats = cost, []
                sign = next(iter(fixed_signs), 0)
                best_flats.append((x, subspace, new_data, sign))
            if best == n + 1:
                raise AssertionError('no monochromatic address flat')
            for x, subspace, new_data, sign in best_flats:
                compact_basis = basis(space_key(subspace))
                next_address_core = frozenset(q for q in range(address_limit)
                                              if all(not dot(q,v) for v in compact_basis))
                combined_data = known_data | set(new_data)
                query_rows = basis(space_key(next_address_core)) + tuple(1 << (core_dim+i) for i in sorted(combined_data))
                next_queried = frozenset({0})
                for q in query_rows:
                    next_queried = extend_space(next_queried, q)
                new_bins = quotient_bins(coefficients,next_queried)
                witness = (shift ^ address_shift) | x | sum(sign << (core_dim+i) for i in new_data)
                if evaluate_bins(new_bins, witness) != {0: (1-2*sign)*denominator}:
                    raise AssertionError('invalid exact monochromatic address witness')
                candidates.append((best,len(new_bins),shift,x,space_key(subspace),tuple(new_data),
                                   next_queried,next_address_core,witness))
        kappa,next_L,shift,_,_,new_data,next_queried,address_core,witness = min(candidates,key=lambda row: row[:6])
        actual_rank = rank_space(next_queried) - rank_space(queried)
        if 2 * (L - next_L) < m or actual_rank != kappa:
            raise AssertionError(('structured certificate scale budget', L, next_L, m, kappa, actual_rank))
        new_rows = basis_label_subset([*basis(space_key(queried)),*basis(space_key(next_queried))])
        added = new_rows[rank_space(queried):]
        band = 0
        while (1 << (2 * (band + 1))) * original_s <= L ** 2:
            band += 1
        trace.append(dict(L=L,m=m,kappa=kappa,actual_added_rank=actual_rank,delta_L=L-next_L,
                          scale_band=band,selected_shift=shift,adding=added,terminal=False,
                          flat_witness=witness, inherited_bin_histogram=sorted(Counter(map(len,bins.values())).items())))
        known_data.update(new_data)
        queried = next_queried
    if rank_space(queried) != len(basis(sorted(coefficients))):
        raise AssertionError('structured common queries did not generate the spectral space')
    return trace


def mock_rejection(m):
    n = m + 1
    support = sorted(x ^ extra for x in (1 << i for i in range(1, n)) for extra in (0, 1))
    by_difference = defaultdict(list)
    for a, b in combinations(support, 2):
        by_difference[a ^ b].append((a, b))
    if min(map(len, by_difference.values())) != 2:
        raise AssertionError('unexpected mock pair multiplicity')
    witness_differences = [2 ^ 4, 2 ^ 8, 4 ^ 8]
    parity_rows = []
    for q in witness_differences:
        pairs = by_difference[q]
        row = 0
        for a,b in pairs:
            row ^= (1 << support.index(a)) ^ (1 << support.index(b))
        if len(pairs) != 2:
            raise AssertionError('sign equation is not binary')
        parity_rows.append(row)
    if parity_rows[0] ^ parity_rows[1] ^ parity_rows[2]:
        raise AssertionError('missing signed triangle contradiction')
    # Explicit two-independent-set partition; an intersecting line contains two points.
    parts = [tuple(x for x in support if x % 2 == parity) for parity in (0,1)]
    if any(len(basis(part)) != len(part) for part in parts) or len(basis(support)) >= len(support):
        raise AssertionError('invalid c certificate')
    return dict(m=m, d=len(basis(support)), s=len(support), c=2, independent_parts=parts,
                difference_multiplicity_histogram=sorted(Counter(map(len,by_difference.values())).items()),
                inconsistent_convolution_rows=[dict(difference=q, products=by_difference[q],
                                                    sign_variables_bitmask=row, rhs=1)
                                               for q,row in zip(witness_differences,parity_rows)])


def process(name, n, spectrum, band_core=None, depth=None):
    coefficients, denominator = spectrum[:2]
    print('checking ' + name, file=__import__('sys').stderr, flush=True)
    d = len(basis(sorted(coefficients)))
    density = density_general(n, coefficients) if band_core is None else density_bands(band_core,spectrum[2],coefficients)
    if n <= 7 and band_core is not None:
        if density['density'] != density_general(n, coefficients)['density']:
            raise AssertionError('band density solver disagrees with exhaustive ambient solver')
    default_depth = d // 2
    return dict(name=name, n=n, d=d, s=len(coefficients), density=density,
                convolution=audit_convolution(coefficients, denominator),
                all_tied_prefixes=all_tied_frontier(n, coefficients, denominator,
                                                    default_depth if depth is None else depth, n <= 4),
                scale_trace=scale_trace(n, coefficients, denominator) if n <= 7 else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('small', 'larger', 'mock'))
    args=parser.parse_args()
    wall, cpu = time.perf_counter(), time.process_time()
    if args.mode == 'mock':
        records = [mock_rejection(m) for m in (3, 4, 5, 8, 16, 32)]
    elif args.mode == 'small':
        rng = random.Random(780828)
        records = []
        for bits in (1,2):
            cells = address_cells(bits)
            records.append(process('address-' + str(bits), bits + len(cells), indicator_bands(cells,bits),bits))
        for name, construct in (('R49',r49_cells),('Tarannikov',weak_folding_cells)):
            core, cells = construct(1)
            records.append(process(name+'-1',core+len(cells),indicator_bands(cells,core),core))
        for index in range(8):
            truth = [1-2*rng.randrange(2) for _ in range(16)]
            records.append(process('random-n4-'+str(index),4,truth_spectrum(truth)))
    else:
        records = []
        for name, construct in (('R49',r49_cells),('Tarannikov',weak_folding_cells)):
            core, cells = construct(2)
            spectrum = indicator_bands(cells,core)
            result = process(name+'-2',core+len(cells),spectrum,core)
            result['scale_trace'] = structured_scale_trace(core,cells,*spectrum[:2])
            records.append(result)
    wall_seconds, cpu_seconds = time.perf_counter()-wall, time.process_time()-cpu
    print(json.dumps(dict(status='pass',mode=args.mode,host=platform.node(),python=platform.python_version(),
                          wall_seconds=wall_seconds,process_cpu_seconds=cpu_seconds,
                          mean_cpu_percent_of_one_core=100*cpu_seconds/wall_seconds,records=records),
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
