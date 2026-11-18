"""Exact checks for the rectangular splitting, alternating decomposition and Walsh cuts."""
from collections import Counter
from fractions import Fraction
from functools import cache
import itertools
import json
import random
import sys

SEED = 8607


def bits(vector):
    while vector:
        least = vector & -vector
        yield least.bit_length() - 1
        vector ^= least


def independent(vectors):
    chosen = []
    pivots = {}
    for vector in vectors:
        residual = vector
        while residual:
            key = residual.bit_length()
            if key not in pivots:
                pivots[key] = residual
                chosen.append(vector)
                break
            residual ^= pivots[key]
    return chosen


def rank(vectors):
    return len(independent(vectors))


def combinations(basis):
    span = [0]
    for x in basis:
        span.extend(y ^ x for y in span[:])
    return span


def apply(columns, x):
    result = 0
    for i in bits(x):
        result ^= columns[i]
    return result


def transpose(columns, row_count):
    return tuple(sum(((column >> row) & 1) << index for index, column in enumerate(columns))
                 for row in range(row_count))


def full_rectangle(family, a, b):
    return (rank(column for matrix in family for column in matrix) == a and
            rank(row for matrix in family for row in transpose(matrix, a)) == b)


def walsh(values):
    coefficients = list(values)
    stride = 1
    while stride < len(coefficients):
        for start in range(0, len(coefficients), 2 * stride):
            for i in range(start, start + stride):
                x, y = coefficients[i], coefficients[i + stride]
                coefficients[i], coefficients[i + stride] = x+y, x-y
        stride *= 2
    return coefficients


@cache
def subspaces(dimension):
    all_points = (1 << (1 << dimension)) - 1
    current = {1}
    result = []
    for d in range(dimension + 1):
        result.extend((space, d) for space in sorted(current))
        expanded = set()
        for space in current:
            outside = all_points ^ space
            elements = tuple(bits(space))
            while outside:
                x = (outside & -outside).bit_length() - 1
                coset = sum(1 << (y ^ x) for y in elements)
                expanded.add(space | coset)
                outside ^= coset
        current = expanded
    return result


def normalize(labels, weights, dimension):
    total = sum(weights)
    L = Fraction(0)
    for subspace, d in subspaces(dimension):
        if d == dimension:
            continue
        inside = sum(w for x, w in zip(labels, weights) if (subspace >> x) & 1)
        L = max(L, Fraction(dimension - d + 1, total - inside))
    return tuple(L*w for w in weights)


def weight_case(old_dimension, rng):
    affine_marker = 1 << (old_dimension - 1)
    required = {affine_marker, *(affine_marker | 1 << i for i in range(old_dimension - 1))}
    labels = tuple(x for x in range(affine_marker, 2*affine_marker)
                   if x in required or rng.getrandbits(1))
    weights = tuple(rng.randrange(1, 10) for _ in labels)
    return labels, normalize(labels, weights, old_dimension)


def verify_two_functions(family, weights, a, b, half_rank=False):
    ranks = [rank(matrix) for matrix in family]
    if half_rank:
        assert all(r % 2 == 0 for r in ranks)
        ranks = [r//2 for r in ranks]
    assert sum(w*((1<<r)-1) for w,r in zip(weights, ranks)) >= a+b
    assert sum(w*((1<<r)-Fraction(1,1<<r)) for w,r in zip(weights,ranks)) >= a+b+1


def audit_rectangle(family, weights, a, b, counts):
    assert full_rectangle(family, a, b)
    verify_two_functions(family, weights, a, b)
    if a < b:
        family = tuple(transpose(matrix, a) for matrix in family)
        a, b = b, a
        counts["rectangular_transpositions"] += 1
    for x in range(1, 1 << b):
        image_basis = independent(apply(matrix, x) for matrix in family)
        if len(image_basis) < a:
            break
    else:
        counts["rectangular_surjective_terminals"] += 1
        for x in range(1, 1 << b):
            assert sum(w for w,M in zip(weights,family) if apply(M,x)) >= a+1
        return
    counts["rectangular_nontrivial_splits"] += 1
    c = len(image_basis)
    C = set(combinations(image_basis))
    kernel_basis = independent(y for y in range(1 << b)
                               if all(apply(M,y) in C for M in family))
    d = len(kernel_basis)
    assert 0 < c < a and 0 < d < b
    image_extended = independent([*image_basis, *(1 << i for i in range(a))])
    lookup = {point: coordinate for coordinate,point in enumerate(combinations(image_extended))}
    kernel_extended = independent([*kernel_basis, *(1 << i for i in range(b))])
    small = tuple(tuple(lookup[apply(M,y)] for y in kernel_basis) for M in family)
    quotient = tuple(tuple(lookup[apply(M,y)] >> c for y in kernel_extended[d:]) for M in family)
    assert all(rank(M) >= rank(Ms)+rank(Mq) for M, Ms, Mq in zip(family, small, quotient))
    audit_rectangle(small, weights, c, d, counts)
    audit_rectangle(quotient, weights, a-c, b-d, counts)


def change_form_basis(family, chosen):
    n = len(chosen)
    lookup = {point: index for index, point in enumerate(combinations(chosen))}
    inverse = tuple(lookup[1 << i] for i in range(n))
    dual_basis = transpose(inverse, n)
    transformed = tuple(tuple(lookup[apply(M, y)] for y in dual_basis) for M in family)
    assert all(M == transpose(M, n) for M in transformed)
    return transformed


def audit_alternating(family, weights, n, counts):
    assert rank(column for M in family for column in M) == n
    verify_two_functions(family, weights, n, 0, half_rank=True)
    counts["alternating_splits"] += 1
    S = independent(M[0] for M in family)
    assert S and all((x & 1) == 0 for x in S)
    extended = independent([*S, *(1 << i for i in range(1, n))])
    H = extended[len(S):]
    m = len(H)
    base = change_form_basis(family, [*H, 1, *S])
    mask = (1<<m)-1
    core_images = independent(column & mask for M in base for column in M[:m])
    a = len(core_images)
    residual = independent([*core_images, *(1 << i for i in range(m))])[a:]
    newer = change_form_basis(base, [*core_images, *residual, *(1 << i for i in range(m,n))])
    core = tuple(tuple(col & ((1<<a)-1) for col in M[:a]) for M in newer)
    mid = (1 << (m+1-a))-1
    cross = tuple(tuple((col>>a)&mid for col in M[m+1:]) for M in newer)
    assert all(rank(M) >= rank(A)+2*rank(B) for M,A,B in zip(family, core, cross))
    if a:
        audit_alternating(core, weights, a, counts)
    else:
        counts["zero_alternating_cores"] += 1
    audit_rectangle(cross, weights, m+1-a, len(S), counts)


def linear_family(generators, labels, columns):
    return tuple(tuple(apply(tuple(M[i] for M in generators), x) for i in range(columns))
                 for x in labels)


def alternating_matrix(n, rng):
    columns = [0]*n
    for i,j in itertools.combinations(range(n),2):
        if rng.getrandbits(1):
            columns[i] ^= 1<<j
            columns[j] ^= 1<<i
    return tuple(columns)


def abstract_tests(rng):
    counts = Counter()
    for v in range(2, 6):
        for _ in range(12):
            labels, weights = weight_case(v,rng)
            counts["fractional_weight_systems"] += any(w < 1 for w in weights)
            for a,b in itertools.product(range(1,6),repeat=2):
                generators = [tuple(rng.getrandbits(a) for _ in range(b)) for _ in range(v)]
                family = linear_family(generators,labels,b)
                if full_rectangle(family,a,b):
                    audit_rectangle(family,weights,a,b,counts)
                    counts["root_rectangular_families"] += 1
            for n in range(2,13):
                generators = [alternating_matrix(n,rng) for _ in range(v)]
                family = linear_family(generators,labels,n)
                if rank(col for M in family for col in M) == n:
                    audit_alternating(family,weights,n,counts)
                    counts["root_alternating_families"] += 1
    return dict(counts)


def quadratic(M, address):
    return sum(((M[i] & address) >> (i+1)).bit_count() for i in bits(address)) & 1


def sampled_cut_value(coefficients, rng, trials):
    labels = [i | len(coefficients) for i,c in enumerate(coefficients) if c]
    effective_basis = independent(labels)
    dimension = len(effective_basis)
    lookup = {point:i for i,point in enumerate(combinations(effective_basis))}
    terms = [(lookup[label], abs(coefficients[label ^ len(coefficients)])) for label in labels]
    if trials is None:
        subsets = subspaces(dimension)
    else:
        subsets = []
        for _ in range(trials):
            d = rng.randrange(dimension)
            chosen = []
            while len(chosen) < d:
                vector = rng.getrandbits(dimension)
                if rank([*chosen,vector]) > len(chosen):
                    chosen.append(vector)
            subsets.append((sum(1<<v for v in combinations(chosen)), d))
        subsets.append((1,0))
    best = Fraction(0)
    count = 0
    for space,d in subsets:
        if d < dimension:
            outside = sum(weight for label,weight in terms if not ((space>>label)&1))
            best = max(best,Fraction((dimension-d+1)*len(coefficients),outside))
            count += 1
    return best, count


def walsh_tests(rng):
    checked_cuts = 0
    parents = 0
    for v in (2,3,4):
        for h in (2,3,5,7):
            for _ in range(8):
                data = tuple((-1)**rng.getrandbits(1) for _ in range(1<<v))
                base, _ = sampled_cut_value(walsh(data), rng, None)
                generators = tuple(alternating_matrix(h,rng) for _ in range(v+1))
                signs = []
                for point in range(1 << (v+h)):
                    u = point & ((1<<h)-1)
                    x = point >> h
                    x ^= sum(quadratic(form,u)<<i for i,form in enumerate(generators[1:]))
                    signs.append((-1)**quadratic(generators[0],u)*data[x])
                cut_value, count = sampled_cut_value(walsh(signs),rng,96)
                assert cut_value <= base
                checked_cuts += count
                parents += 1
    return {"true_parents":parents, "sampled_parent_cuts":checked_cuts}


if __name__ == "__main__":
    rng = random.Random(SEED)
    print(json.dumps({"python":sys.version.split()[0], "seed":SEED,
                      "abstract":abstract_tests(rng), "real":walsh_tests(rng)}, indent=2))
