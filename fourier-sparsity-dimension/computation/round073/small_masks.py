"""Produce affine-normalized representative negatives of size at most eight."""
from itertools import combinations

masks = {0}
for rank in range(6):
    basis = {0, *(1 << i for i in range(rank))}
    available = sorted(set(range(1 << rank)) - basis)
    for extra in range(8 - len(basis) + 1):
        for points in combinations(available, extra):
            masks.add(sum(1 << x for x in basis.union(points)))
for mask in sorted(masks):
    print(mask)
