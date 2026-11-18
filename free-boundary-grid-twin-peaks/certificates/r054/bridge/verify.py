from collections import defaultdict, deque
from functools import cache
from itertools import permutations
from math import perm


def grid(n):
    return {
        (i, j): {(a, b) for a, b in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
                 if 0 <= a < n and 0 <= b < n}
        for i in range(n) for j in range(n)
    }


def legal_suffix(adj, seeds, suffix):
    seen = set(seeds)
    for v in suffix:
        if v in seen or not (adj[v] & seen):
            return False
        seen.add(v)
    return len(seen) == len(adj)


def ordered_path(adj, seeds, pair):
    paths = {v: (v,) for v in seeds}
    queue = deque(sorted(seeds))
    while queue:
        v = queue.popleft()
        if v in pair:
            other, = set(pair) - {v}
            return paths[v][1:-1], v, other
        for w in sorted(adj[v] - paths.keys()):
            paths[w] = paths[v] + (w,)
            queue.append(w)


def image_suffix(block, far, suffix):
    promoted = set(block[:-2]) | set(far)
    return block + tuple(v for v in suffix if v not in promoted)


def exact_small(n):
    adj = grid(n)
    all_vertices = tuple(sorted(adj))
    bits = {v: 1 << i for i, v in enumerate(all_vertices)}
    neighbors = {v: sum(bits[w] for w in adj[v]) for v in all_vertices}
    full = (1 << len(adj)) - 1

    @cache
    def h(mask):
        if mask == full:
            return 1
        return sum(h(mask | bits[v]) for v in all_vertices
                   if not (mask & bits[v]) and (mask & neighbors[v]))

    prefixes = {frozenset({v}) for v in all_vertices}
    checked = 0
    for k in range(1, 4):
        for K in prefixes:
            mask = sum(bits[v] for v in K)
            shell = set.union(*(adj[v] for v in K)) - K
            U = adj.keys() - K - shell
            bad = sum(h(mask | bits[q] | bits[u]) for q in U for u in adj[q] & shell)
            second_into_U = sum(h(mask | bits[u] | bits[q])
                                for u in shell for q in (adj[u] | set.union(*(adj[v] for v in K))) & U)
            assert bad == second_into_U <= h(mask)
            checked += 1
        prefixes = {K | {v} for K in prefixes
                    for v in set.union(*(adj[w] for w in K)) - K}

    p = (1, 1)
    q = (0, 0)
    assert q not in adj[p] and adj[q] <= adj[p]
    trap_weight = sum(h(bits[p] | bits[q] | bits[u]) for u in adj[q])
    print(f"n={n}: {checked} K-shell identities, singleton trap weight={trap_weight}")


def exhaustive_n3():
    adj = grid(3)
    seed = (0, 0)
    U = adj.keys() - {seed} - adj[seed]
    images = 0
    largest_fiber = 0
    for y in U:
        for z in adj[y] & U:
            fibers = defaultdict(int)
            for q in U:
                for u in adj[q]:
                    inside, a, b = ordered_path(adj, {seed}, (q, u))
                    block = inside + (a, b)
                    if set(block) & {y, z}:
                        continue
                    remaining = adj.keys() - {seed, q, u}
                    for suffix in permutations(sorted(remaining)):
                        if not legal_suffix(adj, {seed, q, u}, suffix):
                            continue
                        out = image_suffix(block, (y, z), suffix)
                        assert legal_suffix(adj, {seed, y, z}, out)
                        fibers[(len(inside) + 1, out)] += 1
                        images += 1
            for (distance, _), count in fibers.items():
                directions = 1 if distance == 1 else 2
                assert count <= directions * perm(len(adj) - 3, distance + 1)
                largest_fiber = max(largest_fiber, count)
    print(f"n=3: {images} actual maps, largest observed joint fiber={largest_fiber}")


def fiber_family(n, distance):
    adj = grid(n)
    p = (0, 0)
    q, u = ((0, 2), (0, 1)) if distance == 1 else ((0, distance), (0, distance + 1))
    y, z = (n - 2, n - 2), (n - 2, n - 1)
    inside, a, b = ordered_path(adj, {p}, (q, u))
    assert len(inside) + 1 == distance
    expected_inside = tuple((0, j) for j in range(1, distance))
    assert inside == expected_inside
    block = inside + (a, b)
    fixed_prefix = (p, q, u) + tuple((1, j) for j in range(n)) + tuple(
        (i, n - 1) for i in range(2, n - 2)
    ) + ((n - 3, n - 2),)
    assert len(fixed_prefix) == 2 * n
    assert q in adj[u] and q not in adj[p]
    fixed_suffix = fixed_prefix[3:]
    moved = inside + (y, z)
    seen = set(fixed_prefix)
    assert not (set(moved) & seen)
    assert all(adj[v] & seen for v in moved)
    rest = adj.keys() - seen - set(moved)
    unused_backbone = []
    while rest:
        v = min(v for v in rest if adj[v] & seen)
        rest.remove(v)
        seen.add(v)
        unused_backbone.append(v)
    all_slots = tuple(range(len(adj) - 2 * n))
    output = image_suffix(block, (y, z), fixed_suffix + moved + tuple(unused_backbone))
    assert legal_suffix(adj, {p, y, z}, output)

    count = 0
    for selected in permutations(all_slots, len(moved)):
        positions = dict(zip(selected, moved))
        background = iter(unused_backbone)
        tail = tuple(positions[i] if i in positions else next(background) for i in all_slots)
        suffix = fixed_suffix + tail
        assert legal_suffix(adj, {p, q, u}, suffix)
        assert image_suffix(block, (y, z), suffix) == output
        count += 1
    assert count == perm(len(adj) - 2 * n, distance + 1)
    print(f"n={n}, D={distance}: {count} distinct verified preimages, same output")


if __name__ == "__main__":
    exact_small(3)
    exact_small(4)
    exhaustive_n3()
    for d in (1, 2, 3):
        fiber_family(6, d)
