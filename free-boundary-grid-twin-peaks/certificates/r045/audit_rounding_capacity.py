from collections import defaultdict
from itertools import combinations, product


def profiles(n, lower, upper, prefix=()):
    if len(prefix) == n:
        yield prefix
        return
    choices = range(lower, upper + 1)
    if prefix:
        choices = range(max(lower, prefix[-1] - 1), min(upper, prefix[-1] + 1) + 1)
    for width in choices:
        yield from profiles(n, lower, upper, prefix + (width,))


def lower_round(profile, first_parity):
    return tuple(width - (width - first_parity - row) % 2 for row, width in enumerate(profile))


def grid_data(n):
    points = frozenset((x, y) for y in range(n) for x in range(n))
    neighbors = {
        (x, y): points & {(x-1, y), (x+1, y), (x, y-1), (x, y+1)}
        for x, y in points
    }
    return points, neighbors


def left_set(profile):
    return frozenset((x, y) for y, width in enumerate(profile) for x in range(width))


def component_test(left_a, left_b, points, neighbors):
    common_left = left_a & left_b
    common_right = points - (left_a | left_b)
    for difference in (left_a - left_b, left_b - left_a):
        remaining = set(difference)
        while remaining:
            current = {remaining.pop()}
            component = set(current)
            while current:
                current = set.union(*(set(neighbors[v]) for v in current)) & remaining
                remaining -= current
                component |= current
            if not any(neighbors[v] & common_left and neighbors[v] & common_right for v in component):
                return False
    return True


def exact_feasibility(left_pieces, points, neighbors, roots):
    left_root, right_root = roots
    assert all(left_root in left and right_root not in left for left in left_pieces)
    assert all(neighbors[left_root] <= left and not (neighbors[right_root] & left) for left in left_pieces)
    seen = {left_root, right_root}
    required_neighbors = {
        v: [neighbors[v] & (left if v in left else points - left) for left in left_pieces]
        for v in points - seen
    }
    while len(seen) < len(points):
        batch = {v for v in points - seen if all(local & seen for local in required_neighbors[v])}
        if not batch:
            return False
        seen |= batch
    return True


def audit_compatibility():
    compatibility_checks = 0
    class_checks = 0
    widths_checked = 0
    for n in range(4, 8):
        points, neighbors = grid_data(n)
        roots = ((0, 0), (n - 1, 0))
        options = list(profiles(n, 1, n - 1, (n // 2,)))
        sets = {profile: left_set(profile) for profile in options}
        groups = defaultdict(list)
        for h in options:
            g = lower_round(h, h[0] % 2)
            assert g[0] == h[0]
            assert all(abs(a-b) == 1 for a,b in zip(g, g[1:]))
            assert all(width - rounded in (0, 1) for width, rounded in zip(h, g))
            groups[g].append(sets[h])
            widths_checked += 1
        for group in groups.values():
            assert exact_feasibility(group, points, neighbors, roots)
            class_checks += 1
        if n <= 6:
            for a, b in combinations(sets.values(), 2):
                assert component_test(a, b, points, neighbors) == exact_feasibility([a, b], points, neighbors, roots)
                compatibility_checks += 1
    return widths_checked, class_checks, compatibility_checks


def count_by_transfer(n, lower, upper):
    counts = {h: int(h % 2 == 0) for h in range(lower, upper + 1)}
    for _ in range(n-1):
        counts = {h: counts.get(h-1, 0) + counts.get(h+1, 0) for h in counts}
    return sum(counts.values())


def count_by_walk_range(n, lower, upper):
    total = 0
    for steps in product((-1, 1), repeat=n-1):
        value = minimum = maximum = 0
        for step in steps:
            value += step
            minimum = min(value, minimum)
            maximum = max(value, maximum)
        a = lower - minimum
        b = upper - maximum
        total += max(0, b // 2 - (a - 1) // 2)
    return total


def audit_capacity():
    cases = 0
    for n in range(7, 15):
        root_band_width = (n + 3) // 4
        lower = root_band_width + 1
        upper = n - root_band_width - 1
        walk_count = count_by_transfer(n, lower, upper)
        assert walk_count == count_by_walk_range(n, lower, upper)
        upper_count = (upper // 2 - (lower - 2) // 2) * 2**(n-1)
        assert walk_count <= upper_count
        print(f"n={n}: confined all-origin parity walks={walk_count}, rounding upper={upper_count}")
        cases += 1
    return cases


if __name__ == "__main__":
    print("input curves / jointly feasible rounded classes / exact pair tests:", audit_compatibility())
    print("transfer-versus-range count cases:", audit_capacity())
