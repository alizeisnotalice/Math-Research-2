from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb, prod
from random import Random
import sys


def neighbors(n):
    return {
        (x, y): frozenset(
            (u, v)
            for u, v in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1))
            if 0 <= u < n and 0 <= v < n
        )
        for x in range(n)
        for y in range(n)
    }


def higher_neighbors(word, graph):
    seen = set()
    higher = {}
    for vertex in word:
        higher[vertex] = graph[vertex] & seen
        seen.add(vertex)
    return higher


def cut_edges(graph, left):
    return sum(len(graph[vertex] - left) for vertex in left)


def theta(higher, left):
    roots = [vertex for vertex, choices in higher.items() if not choices]
    if len(roots) != 2 or ((roots[0] in left) == (roots[1] in left)):
        return Fraction(0)
    return prod(
        Fraction(sum((parent in left) == (vertex in left) for parent in choices), len(choices))
        for vertex, choices in higher.items()
        if choices
    )


def all_parent_cut_counts(word, higher):
    choices = [tuple(higher[v]) or (None,) for v in word]
    high_root = word[0]
    counts = Counter()
    for assignment in product(*choices):
        root_of = {}
        for vertex, parent in zip(word, assignment):
            root_of[vertex] = vertex if parent is None else root_of[parent]
        counts[frozenset(vertex for vertex in word if root_of[vertex] == high_root)] += 1
    return counts


def audit_parent_words(n, words):
    graph = neighbors(n)
    examined = 0
    total = 0
    for word in words:
        higher = higher_neighbors(word, graph)
        roots = [vertex for vertex in word if not higher[vertex]]
        if len(roots) != 2:
            continue
        free_vertices = [vertex for vertex in word if vertex not in roots]
        counts = all_parent_cut_counts(word, higher)
        assignments = sum(counts.values())
        total += assignments
        examined += 1
        summed = Fraction(0)
        for bits in product((False, True), repeat=len(free_vertices)):
            left = frozenset([roots[0]] + [v for v, choose in zip(free_vertices, bits) if choose])
            actual = theta(higher, left)
            assert actual == Fraction(counts[left], assignments)
            assert actual <= Fraction(3, 4) ** cut_edges(graph, left)
            summed += actual
        assert summed == 1
    return examined, total


def sharp_word(n, k):
    left_core = [(x, y) for x in range(k - 1) for y in range(n)]
    right_core = [(x, y) for x in reversed(range(k + 1, n)) for y in range(n)]
    interface_early = [(k - 1 + y % 2, y) for y in range(n)]
    interface_late = [(k - y % 2, y) for y in range(n)]
    return left_core + right_core + interface_early + interface_late


def audit_sharp_cuts():
    examined = 0
    for n in range(4, 10):
        graph = neighbors(n)
        for k in range(2, n - 1):
            word = sharp_word(n, k)
            for axis in (0, 1):
                rotated = [(y, x) for x, y in word] if axis else word
                assert len(rotated) == len(set(rotated)) == n * n
                higher = higher_neighbors(rotated, graph)
                roots = [vertex for vertex in rotated if not higher[vertex]]
                assert len(roots) == 2
                assert sum(abs(a - b) for a, b in zip(*roots)) == n - 1
                left = frozenset(vertex for vertex in graph if vertex[axis] < k)
                assert graph[roots[0]] <= left
                assert not (graph[roots[1]] & left)
                assert cut_edges(graph, left) == n
                assert theta(higher, left) == Fraction(2, 3) ** 2 * Fraction(3, 4) ** (n - 2)
                examined += 1
    return examined


def profiles(n):
    def visit(prefix):
        if len(prefix) == n:
            yield prefix
            return
        for k in range(max(1, prefix[-1] - 1), min(n - 1, prefix[-1] + 1) + 1):
            yield from visit(prefix + (k,))

    for k in range(1, n):
        yield from visit((k,))


def audit_profile_counts():
    examined = 0
    for n in range(2, 8):
        graph = neighbors(n)
        counts = Counter()
        for profile in profiles(n):
            left = frozenset((x, y) for x, y in graph if x < profile[y])
            jumps = sum(a != b for a, b in zip(profile, profile[1:]))
            assert cut_edges(graph, left) == n + jumps
            counts[jumps] += 1
            examined += 1
        for jumps, count in counts.items():
            assert count <= (n - 1) * comb(n - 1, jumps) * 2**jumps
    for width, tilt in ((20, 32), (14, 20)):
        base_rational_part = Fraction(3, 4) * (1 + Fraction(3, 2 * tilt))
        base_to_power = tilt * base_rational_part**width
        assert base_to_power < 1
        for n in range(2, 201):
            cutoff = n // width
            weighted_count = sum(comb(n - 1, j) * Fraction(3, 2) ** j for j in range(cutoff + 1))
            assert weighted_count <= tilt**cutoff * (1 + Fraction(3, 2 * tilt)) ** (n - 1)
        print(f"tilted base^{width} = {base_to_power}; cutoff checks n=2..200")
    return examined


def audit_corridors():
    examined = 0
    for n in range(5, 10):
        for center in range(2, n - 2):
            graph = neighbors(n)
            left = [(x, y) for x in range(center) for y in range(n)]
            right = [(x, y) for x in reversed(range(center + 1, n)) for y in range(n)]
            corridor = [(center, y) for y in range(n)]
            higher = higher_neighbors(left + right + corridor, graph)
            assert [vertex for vertex, choices in higher.items() if not choices] == [left[0], right[0]]
            assert graph[left[0]] <= set(left)
            assert graph[right[0]] <= set(right)
            assert right[0][0] - left[0][0] == n - 1
            total = Fraction(0)
            for bits in product((False, True), repeat=n):
                left_basin = frozenset(left + [v for v, choose in zip(corridor, bits) if choose])
                mass = theta(higher, left_basin)
                assert mass > 0
                total += mass
                examined += 1
            assert total == 1
    return examined


def main():
    print(sys.version)
    print("n=2 all words and their K=2 parent selections:", audit_parent_words(2, permutations(neighbors(2))))
    rng = Random(20261118)
    points = list(neighbors(3))
    sample = [tuple(rng.sample(points, len(points))) for _ in range(200)]
    print("n=3 seeded sample of 200 words, K=2 and their parent selections:", audit_parent_words(3, sample))
    print("n=4..9 sharp internal cuts, both axes:", audit_sharp_cuts())
    print("n=2..7 complete Lipschitz profiles:", audit_profile_counts())
    print("n=5..9 width-one corridor cuts:", audit_corridors())


if __name__ == "__main__":
    main()
