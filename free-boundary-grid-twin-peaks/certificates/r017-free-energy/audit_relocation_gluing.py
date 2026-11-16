#!/usr/bin/env python3
"""Independent finite audits for root relocation, interface gluing, and tree scales."""

from itertools import combinations, permutations
from math import ceil, comb


def vertices(a, b):
    return tuple((x, y) for y in range(b) for x in range(a))


def adjacency(a, b):
    V = vertices(a, b)
    return {v: {(v[0] + dx, v[1] + dy) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
                if (v[0] + dx, v[1] + dy) in V} for v in V}


def onepeak(word, adj, root=None):
    seen = set()
    for j, v in enumerate(word):
        if j == 0:
            if root is not None and v != root:
                return False
        elif not (adj[v] & seen):
            return False
        seen.add(v)
    return True


def path_to_corner(root):
    x, y = root
    # Deterministic shortest path: horizontal first, then vertical.
    return tuple((xx, 0) for xx in range(x + 1)) + tuple((x, yy) for yy in range(1, y + 1))


def relocate(word):
    path = path_to_corner(word[0])
    P = set(path)
    return path + tuple(v for v in word if v not in P)


def audit_relocation(a, b):
    V, adj = vertices(a, b), adjacency(a, b)
    fibres = {}
    source = 0
    for w in permutations(V):
        if onepeak(w, adj):
            source += 1
            image = relocate(w)
            assert onepeak(image, adj, root=(0, 0))
            fibres[image] = fibres.get(image, 0) + 1
    m = a * b
    assert max(fibres.values()) <= m ** (a + b)
    return {'rectangle': [a, b], 'source_onepeak_words': source,
            'corner_images': len(fibres), 'max_fibre': max(fibres.values()),
            'proved_fibre_bound': m ** (a + b)}


def onepeak_root_words(a, b, root):
    V, adj = vertices(a, b), adjacency(a, b)
    return [w for w in permutations(V) if onepeak(w, adj, root=root)]


def audit_interface_2x2():
    # Left R_{2,2} and right R_{2,2}; roots (1,0) and (2,0) touch.
    left = onepeak_root_words(2, 2, (1, 0))
    right0 = onepeak_root_words(2, 2, (0, 0))
    right = [tuple((x + 2, y) for x, y in w) for w in right0]
    V, adj = vertices(4, 2), adjacency(4, 2)
    outputs = set()
    for u in left:
        for v in right:
            for slots in combinations(range(1, 8), 3):  # first global point is u[0]
                slots = set(slots)
                iu, iv = 1, 0
                word = [u[0]]
                for t in range(1, 8):
                    if t in slots:
                        word.append(u[iu]); iu += 1
                    else:
                        word.append(v[iv]); iv += 1
                word = tuple(word)
                assert onepeak(word, adj, root=(1, 0))
                assert word not in outputs
                outputs.add(word)
    expected = comb(7, 3) * len(left) * len(right)
    assert len(outputs) == expected
    return {'left_interface_corner_words': len(left),
            'right_interface_corner_words': len(right),
            'good_interleavings_per_pair': comb(7, 3),
            'outputs': len(outputs), 'exact_injection_verified': True}


def audit_quadtree_coefficients():
    checked = 0
    for q in range(1, 513):
        J = (q - 1).bit_length()
        for ell in range(J):
            t = J - ell
            # ceil(q/2^ell) <= 2^t and 4^J < 4q^2.
            assert ceil(q / (2 ** ell)) <= 2 ** t
            assert 4 ** ell * ceil(q / (2 ** ell)) * 2 ** t <= 4 * q * q
            checked += 1
    return {'q_range': [1, 512], 'level_coefficients_checked': checked,
            'dyadic_reindexing_bound_verified': True}


if __name__ == '__main__':
    print({'relocation': [audit_relocation(2, 2), audit_relocation(3, 2)],
           'interface': audit_interface_2x2(),
           'quadtree': audit_quadtree_coefficients()})
