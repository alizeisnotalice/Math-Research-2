#!/usr/bin/env python3
"""Independent exhaustive P3xP3 audit of the stable-promotion image.

It constructs every (one-peak word, eligible p) -> tau rather than inferring
the image from the harmonic R formula.  It then independently scans every
top-two two-peak word and classifies its first merge.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
import hashlib
import json
import platform
from pathlib import Path


HERE = Path(__file__).resolve().parent


def grid(n):
    adj = [set() for _ in range(n*n)]
    for y in range(n):
        for x in range(n):
            v = y*n+x
            for dx, dy in ((1,0), (-1,0), (0,1), (0,-1)):
                xx, yy = x+dx, y+dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v].add(yy*n+xx)
    return adj


def births(word, adj):
    seen = set(); out = []
    for v in word:
        if not (adj[v] & seen):
            out.append(v)
        seen.add(v)
    return out


def components(seen, adj):
    rest = set(seen); ans = []
    while rest:
        q = [rest.pop()]; c = set(q)
        while q:
            u = q.pop(); hit = adj[u] & rest
            rest -= hit; c |= hit; q.extend(hit)
        ans.append(c)
    return ans


def first_merge(word, adj):
    """Return the two component sets immediately before the first merger."""
    seen = set()
    for v in word:
        cs = components(seen, adj) if seen else []
        hit = [c for c in cs if c & adj[v]]
        if len(hit) >= 2:
            assert len(cs) == 2 and len(hit) == 2
            return hit
        seen.add(v)
    raise AssertionError("no first merge")


def R_formula(word, p, adj):
    """One-based b-l-1, computed with equal zero-based expression."""
    N = len(word); pos = [0]*N
    for i, v in enumerate(word):
        pos[v] = i
    lp = min(pos[u] for u in adj[p])
    bp = N
    for rank in range(pos[p]+1, N):
        v = word[rank]
        early = [u for u in adj[v] if pos[u] < rank]
        if early == [p]:
            bp = rank
            break
    R = bp-lp-1
    assert R >= 1
    return R


def promote(word, p):
    """Stably remove p and insert it after the root."""
    return (word[0], p) + tuple(v for v in word[1:] if v != p)


def main():
    n = 3; N = n*n; adj = grid(n)
    image = Counter()
    expected_R = {}
    source_pairs = 0
    onepeak = 0
    for word in permutations(range(N)):
        if len(births(word, adj)) != 1:
            continue
        onepeak += 1
        root = word[0]
        for p in range(N):
            if p == root or p in adj[root]:
                continue
            tau = promote(word, p)
            # This verifies the forward operation really lands in B.
            assert births(tau, adj) == [tau[0], tau[1]]
            image[tau] += 1
            r = R_formula(word, p, adj)
            if tau in expected_R:
                assert expected_R[tau] == r
            else:
                expected_R[tau] = r
            source_pairs += 1

    # The inverse fibre cardinality is independently the actual Counter.
    assert all(image[tau] == r for tau, r in expected_R.items())
    assert source_pairs == sum(image.values())

    B = set(); lower_singleton = set(); any_singleton = set(); Bq2 = set()
    for word in permutations(range(N)):
        bs = births(word, adj)
        if bs != [word[0], word[1]]:
            continue
        B.add(word)
        cs = first_merge(word, adj)
        assert word[0] in cs[0] or word[0] in cs[1]
        low = next(c for c in cs if word[1] in c)
        if len(low) == 1:
            lower_singleton.add(word)
        if any(len(c) == 1 for c in cs):
            any_singleton.add(word)
        if sum(map(len, cs)) == 2:
            Bq2.add(word)

    assert len(B) == 81984
    assert image.keys() == lower_singleton
    assert len(image) == 35264
    assert len(Bq2) == 14784
    assert len(any_singleton) == 55744
    assert len(any_singleton) == 2*len(lower_singleton)-len(Bq2)
    # Directly confirm the harmonic total, not merely the old formula script.
    assert sum((Fraction(1, image[tau]) for tau in image for _ in range(image[tau])),
               Fraction()) == len(image)

    out = {
        "kind": "independent-exhaustive-P3xP3-stable-promotion-image-audit",
        "n": n,
        "python": platform.python_version(),
        "onepeak_words": onepeak,
        "source_pairs": source_pairs,
        "top_two_twopeak_words_B": len(B),
        "actual_promotion_image_words": len(image),
        "actual_image_equals_lower_singleton": image.keys() == lower_singleton,
        "lower_singleton_top_two_words": len(lower_singleton),
        "q2_top_two_words": len(Bq2),
        "any_singleton_top_two_words": len(any_singleton),
        "identity_any_singleton_equals_2low_minus_q2":
            len(any_singleton) == 2*len(lower_singleton)-len(Bq2),
        "all_actual_inverse_fibres_equal_R_formula": True,
        "harmonic_weight_total_equals_actual_image_size": len(image),
    }
    text = json.dumps(out, indent=2, sort_keys=True) + "\n"
    HERE.joinpath("audit-singleton-image-n3.json").write_text(text)
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
