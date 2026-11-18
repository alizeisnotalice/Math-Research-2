"""R72 exact low-rank positive coverage tests in dimension three."""
import functools
import itertools
import json
import platform
import random
from check import span, inside, dot


@functools.cache
def plus(W, v):
    return W | frozenset(z ^ v for z in W)


def main():
    X = range(8)
    spaces = {span(J) for k in range(4) for J in itertools.combinations(range(1, 8), k)}
    normal = {}
    for W in spaces:
        for y in X:
            mask = sum(1 << x for x in X if inside(x, y, W))
            normal[mask] = W

    def partitions(mask):
        if not mask:
            yield ()
            return
        first = mask & -mask
        for flat in sorted(normal):
            if flat & first and flat & mask == flat:
                for tail in partitions(mask ^ flat):
                    yield (flat, *tail)

    all_parts = list(partitions(255))
    rng = random.Random(72027)
    families = [[(f, 1) for f in part] for part in all_parts]
    for _ in range(1000):
        a, b = rng.choices(all_parts, k=2)
        fl = [(f, 3) for f in a] + [(f, 5) for f in b]
        if rng.randrange(2):
            fl.append((rng.choice(sorted(normal)), 4))
        families.append(fl)
    anchors = 0
    for weighted in families:
        D = [sum(w for F,w in weighted if F & (1 << x)) for x in X]
        U = span(v for v in X if sum(d*(-1)**dot(x,v) for x,d in enumerate(D)))
        S = span(v for F,_ in weighted for v in normal[F])
        target = span([*S,*U])
        for x in X:
            reached = {U}
            for flat,_ in weighted:
                if flat & (1 << x):
                    continue
                reached |= {plus(W,v) for W in reached for v in normal[flat]}
                if target in reached:
                    break
            assert target in reached
            anchors += 1
    print(json.dumps(dict(python=platform.python_version(), seed=72027,
        dimensions=3, all_unweighted_partitions=len(all_parts),
        mixed_weighted_covers=1000, covers=len(families), anchors=anchors,
        check="actual spanning bases with at most one representative per inactive indexed flat"),indent=2))


if __name__=='__main__':
    main()
