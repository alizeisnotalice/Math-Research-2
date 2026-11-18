"""Finite nearest-neighbour calibration of the ordered hitting determinant."""
from collections import defaultdict
from fractions import Fraction
from itertools import product, permutations
import json

STARTS = (-2, 0, 2)
HEIGHT = 8

def barrier(x):
    return min(5, abs(x - 1) // 2 + int(x >= 3))

def block(x):
    return 0 if x < -1 else (1 if x < 2 else 2)

def transition(x, time, step):
    return x if time <= barrier(x) else x + step

def endpoint_counts(start):
    law = {start: 1}
    for time in range(HEIGHT, 0, -1):
        after = defaultdict(int)
        for x, weight in law.items():
            for step in (-1, 1):
                after[transition(x, time, step)] += weight
        law = after
    return [sum(weight for x, weight in law.items() if block(x) == j)
            for j in range(3)]

def ordered_path_count():
    law = {STARTS: 1}
    for time in range(HEIGHT, 0, -1):
        after = defaultdict(int)
        for positions, weight in law.items():
            for steps in product((-1, 1), repeat=3):
                targets = tuple(transition(x, time, step)
                                for x, step in zip(positions, steps))
                if targets[0] < targets[1] < targets[2]:
                    after[targets] += weight
        law = after
    return sum(weight for pos, weight in law.items()
               if tuple(map(block, pos)) == (0, 1, 2))

A = list(map(endpoint_counts, STARTS))
det = 0
for p in permutations(range(3)):
    inversions = sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
    det += (-1)**inversions * A[0][p[0]] * A[1][p[1]] * A[2][p[2]]
ordered = ordered_path_count()
assert det == ordered and det > 0
print(json.dumps(dict(height=HEIGHT, starts=STARTS, counts=A,
                      determinant_count=det, ordered_count=ordered,
                      probability=str(Fraction(det, 2**(3*HEIGHT)))),indent=2))
