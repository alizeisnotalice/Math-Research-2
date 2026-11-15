#!/usr/bin/env python3
"""Verify a Bezout identity by standard-library rational polynomial arithmetic."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from verify_positive_certificate import add, mul, scale, const, variable

I = [0, 2, 4, 7]
J = [0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 15, 17, 20]
one = const(1, 2); p = variable(0, 2); q = variable(1, 2)
A = {0: one, 2: p, 4: q, 7: one}
B = {0: one, 2: add(one, p, -1), 9: p,
     15: add(one, p, -1), 17: add(one, q, -1), 20: one}
B[4] = add(add(one, q, -1), mul(p, B[2]), -1)
B[6] = add(add(one, mul(p, B[4]), -1), mul(q, B[2]), -1)
B[8] = add(add(one, mul(p, B[6]), -1), mul(q, B[4]), -1)
B[10] = add(add(one, mul(p, B[8]), -1), mul(q, B[6]), -1)
B[11] = add(add(one, B[4], -1), mul(p, B[9]), -1)

def residual(k):
    answer = scale(one, -1)
    for i in I:
        if k-i in J:
            answer = add(answer, mul(A[i], B[k-i]))
    return answer

for k in [2, 4, 6, 8, 9, 10, 11, 22, 24]:
    assert residual(k) == {}
f = [residual(k) for k in [13, 15, 17]]
# Numerators from exact rational linear algebra, independently expanded here.
D = 21979860606
numerators = [
 {(4,0):9112363463,(3,1):-10697186667,(3,0):-47527989905,
  (2,2):-19566641178,(2,1):-25361243740,(2,0):8384766323,
  (1,2):26777520321,(1,1):164844570717,(1,0):-10782917857,
  (0,3):47518985718,(0,2):-6357891314,(0,1):8499542887,(0,0):29823847673},
 {(2,1):12233018508,(2,0):-78123737593,(1,2):13976172270,
  (1,1):-3409169126,(1,0):-14845469524,(0,2):-12578555043,
  (0,1):57971452097,(0,0):5339275731},
 {(2,0):-27337090389,(1,1):44324578509,(1,0):55347868659,
  (0,2):72676095804,(0,1):10472841057,(0,0):16640584875}]
Ps = [{m: Q(c, D) for m, c in P.items()} for P in numerators]
answer = {}
for P, F in zip(Ps, f):
    answer = add(answer, mul(P, F))
assert answer == one

def encode(poly):
    return [{'powers': list(m), 'coefficient': [c.numerator, c.denominator]}
            for m, c in sorted(poly.items())]
print(json.dumps({'I': I, 'J': J, 'field': 'Q', 'variables': ['p=a2', 'q=a4'],
                  'eliminated_diagonals': [2,4,6,8,9,10,11,22,24],
                  'residual_diagonals': [13,15,17],
                  'f': [encode(F) for F in f], 'P': [encode(P) for P in Ps],
                  'exact_sum_P_f': encode(answer),
                  'conclusion': 'No complex solution of the normalized support equations',
                  'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
