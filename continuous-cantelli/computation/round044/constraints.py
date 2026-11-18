"""Exact finite-constraint Jacobian and the second-derivative control."""
from fractions import Fraction as Q
from itertools import permutations
import json

# Nodes x=log(3),log(4),log(5) belong to (1,2); parameters 0,1,2.
def entry(t, k):
    return (Q(k)**t + Q(k)**(-t))/2

a = [[entry(t, k) for k in (3,4,5)] for t in range(3)]
det = Q(0)
for perm in permutations(range(3)):
    inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1,3))
    det += (-1)**inversions * a[0][perm[0]] * a[1][perm[1]] * a[2][perm[2]]
assert det > 0
assert 12*Q(1,24)-1 == -Q(1,2) and 12*Q(1,8)-1 == Q(1,2)
print(json.dumps(dict(evaluation_determinant=str(det), curvature_margins=["-1/2","1/2"])))
