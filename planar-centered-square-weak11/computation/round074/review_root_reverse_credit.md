# R74 review of the root reverse-credit numerical scope

Status: PASS for the stated exploratory finite class and the exact selected
two-atom equality check. No1200-case area replay or general theorem is claimed.

The reviewed root script is `.local/round074-root/probe_reverse_credit.py`,
SHA8f721f9a0cb5816131eb9f2e242978f96a6c69b5a66c7631361186ff8232b5f9.
Its saved output and receipt hashes agree. The imported complete-closure/
area engine has SHA3305e4006ff17476b235da0d1c4bd4b20ab4e9a883a30b78322378abbe174671.
The supplement manifest pins all four actual source artifacts.

The script tests positive finite measures supported in the **closed** square
Q=[-1/2,1/2]^2, with no atom at its center. Counts2,3,4,6,8,12 each have200
fixed-seed draws. Positions lie on the1/20 lattice, using three modes:
the full square, its boundary sides, and four corner subregions. Integer
weights1,...,20 are scaled along their common weight ray. This is neither
exhaustive position enumeration nor the entire center-feasible mass polytope.

For integer distance d and cumulative raw weight W_d, the chosen scale is

    min_d 4d^2/(400 W_d).

Thus every complete closed L-infinity distance prefix has mass at most
4(d/20)^2. Ties are grouped before cumulation, and at least one prefix
saturates. Since no atom is at0, these exact rational conditions are
equivalent to0 being outside the strict halo. The largest-distance row
also gives total mass<=1. The normalization is correct.

Calling `geom.closures(pts,[0]*count)` is safe here: that function returns
every nonempty coordinate-bbox closure and does not discard rows marked
ineligible. The probe ignores its zero-mass eligibility flag and computes
actual rational mass per row. Its comparison

    integer_span_squared < 400*actual_mass

is exactly physical_span_squared<actual_mass. Equality is correctly excluded.
For each eligible closure the receiver rectangle has endpoints
maxsource-sqrt(mass)/2 and minsource+sqrt(mass)/2 in both axes. Their union
is the complete strict halo for a finite positive source measure. The clip
to Q and call to the union-area routine therefore represent the intended
localized halo, but square roots, endpoints, merging and areas are **floating**.
There is no directed error bound, so the ratios only screen this finite sample.
The printed key `number_above_4_plus_1e6` uses the actual expression4+1e-6;
its label omits the minus sign, a cosmetic ambiguity rather than a change
to the tested expression. Negative search results do not prove coefficient4.

## Exact top-ratio case

Only seed742090 is independently reconstructed. It gives raw weights8,7,
scale1/15, and atoms

    (-1/2,-1/2) of mass8/15,
    ( 1/2, 1/2) of mass7/15.

There is one center distance prefix, t=1/2, with total mass1=4t^2;
the center is outside the strict halo. The complete nonempty subset
catalog has three masks. The pair has mass=span_squared=1, so it is
ineligible. Only the two singleton halos contribute. Their squared radii
are2/15 and7/60, both below1/4, so the two squares are disjoint. Inside
closed Q, each contributes exactly its corner quarter, area m_i/4.
Consequently

    nu(Q)=1,  |H(nu) intersect Q|=1/4,  ratio=4.

The root float ratio4.000000000000001 is a rounding overshoot, not a strict
counterexample. The square's boundary matters for the numerator: with
open Q both atoms would be excluded and nu(Q)=0. Lebesgue clipped area
does not change on that boundary, but this particular equality argument
must retain closed-Q mass semantics.

The deterministic checker uses Python3.13.18/Fraction, independently
reconstructs only this seed, checks all three subsets and the complete
center prefix, and verifies the exact quarter areas. It exited0; PID94189
was verified absent. Source receipt/status inspection is kept distinct
from independently recomputing all sample halo areas. No agent, MATLAB,
optimization, general reverse-credit proof, or Borel conclusion was added.

Reproduce with the root-owned three probe artifacts in one directory:

    python computation/round074/check_reverse_credit_diagonal.py --probe-root /absolute/probe/directory --out computation/results/round074-reverse-credit-review/diagonal-exact-check.json

All earlier R74 main/preliminary and R73 bundles remain unchanged.
