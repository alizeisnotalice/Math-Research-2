# R21--22 rational interval audit

Audited commit: `7ae44feb6a38e9e69207b234d494af307138f546`.  This audit covers
the arithmetic and formula evaluation only, not the geometric/topological
claim that the cuts form a universal cover.

## Verdict

No arithmetic, Taylor-factorial, midpoint-remainder, square-root, arctangent,
or selected-root defect was found which changes the current conditional strict
formula conclusion

`V(1.5494 deg) < 0.8440935944`.

The exact replay of `certify_cut_area_formula.py` returned

```
[422046797191201692257571/500000000000000000000000,
 42204679719529451011179/50000000000000000000000]
```

whose upper-end margin to the target is
`470548988821/50000000000000000000000 > 0`.

## Checks made

* `I` uses exact `Fraction` endpoints; `rounded()` floors lowers and ceils
  uppers.  Division rejects intervals meeting zero.  `sqrt()` has the correct
  integer-square lower/upper tests (the extra upper unit is conservative).
* `rational_trig.py`: the Machin interval has the correct alternating-series
  endpoint direction.  Its sine degree-81/cosine degree-80 polynomial uses a
  conservative degree-82 Lagrange bound.  The one-step displacement `rad` is
  valid since both derivatives have absolute value at most one.
* Both 32-term arctangent implementations end at degree 63 and use the
  conservative `|z|^65/65` tail; each actual call passed the strict
  `|z| < 1/2` and positive-dot-product guards.  Thus the short-arc `atan2`
  branch is fixed.
* The jet stores derivatives divided by factorials.  The inverse and square
  root coefficient recurrences are the standard coefficient recurrences, and
  the density formula uses the coefficient of `Q_x Q_y' - Q_y Q_x'` through
  order two.  Its order-three jet is sufficient.  Exact scalar validation of
  `integral[-Delta/2,Delta/2] t^2 dt = Delta^3/12` is recorded in
  `interval-invariants.json`; hence the code's `f2/(12*N**3)` is correct when
  `f2` encloses `f''/2`.
* The additional exact interval audit confirms the signs/labels of the path
  circle/tangent roots (`G` right, `Z` upper, `T` the upper tangent, `Y`
  upper) and of the R22 roots: `S` is strictly lower than its mate, `T`
  strictly left and strictly closer to `J` than its mate, and `U` strictly
  upper.  It does not prove that this coordinate realization is the intended
  geometric construction.

## Actual runs and cross-checks

All runs used `/Users/zzh/.venvs/main/bin/python3.13` (Python 3.13.18,
`mpmath 1.3.0` where used), and exited zero.

* The 2048-cell AH rerun reproduced every mathematical field of the frozen
  `022-ah-area-refined.json`, including exactly
  `[20379591204838001029/250000000000000000000000,
    81518373003826865627/1000000000000000000000000]`.
  Byte strings differ only because the runtime-version field differs.
* An independent 60-digit mpmath adaptive integral returned
  `0.00008151836891801715773478188412529156023780697735642674438`, inside
  the rational 2048-cell interval.  This is a numerical cross-check only.
* A fresh 4096-cell rational run gave
  `[20379592100740497583/250000000000000000000000,
    40759184716134431033/500000000000000000000000]`, width
  `514653435867/500000000000000000000000`, and also contains that mpmath
  value.  It is a consistency refinement, not needed for the reported R22
  bound.

## Finding: replay dependency is incompletely bound (repair recommended)

`certify_cut_area_formula.py` accepts `022-ah-area-refined.json` after only
checking its `cells` field and the hash of `certify_ah_area.py`.  It neither
reruns the AH computation nor checks that the JSON interval equals its actual
output.  In addition, the AH calculation dynamically depends on
`certify_ah_path.py` and `rational_trig.py`, whose hashes are not recorded in
the AH JSON or checked by R22.  A hand-edited/foreign JSON with a current
`script_sha256` could therefore make the standalone R22 verifier accept a
spurious AH interval.

This did **not** alter this audit's result: fresh actual execution reproduced
the frozen mathematical fields.  Still, for a self-checking certificate,
have R22 rerun AH (or verify an output hash and equality of every required
field) and record/check a dependency manifest for the path and trig helpers.

## Files and SHA-256

Source scripts:

```
98f50387de75d0a331f053fbb8d194f738fd8c68cfc0dcfec2915a94eb192b9a  certify_ah_path.py
97e0076b911f51a9c532e14cacb9fb6e3ae7e68707a3c5af77048f9c6623e87f  rational_trig.py
a348d72d3905ead1f5a278f2c2c382e3a073e4730c018083039b65c415183fd3  certify_ah_area.py
ea3632bfdb6b5934da341d0be818935967daf9e9525ff2d0b2f3719fd835bc4a  certify_cut_area_formula.py
```

Audit artifacts:

```
4518eed1637c4a9556c6b70a3c739b17f4eeb55f1cc5c80db81557cf1546a557  ah-rerun-2048.json
83606adb7d09d5f508eef00f67031760f130e74a4900e912e0fab1fca620ad1d  ah-rerun-4096.json
035e241497fb4ee6a928bf9fdba35098e58c745e7c239581ec0cf6ec035a04b2  cut-rerun.json
b911382e5bd070e8b14bd534dabcc937d602aad0c8df7593bab2b8234da0e71c  path-rerun.json
97253e8d401fc397bf7bfb371ed72824e45696ba64e0649fc64c466ac2aa65ea  ah-independent-mpmath.json
d8b0f3e9049bde3aa32005ca489d86fb222c03d81c8d2afd9cee5689dbde76ce  validate_interval_invariants.py
9c3fc446703eb9ada687a021cab3d84ae509265f80dab5499202ace8ec31a7e6  interval-invariants.json
```
