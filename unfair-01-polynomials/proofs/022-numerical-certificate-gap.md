# R8--L9: feasibility review for a strict exclusion certificate

## Materials and scope

Read-only inputs reviewed:

- completed DS checkpoint `round021-R8-L9-ds/checkpoint.json`, SHA-256
  `7d45df172121371063c3f681f05f928424ab948e6e8be1ecbb2f7804c582938e`;
- exact-target/root input `secondM8-L9-inputs.json`, SHA-256
  `92bbaebd9ae148c1cb58bb1759a3cc9ce45efd96baf1d6f3cb354a84efc4b88d`;
- prepared option input `secondM8-L9-block-options.json`, SHA-256
  `81e2d55bb9736e17b023bf493c3e4fd95c408130af8468a335d41e026aecd85a`.

The checkpoint claims complete contiguous coverage of all
`607041386` complement-canonical balanced selections, and correctly records
the v3 source hash and input hashes.  It explicitly states that MPS DS ranking
is not a proof.  This review agrees with that boundary.

## What the completed numerical scan says

For each E, the table reports the first retained top64 row, independently
reconstructed on CPU float64 from individual root factors.  `score` is the
sum of distances outside `[0,1]` over the 78 nonendpoint coefficients of the
two factors; `min` and `max` are the CPU reconstructed extrema.

| E | top64 minimum score | CPU min | CPU max | score/78 |
|---:|---:|---:|---:|---:|
| 1 | 11.899330 | -1.313070 | 1.622834 | 0.152556 |
| 2 | 10.023004 | -1.513852 | 1.122530 | 0.128500 |
| 3 | 6.489216 | -0.965508 | 1.197993 | 0.083195 |
| 4 | 4.263879 | -0.609848 | 1.187025 | 0.054665 |
| 5 | 4.445771 | -0.650548 | 1.360330 | 0.056997 |
| 6 | 2.666022 | -0.829490 | 1.051004 | 0.034180 |
| 7 | 5.426912 | -0.949115 | 1.199150 | 0.069576 |
| 8 | 9.665456 | -0.998258 | 1.893836 | 0.123916 |

For E=9, code 0 is the expected fair mode: its reconstructed factors are
numerically `H(x)=1+x+x^2+x^3+x^4+x^5+x^8` and `H(x^9)`.  Its tiny observed
residuals are floating-point artifacts, not a candidate to exclude.  The
seven other E=9 modes have smallest stored score `9.962337`, hence stored
`score/78 = 0.127722`.

The CPU/GPU score discrepancy over retained rows is empirically at most about
`1.94e-6` for E=1--8 and `7.68e-6` in the E=9 top64 (the latter includes the
near-exact fair row).  This is useful engineering evidence only.

If a score lower bound were rigorous, score/78 would certify that at least one
coefficient has distance at least that amount from `[0,1]`.  The smallest
stored nonfair value is therefore numerically very well separated from zero.

## Why the checkpoint cannot itself be upgraded to a proof

It has no independently checkable all-domain margin ledger.  In particular it
does not retain, for every mixed-radix range, a lower bound on every score or a
witness coefficient and sign.  It retains only 64 GPU-selected rows per E.
The summary count `closed_unit_interval_fractional_masks=0` is likewise a GPU
classification count, not a directed interval statement.  Consequently:

1. a missed/incorrect GPU comparison, shader arithmetic error, or top64
   retention error could hide a smaller score outside the saved rows;
2. float64 root reconstruction of saved rows does not independently recompute
   all 607 million selections;
3. the decimal roots/options are explicitly `numerical_prepared_not_algebraically_certified`.

Thus **strict exclusion is not presently justified**, despite the encouraging
numerical margins.

## Coarse error-scale calculation

Let `u=2^-24` be float32 unit roundoff.  For each E and each side, I bounded a
selected four-block product in 1-norm by the product of the maximum option
1-norm in each block.  The worst side for E=3, for example, is
`M_3=1.1047e12`, hence `u^2 M_3=3.92e-3`.  If one had a *proved* DS forward
bound of the form

```
absolute coefficient error <= K u^2 M_E,
```

then comparison to the stored score/78 margins would require at least the
following upper limits on `K` for the worst side of each E:

| E | worst `u^2 M_E` | stored margin | required `K <` |
|---:|---:|---:|---:|
| 1 | 2.55e-6 | 0.152556 | 5.98e4 |
| 2 | 2.01e-3 | 0.128500 | 63.9 |
| 3 | 3.92e-3 | 0.083195 | 21.2 |
| 4 | 2.17e-3 | 0.054665 | 25.2 |
| 5 | 9.86e-5 | 0.056997 | 578 |
| 6 | 6.08e-5 | 0.034180 | 562 |
| 7 | 1.77e-5 | 0.069576 | 3.93e3 |
| 8 | 1.18e-7 | 0.123916 | 1.05e6 |

This is **not** a DS certificate: no such directed/error-constant theorem has
been supplied for the Metal shader, and it also omits error in numerically
prepared roots.  It quantifies the obstruction.  A loose but seemingly
"safe" multiplier such as `K=1e5` fails badly (E=3 would allow about 392,
versus stored margin 0.083).  A proof needs a tight interval calculation, not
a post-hoc appeal to double-single precision.

## Exact algebraic starting point and feasible certificate route

The exact mother polynomial is

```
H(x)=1+x+x^2+x^3+x^4+x^5+x^8.
```

I checked exactly that `gcd(H,H')=1` and, by Sturm count, H has zero real
roots.  Therefore H and H(x^9) have simple nonreal roots.  Numerical-only
diagnostics give minimum root separation about 0.5733 and minimum
`|H'(root)|` about 4.938; they support, but do not replace, interval isolation.

A realistic strict route is:

1. Use rational complex interval Newton/Krawczyk (or Arb-style validated
   complex balls) to isolate the four upper-half-plane roots of exact H, prove
   one simple root per disk, and form conjugate-pair quadratic coefficient
   balls.  Derive all ninth-root lifts with certified wedge/disc tests, rather
   than a nonvalidated complex power.
2. Prove from the exact target's binary coefficients and `A(0)=B(0)=1` that
   every nonnegative factor has all coefficients in `[0,1]` and leading
   coefficient one.  Then prove that real factors correspond to conjugate
   root-group choices and that the balance/canonical rules enumerate every
   permitted choice exactly once (apart from the stated complement convention).
3. Perform a separate deterministic **ball-arithmetic branch-and-bound** over
   the 607041386 canonical mixed-radix codes.  A node is discarded only when a
   coefficient ball for A or B lies wholly below zero or wholly above one.
   Record a Merkle/hashable range ledger: range, factor/witness index, rational
   interval, sign, input hash, precision, and coverage tree.  If a node cannot
   be discarded, split it or raise precision.  This is not a repeat GPU scan.
4. Verify code E=9/0 symbolically as `(H(x),H(x^9))`.  Every other leaf/range
   needs an interval witness.  A final verifier checks exact range coverage,
   no overlaps/gaps, and every witness inclusion.

## Concrete next proposition

The appropriate next target is the following fully checkable statement, not a
claim about the current checkpoint:

> For every canonical balanced root selection for R8--L9 other than E=9,
> code 0, certified coefficient balls contain a witness `c` satisfying
> `c <= -2^-7` or `c >= 1+2^-7`; E=9, code 0 is exactly
> `(H(x),H(x^9))`.

The numerical scan makes `2^-7` plausible (its smallest *stored* score/78 is
about 0.034), but only the proposed exact interval ledger can establish this
over the full selection set.  If that proposition is proved, nonnegativity and
the binary target immediately exclude every nonfair selection.
