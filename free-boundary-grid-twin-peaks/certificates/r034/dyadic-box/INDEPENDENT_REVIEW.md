# Independent adversarial review: dyadic-box harmonic bootstrap

Scope: review of the stated r036 argument only.  No computation was run and
no source proof was modified.

## Result

The proposed deterministic inequality is valid, conditional only on the
already-audited reinsertion identity `R=b-l-1` and its harmonic image
inequality.  It yields the advertised lower bound for every one-peak word and
therefore `B/A1 -> infinity`; it does **not** settle the original two-peak
distance question.

## Checks

Let `v_k` be the first vertex at L-infinity shell `k`.  Before it appears no
vertex at shell at least `k` has appeared.  A neighbour of `v_k` that is
earlier is therefore in shell `k-1`.  If two coordinate offsets of `v_k-r`
have absolute value `k`, every neighbour remains in shell `k` or exits it, so
such a vertex cannot be first in its shell.  With exactly one active maximum
coordinate it has exactly one shell-`k-1` neighbour, denoted `p_k`.  Thus its
early-neighbour set is exactly `{p_k}`.  This also works for a clipped free
boundary box: clipping removes outward neighbours but creates no additional
strictly inward one.

Every neighbour of `p_k` has shell at least `k-2`, so `l_{p_k} >= T_{k-2}`;
`v_k` is a sole-parent child of `p_k`, so `b_{p_k} <= T_k`.  Hence
`R_{p_k} <= T_k-T_{k-2}-1`.  For `k=m+1,...,2m`, summation telescopes to
`T_{2m}+T_{2m-1}-T_m-T_{m-1}-m <= 2T_{2m}`.  The first hit of shell `2m`
occurs after at most all vertices of the strict L-infinity ball plus one, so
`T_{2m} <= (4m-1)^2+1`, even at a free boundary.  Consequently the claimed
`sum R <= 32m^2` and Cauchy lower bound `sum 1/R >= 1/32` hold.

For `m=2^j`, the selected parents lie in the mutually disjoint shell ranges
`m,...,2m-1`; for `k>=3` they are neither the root nor a root-neighbour.  The
count of eligible dyadic blocks is `(floor(log2 D)-1)_+`.  Finally
`D >= ceil((n-1)/2)` for every root coordinate, giving the stated uniform
bound.

No quantifier reversal, endpoint omission, or time-index off-by-one was found.
The conclusion remains only a divergence of `B/A1`; it supplies no direct
control of `B/A2`, of the L-change of measure, or of far-peak probability.
