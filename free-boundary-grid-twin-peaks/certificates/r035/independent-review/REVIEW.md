# r035 independent mathematical review

Scope: read-only review of the four staged proof packages.  The conclusions
reviewed below remain auxiliary statements and do not settle the original
twin-peaks question.

## Result

No genuine mathematical or evidence gap was found in the stated ranges.

* In `gap-mixture`, root swapping is a fixed-point-free involution on the
  top-two words because the first two roots are nonadjacent and every suffix
  sees the same set.  Formula (7) normalizes correctly; the `q=2` and one-
  singleton cases are separately and correctly handled.  The fixed-
  `t` approximation uses `q -> infinity` in B-law, not an illicit infinite
  sum interchange.  Endpoint atoms account for exactly half their mass in
  finite-gap limits, yielding the tightness criterion stated.  This says
  nothing about uniform integrability or `A2/B`.
* In `weighted-reinsertion`, inserting after the first low-root neighbour
  makes its rank `q`, hence the fibre weight is `(q-1)/R`; the `q-1` possible
  low-root ranks count the unrestricted `B_+` fibre, not merely top-two
  words.  The rank/label conversion gives `D=R` and denominator `N-j-1` as
  written.  The strengthened dyadic bound is only conditional on the actual
  shell time and is not misused as an area lower bound.
* In `completion-transfer`, the reverse stable move has at most `M` slots,
  so the direction of (6) is correct.  `omega_k` has denominators `M+j`,
  giving `omega_k <= k/M`; the revised factor and the algebra in (10) have
  the correct inequality direction.  The conclusion is explicitly averaged,
  not pointwise.
* In the odd comb lemma (`n>=5`), the selected range `2<=j<=L-1` excludes
  root/root-neighbour hazards.  Delayed odd-row vertices cannot preempt the
  horizontal sole-child, so (2) survives arbitrary tooth interleaving.  The
  toothwise telescope and the `T=h(L+2)-1` identity are correct; the result is
  properly confined to this comb family.

The existing finite scripts are suitable index checks only; the general claims
rest on the displayed fibre/transfer arguments, not on those experiments.
