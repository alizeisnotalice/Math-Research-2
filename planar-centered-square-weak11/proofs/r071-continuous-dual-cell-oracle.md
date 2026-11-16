# R71: exact all-position rectangle dual checking

2026-11-16. This appendix concerns a fixed finite collection of receiver
witnesses and valid capacity constraints. It supplies no new general bound
for the centered-square weak constant C2.

## The finite arrangement theorem

Let W_j be bounded axis-parallel rectangles, each explicitly open or closed,
and let C_i be bounded closed axis-parallel rectangles. Degenerate rectangles
are allowed: a closed point or segment is retained, whereas a rectangle open
in both coordinates with a degenerate side is empty. All endpoints and the
coefficients lambda_j, gamma_i are rational, with coefficients nonnegative.
Consider the desired inequality

    sum_j lambda_j 1_{W_j}(q) <= sum_i gamma_i 1_{C_i}(q).       (1)

Collect every x and y endpoint from all rectangles. Intersect the resulting
arrangement with the closed coordinate hull of the witnesses; endpoints of
that hull are included because they themselves are witness endpoints. Each
axis is a finite disjoint union of singleton points and open intervals
between consecutive endpoints. For a product cell choose the endpoint for
a singleton and the rational midpoint for an open interval. Every indicator
in (1) is constant on each product cell. This includes all zero-dimensional
corners and one-dimensional edges; no convention of ignoring null boundaries
is admissible for an atomic outside measure.

Consequently (1) holds everywhere in the witness hull if and only if it holds
at those finitely many representatives. Cells outside the witness union can
be omitted: the left side is zero there and the right side is nonnegative.
Points outside the witness hull likewise satisfy (1). Thus passing this exact
finite check certifies (1) for **every** q in the plane, and integrating against
any nonnegative outside measure nu gives

    sum_j lambda_j nu(W_j) <= sum_i gamma_i nu(C_i).            (2)

This is a precise extension of a dual over a finite candidate source grid.
Checking the dual only on that grid does not prove (1). A negative cell gives
an exact separating column: its representative q, the active witness/capacity
indices, and the rational slack

    s(q)=sum_i gamma_i 1_{C_i}(q)-sum_j lambda_j 1_{W_j}(q)<0.

The zero-dimensional and positive-dimensional minima are separately recorded.
The most negative two-dimensional cell is often a useful next LP column,
but the eventual certificate must check all dimensions. Adding a new source
column may require adding the corresponding bounding-box capacity menu; the
checker does not manufacture missing primal constraints.

## Permitted measure-zero exclusions

The default certificate is pointwise on the whole plane, including the source
locations S. Optionally the caller may assert nu({p})=0 for specified finitely
many points. Those coordinates are added to both axis arrangements, so each
excluded point is exactly a singleton product cell and cannot hide an open
region. The assertion must be justified by the mathematical model; the checker
does not infer it from a finite candidate grid.

An optional second exclusion lists supplied closed capacities C_i with **exact
right-hand side zero**. For any primal-feasible nonnegative measure,

    0 <= nu(C_i) <= 0,

so their union can be omitted when checking (2). The code rejects this exemption
if an id is missing or its exact supplied right-hand side is nonzero. It still
does not independently prove the supplied capacity is mathematically valid.
Passing with any exclusion yields (1) nu-almost everywhere under the declared
assumptions, rather than an unqualified pointwise theorem.

## Input and certificate

`computation/round071/continuous_dual.py` exports `check_dual(data)`. The CLI is

    python3 -B computation/round071/continuous_dual.py --input INPUT --output OUTPUT

The schema is `r071-continuous-dual-v1`, with these fields:

* `witnesses`: objects `{id,bounds,weight,boundary}`; boundary is `open` or
  `closed` and bounds are `[xmin,xmax,ymin,ymax]`.
* `caps`: the same objects, with boundary `closed`; optional `kind` distinguishes
  the intended bbox and residual constraints without changing their geometry.
* Coordinates, weights and optional `rhs` values are JSON integers or exact
  rational strings. JSON floating-point numbers and negative coefficients
  are rejected. A supplied negative capacity RHS is also rejected, because
  the nonnegative-measure primal would already be inconsistent.
* Optional `excluded_points` and `zero_capacity_cap_ids` have the meanings above.
* `metadata` is preserved in the input but not used to establish any inequality.

Every rectangle id must be globally unique. If RHS values are supplied for
every rectangle, the program also performs the exact arithmetic sum
`sum(weight*rhs)`. Witness RHS uses the caller's affine convention (in the
current shared LP it is mu_S(Q)-4r²); this arithmetic is **not** a new proof
of those RHS values or of the primal/dual signs.

The output contains the complete axis arrangement, every checked witness-union
cell, its rational representative, dimension, active ids and exact slack.
This is a finite certificate: given the input, a verifier can enumerate every
product cell, check its membership, and check all slacks. The script receipt
records Python/platform and SHA-256 of the precise input and checker. There is
no numerical tolerance, random seed or optimization algorithm in this checker.
`--summary-only` is available for intermediate column generation; archive a
full output for a final certificate. Zero weighted rectangles remain in the
arrangement and do not invalidate the proof.

## Actual regression and early-column evidence

The deterministic regression input is
`computation/inputs/r071-continuous-regression.json`; the runner is
`computation/round071/check_continuous_dual_regression.py`. Nine fixtures cover
open witness boundary, closed capacity boundary, overlapping rectangles,
closed residual R0, empty open degeneracy, positions outside a finite source
window, a closed singleton violation, a rational sliver of width 10^-18,
explicit point exclusion, and zero capacity collar exclusion (the boundary
and overlap checks share one fixture). Four additional invalid inputs are
rejected. Every returned representative is independently recomputed by a
direct membership loop rather than the oracle's cached axis bitmasks. Python
3.15.0 actually ran the fixture suite successfully; the complete output and
receipt are in `computation/results/r071-continuous-regression.json`.

Two early sparse duals from the compute worker were also actually checked:

* The single-receiver source-grid dual had affine constant -201/2500, but
  26 of 35 full-position cells had negative slack, with minimum -1. It does
  not certify an arbitrary-outside obstruction.
* The perimeter-twenty source-grid dual had constant -1/100, but 172 of 241
  cells had negative slack. A two-dimensional negative representative was
  (-3/10,47/200), slack -1/3, active witness ids 8095 and 8112, with no active
  positively weighted capacity from that sparse dual. It likewise does not
  certify an arbitrary-outside obstruction.

These are column-generation diagnostics, not mathematical counterexamples
to the desired shared-source method. Their exact input hashes and complete
cell outputs are in the two accompanying result files. The later integrated
round report must identify any final dual separately and check it with its
own exact input; early negative dual constants cannot be promoted to a result
for arbitrary outside positions.

For the six-active-witness subset the same all-position check was actually
run, giving 172 negative cells among 241, minimum -1/3, and the same
two-dimensional representative (-3/10,47/200), now with active witness ids
8074 and 8091. Its containing open cell is

    (-7/20,-1/4) x (23/100,6/25).

For S=[0,24/25]² this point is outside the zero capacity collar
[-1/25,1]², so the violation cannot be repaired merely by adding a zero cost
collar capacity. The corresponding full certificate table is
`computation/results/r071-six-active-continuous.json`.

The compute worker then added precisely this rational column, rebuilt all
bounding-box capacities and complete residual radius prefixes, and reported
an exact finite-support primal/dual optimum **+101/10000**, replacing the
previous -1/100. Its supplied file
`computation/results/round071-joint/one-column-result.json` was read here and
records 145 candidate outside positions and the added column. The solver
run and rational primal/dual verification are the compute worker's evidence;
this appendix does not claim an independent replay of that LP or of its
geometry. No further column search was launched. Thus the negative finite-grid
obstruction demonstrably needed this all-position check; it was not a valid
certificate for arbitrary outside supports. The integrated audit must separately
establish the precise actual-source consequence of the positive primal.
