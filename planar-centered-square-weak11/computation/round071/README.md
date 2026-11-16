# R71 actual shared outside-mass computation

Actual run date: 2026-11-16, Asia/Shanghai. Python 3.13.18,
NumPy 2.5.4, SciPy 1.18.1 (HiGHS via `linprog`). No MATLAB was
started by this computation agent. Root owns the independent MATLAB runs.

This component fixes S to the 5 by 5 grid `(6i/25,6j/25)`,
`0<=i,j<=4`, with every mass `1/25`. Its total mass is 1 and side is
`24/25`. The initial possible outside positions are the 13 by 13 grid
of the same step, indices -4 through 8, excluding S: 144 positions.
Every LP uses the **same** outside-mass variables for all receivers.

## Exact finite model

For every fixed-position coordinate-bbox closure B, impose
`nu(B) <= max(1,d(B)^2)-muS(B)`. These are the R64 exact finite
mass-polytope constraints for L<=1; fixing S gives L=1. For each
receiver x, impose `nu(closedQ(x,d)) <= 4d^2` at **every distinct
outside-source L-infinity distance d**. These complete distance-prefix
constraints exclude x from the strict residual halo at all radii.
The designated full witness has fixed radius r and the additional row
`-nu(closedQ(x,r))+eta <= muS(closedQ(x,r))-4r^2`.
Maximize eta with nonnegative outside masses and free eta.

The LP captures at `distance <= r`, including boundary atoms. A
strictly positive rational eta permits slightly enlarging the radius
to obtain a genuine open-square full witness. No closed/open convention
is silently exchanged in either the finite model or continuous oracle.
The historical main-driver JSON has an old top-level `scope` template
saying radius 49/100 even for menu/25 models. That string is corrected
explicitly by `scope-addendum.json`; per-model `radius_hundredths`,
exact witness rows and MAT matrices determine the actual instance.
The final code now writes a dynamic radius list, verified by one new
r=1/4 metadata smoke solve. Historical raw results are not rewritten.

Integer geometry in the original grid uses hundredths and RHS
denominator 10000. The augmented model uses units 1/200 and denominator
40000. Floating LP proposals are rationalized and then every cap,
residual prefix and witness is checked by `Fraction`. The dual
weights are rationalized, normalized so witness weights sum to 1,
and any negative finite column deficit is covered by a valid singleton
cap. Final exact primal and dual bounds coincide in the saved models.

## Bounded probes and finite obstruction

There were 66 actual LP solves, including the initial 14 large-radius
singletons, a 32-model explicit radius menu, diagnostic duplicates,
eight prescribed shared combinations, and seven dual-guided six-point
or omission models, and the metadata smoke repeat. These are itemized in `probe-summary.json`; no
random search or unreported large batch was run. The r=.49 singletons
were individually infeasible, and do not provide a shared-mass loss.
The explicit menu r=.15,.25,.37,.49 supplied positive singleton models.

At fixed r=1/4, the three D4 receiver orbits used below have exact
singleton optima 151/2500, 7/100, and 7/100, represented by
`(-1/10,0)`, `(-1/10,6/25)`, `(-1/10,12/25)`.
`three-orbit-singletons.json` retains the three exact primal/dual
certificates. Symmetries of both S and the outside grid transfer these
certificates to all six selected receivers:

    (-1/10,6/25), (-1/10,12/25), (-1/10,24/25),
    (0,53/50), (12/25,53/50), (18/25,53/50).

For these six receivers jointly, the initial 144-position model has
exact optimum eta=-1/100. The dual has weight 1/6 on each witness and
on two bbox caps with residual capacities 9/25 and 11/25; residual-halo
rows have zero dual weight. The six witness RHS sum to -43/50, so the
dual value is `(9/25+11/25-43/50)/6=-1/100`. Omitting receivers in the
listed order gives exact optima `0,1/300,0,0,1/300,0`. This is a finite
position obstruction, not a theorem for arbitrary outside positions.

The continuous exact indicator oracle, implemented and run by the
separate continuous-dual agent, finds uncovered non-grid columns.
For the six-point dual its open 2D negative-cell representative is
`q=(-3/10,47/200)` with column slack `-1/3`. It lies outside the closed
collar `[-1/25,1]^2`. The separately retained single-negative,
perimeter-twenty and six-active dual-geometry JSON files are the exact
inputs of that oracle. Its independent outputs are integrated by root.

## One genuine extra position reverses the sign

Adding only that q, rebuilding all finite bbox caps and every residual
distance prefix, and retaining four globally valid side-one zero-cap
collar squares gives 145 possible outside positions and 7940 rows.
The exact new optimum is eta=101/10000>0. Expansion stopped here.
The selected rational primal has only eight positive outside atoms:

| position | mass |
| --- | --- |
| (-6/25,0) | 599/10000 |
| (-6/25,18/25) | 599/10000 |
| (-6/25,6/5) | 601/5000 |
| (6/25,6/5) | 599/10000 |
| (12/25,6/5) | 9/5000 |
| (18/25,6/5) | 49/625 |
| (24/25,6/5) | 599/10000 |
| (-3/10,47/200) | 401/5000 |

Thus the actual positive input has 33 atoms and total mass 7601/5000.
The actual-support check enumerates 456 distinct complete bbox
closures, of which 177 are eligible. L=1 with two maximum-mass ties;
all atoms lie in one eligible component. Actual span is 63/50 and
`dA^2-(35/33)mu(A)=-371/15000`: this example is in the R53 terminal
size class. It is not a counterexample in the remaining nonterminal
R63 selection domain. The actual 33-atom input, without zero-mass
candidate positions, is `positive33-actual-input.json`.

All six full closed r=1/4 witnesses capture mass 2601/10000 and have
margin 101/10000. They also have genuine open witnesses, for example
the common radius 251/1000. The first five receivers have strictly
positive residual-prefix slack at every positive prefix. The sixth
has equality at radius 7/50 and residual mass 49/625.

## A small actual deletion-region certificate

For a receiver x, let d be a positive residual distance event, M_d its
closed-prefix mass and delta_d=4d^2-M_d. For any epsilon satisfying

    epsilon <= min(d/2, delta_d/(16d)) for every event,

one has `M_d < 4(d-epsilon)^2` when every delta_d>0. If y lies in the
closed L-infinity epsilon-box about x, any residual subset captured at
y has its x-distance at most its y-radius plus epsilon; taking its
farthest x-distance event proves `y notin H(nu)`. This keeps the same
outside mass vector at every receiver and every nearby point.

For each of the first five receivers, epsilon=101/40000 passes all
these exact inequalities. The common **open** full radius 253/1000
is strictly greater than `1/4+epsilon`, and

    2601/10000 - 4(253/1000)^2 = 254/62500 > 0.

Therefore the five closed epsilon-boxes lie in the actual deletion
region `H(muS+nu) \ H(nu)`. They are pairwise disjoint and have total
area at least `10201/80000000 = .0001275125`. The sixth point is not
assigned a symmetric box by this strict-prefix lemma. No full halo
or deletion area was calculated, and this small lower bound supplies
no new lower bound for C2.

## Reproduction and retained scope

Run from the topic root using `/Users/zzh/.venvs/main/bin/python3`, with
TMPDIR/TMP/TEMP and cache directed to `/Volumes/data/CodexTemp/T2`:

    python computation/round071/six_receiver_check.py
    python computation/round071/joint_outside_lp.py --out computation/results/round071-joint/single-positive.json --centres '[[-10,24]]' --radius 25 --export computation/inputs/r071/single-positive
    python computation/round071/augment_support.py --stem computation/results/round071-joint/one-column
    python computation/round071/primal_geometry_certificate.py

MAT schemas are `A`, `b`, `f`, `lb`, `ub`, `scipy_x`; minimize f'x,
with eta the last variable and coefficient -1. The single-positive
and six-active models have 145 variables; one-column has 146.
Root independently re-solves these three numerical models in MATLAB.
The exact six-active model and augmented model JSON are retained
because the independent audit reads them. Repeated 4 MiB diagnostic
matrices and cached stdout/stderr are excluded from the final
whitelist. Two changed entrypoint versions are frozen from the run
helper's actual `input-script.py` snapshots, rather than reconstructed.

All eleven recorded Python child processes exited 0 and were verified
absent. Actual versions, commands, hashes, PIDs and timings are in
`process-receipts.json`. This component did not commit or push.
The component manifest is a precise file whitelist for root integration.

The useful result is an actual shared outside input and an exact
correction of a finite-grid negative certificate. No general shared
obstruction, 3.5 deletion bound, extremizing sequence, or exact C2
was proved. General C2 bounds remain unchanged.
