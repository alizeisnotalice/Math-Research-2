# R79 quadratic continuous-old envelope: interface and implementation

This interface was written before production permission. One authorized
production batch has now completed; see `PRODUCTION.md` for actual runtime,
results and the pending independent-audit status. The R78 components and
their results remain frozen and unchanged.

The receiver coordinates are reflected to X=-x,Y=-y, with domain
[1/400,1/2]^2. The same reflection is applied to outside atom positions.
For n=20, h=21/400, c=400/441, old capture mass is bounded above by
c(r-X+h)(r-Y+h) when r>max(X,Y). For outside subset mass B, define

```
F(X,Y,r)=(4-c)r^2+c(X+Y-2h)r-c(X-h)(Y-h).
```

The continuous quadratic set consists of points admitting some strict
radius R(X,Y)<r<1/2 with B>F(X,Y,r), where R includes the origin and every
outside source in the subset. This set is an **outer state** for actual
old-involving halo witnesses. Being inside this state does not prove actual
old capture or actual halo membership.

Two variants are evaluated on the same leaves and actual residual:

1. `exact_pureold_plus_three_mixed`: the true pure-old Gamma halo, the origin
   rectangle (1/400,1/40)^2, plus the three nonempty-outside quadratic sets.
2. `four_soft_branches`: all four quadratic subsets, including B=0.

The original 402-atom input files are byte copies of the frozen R78 inputs.
Their inherited scope metadata refers to those earlier computations; the
new calculation has the scope stated here. No parameter search is encoded.

## Exact cell classification

Let C=[u,v]x[s,t] and let [a,b]x[d,e] be the bbox of the origin and selected
outside atoms in positive coordinates. Exact radius extrema are

```
Rlo=max((b-a)/2,(e-d)/2,b-v,u-a,e-t,s-d,0)
Rhi=max(b-u,v-a,e-s,t-d,0).
```

The quadratic minimum over [R,1/2] is obtained by rationally clamping
r0=c(2h-X-Y)/(2(4-c)). The coefficients, minimizer and minimum are all exact
Fractions; no root or floating optimizer is used for this classification.

- `out`: Rlo>=1/2, or min F(u,s,r) on [Rlo,1/2] is at least B.
- `in`: Rhi<1/2 and min F(v,t,r) on [Rhi,1/2] is strictly less than B.
- otherwise `unknown`.

These statements concern the continuous outer state. For `in`, the strict
margin allows perturbing a minimizer at an endpoint into the open radius
interval. Equality never certifies `in`. Monotonicity in X,Y on the
admissible radius domain gives both safe directions. The separately owned
theory note is to be reviewed before a production run is authorized.

## Identical actual residual and area bounds

All three genuine residual labels are retained. Their positive-coordinate
center rectangles use exact integer-isqrt root intervals of width at most
1/(400*2^72); ineligible empty receivers are discarded by exact comparison.
This is the same directed construction as the pinned frozen R75/R78 code,
implemented locally without altering or executing those old scripts.

Rectangle area minus the residual union is evaluated exactly by
inclusion-exclusion over at most three rectangles. For every leaf:

- an envelope-inside cell contributes C minus residual-outer to the state
  lower bound and C minus residual-inner to the state upper bound;
- an envelope-outside cell contributes zero;
- an unknown cell contributes zero to the soft-state lower bound and
  C minus residual-inner to its upper bound.

For the exact-pureold variant, the exact intersection with the pure-old
rectangle is retained even when the mixed state is outside or unknown.
The whole-cell upper bound still applies when mixed membership is unknown.
Cells covered by the residual-inner union contribute zero and need no
polynomial tests. This does not assert that the envelope set itself is
empty in those cells.

Thus lower/upper brackets enclose |E_quad minus Hnu|. Only their **upper**
endpoints give an upper bound on actual Gamma deletion. The state lower
endpoint may diagnose a loose state but cannot be used as an actual
deletion lower bound.

## Complete coverage and certificate schema

The initial domain is divided into 2x2 disjoint-interior tiles per input.
Each tile is processed independently, with no nested pool or BLAS threads.
At most eight such tasks run concurrently for the proposed minimal run;
the CLI permits no more than twelve workers.

Each split replaces one closed cell by all four midpoint children, ordered
lowerX/lowerY, upperX/lowerY, lowerX/upperY, upperX/upperY. Boundary overlap
has area zero. A leaf's tile id and base-four path uniquely reconstruct its
rational box. Every leaf remains in the compressed certificate, including
unresolved leaves if a node/depth cap is reached. Full coverage can be
replayed from this complete prefix code; no omitted branch is inferred
from numerical success.

Per-leaf records contain the four branch states, mixed/soft union states,
exact residual-complement areas, exact pure-old remaining areas, both
variant lower/upper contributions, and raw unresolved-cell areas. All
mathematical numbers are rational strings. Every tile also records
created/split/leaf counts, maximum depth, termination reason and totals.
The node identity leaf_count=1+3*split_count and exact domain area sum are
checked before serialization. The final per-case certificate contains all
tiles, source/input/parameter hashes, actual residual root data and the
precise child order.

The adaptive priority is the largest per-leaf area uncertainty across the
two variants. Each tile target is the case target divided by four, so
meeting all tile targets establishes case area gaps <=1/10000 for both
variants. The proposed first-run limits are 20,000 created nodes per tile
and depth 24. A resource limit returns complete remaining leaves and a
truthful non-success status; it never converts an unresolved branch into
an outside cell. Unresolved leaf area and the tighter residual-subtracted
area gap are reported separately.

## Production command

From the topic directory, using the existing Python interpreter:

```
TMPDIR=/Volumes/data/CodexTemp/T2 TMP=/Volumes/data/CodexTemp/T2 TEMP=/Volumes/data/CodexTemp/T2 XDG_CACHE_HOME=/Volumes/data/CodexTemp/T2/cache /Users/zzh/.venvs/main/bin/python3 /Volumes/data/CodexData/skills/python-math/scripts/run_math.py --python /Users/zzh/.venvs/main/bin/python3 --workdir /Volumes/data/CodexProjects/T2/.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11 --log-root /Volumes/data/CodexTemp/T2/r079-quad --timeout 600 computation/round079-quad/quadratic_envelope.py --parameters computation/inputs/r079-quad/parameters.json --output-dir computation/results/round079-quad/production --workers 8
```

The production entry point refuses to replace an existing summary or case
certificate. Runtime versions, actual parent/worker PIDs, timing and hashes
are recorded in the production summary and receipt. Observed process exit
checks and actual-R78 comparisons are in the separate postprocess result.
Syntax-only checks during implementation are not mathematical evidence.

Existing dependencies to retain in the eventual component ledger:

- corrected actual input SHA256
  `d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a`;
- far-face actual input SHA256
  `e343e448b9582a731ddaccca37452adcce1c6f28b4153aa520a0ec2b4a3ef574`;
- frozen R75 directed-radius/union source SHA256
  `a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`;
- frozen R77 unequal-axis caps theorem SHA256
  `3c71513ed72567445d4900b96f819f49d59481a7a6e9a1d1b4c1550cb4a211aa`.

No new agent, MATLAB, commit or push is authorized here.
