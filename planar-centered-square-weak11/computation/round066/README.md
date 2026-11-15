# R66 fixed-theta competition-budget diagnostics

Prepared for actual execution by root. These new modules do not modify frozen
R63/R64 code and do not launch workers themselves. No execution is claimed by
this README before root supplies receipts.

For a prescribed proper maximum-mass strictly eligible block S, let L=m(S),
theta=3/4 and C=bbox(S)+[-sqrt(theta L)/2,sqrt(theta L)/2]^2. High competitors
are all strictly eligible bbox closures T that intersect S and have
m(T)>theta L. **S is retained**: a singleton's E_S can extend outside C.

`theta_competitor_budget.py` records exact rational inner/outer area bounds for

* C itself;
* low old credit |C intersect H(remaining)|;
* high unpaid area |(union high E_T) minus (C union H(remaining))|;
* the hole |C minus H(full)|;
* actual deletion Delta=|H(full) minus H(remaining)|;
* the analytic budget B=area(C)-low old credit+high unpaid.

All geometry aggregation uses rectangle unions. Individual competitor areas
and tiny-width scores are diagnostics only and are not summed as a union.
The exact identity B=Delta+hole is a mathematical interface from the core;
the code additionally checks interval compatibility of its separately enclosed
expressions. Since C is nested, B is nondecreasing in theta. The present theta
does not start an optimization scan. An input-specific successful budget is
not a general selection or compensation theorem.

The inner high-unpaid calculation uses high inner rectangles minus the union
of C outer and remaining outer rectangles. The upper calculation reverses all
enclosure directions. Credit uses C inner intersect remaining inner for its
lower bound and both outers for its upper bound. All radius enclosures are
certified by exact integer squared inequalities, without floating geometry.

For theta=3/4, the module also records a rational enclosure of 7/4-sqrt(3)
and tests the sufficient universal compensation inequality

    (high unpaid - low old credit)/L <= 7/4-sqrt(3).

It separately tests the more precise input-specific budget B/L<=7/2.

`erosion_ring_diagnostic.py` is a separate source version. For each competitor
T it records a=m(T intersect S), b=m(T minus S), and the ring between the
same full-T bbox receiver rectangles at masses m(T) and b. The latter rectangle
is contained in H(remaining). Ring area is at most a: with widths dx,dy fixed,
F(t)=(sqrt(t)-dx)_+(sqrt(t)-dy)_+ has derivative between 0 and 1 wherever
nonzero, so F(m)-F(b)<=m-b=a. Outward geometric upper intervals can exceed a
by endpoint rounding; the exact analytic cap is recorded separately.

Every ring is decomposed into at most four rectangles. The code computes the
union of all rings outside C, which bounds the actual high unpaid union.
Overlapping a's are never charged as if their sum were L. The standalone
module consumes a source-hash-bound PASS budget result and preserves that
dependency hash in its own result.

Root should first execute the R64 final 21-source input, then twelve independent
R63 strata if the pilot is meaningful. The R63 adapter preserves all 199 proper
maximum-block ties from the fixed 72-input family; the twelve improper full
maximum blocks are explicitly skipped. Limited pilots (`--case`, `--max-blocks`)
mark whether all proper ties were evaluated. Each completed block is atomically
checkpointed. Default budgets are 120 seconds per block and 900 seconds per
job. The geometry module is pinned to R64 SHA
`4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246`.
