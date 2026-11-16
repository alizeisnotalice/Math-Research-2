# R73 four-control regional credit diagnostic

Actual date: 2026-11-16, Asia/Shanghai. Python3.13.18, standard-library
Fraction/integer arithmetic. This component reuses the **already certified**
72-bit full/residual halo rectangles of four R72 controls. It does not
re-enumerate old source closures, recompute square roots, or run an old
complete-halo evaluation. New calculations are clipped-region sweeps and
exact unit-window source-mass enumeration.

## Regions and quantities

All four controls have S=S5, mass1 and bbox[0,24/25]^2. Fix

    C=[-1/2,73/50]^2,       |C|=2401/625,
    J=[-1/25,1]^2,         |J|=676/625.

Each axis is partitioned into the three genuinely disjoint intervals
`[Clo,Jlo)`, `[Jlo,Jhi]`, `(Jhi,Chi]`. Their nine Cartesian products
partition closed C exactly. The middle product is closed J, four
side products have one middle coordinate, and four corner products
have no middle coordinate. The area computation may clip at closed
bounds, since the boundary differences have Lebesgue measure zero.

For each piece K, the report encloses

    full = |K intersect H(sigma+nu)|,
    residual = |K intersect H(nu)|,
    deletion = |K intersect (H(sigma+nu) minus H(nu))|,
    missing = |K minus deletion|,
    vacancy = |K minus H(sigma+nu)|.

Actual residual halo is contained in actual full halo; therefore
`missing=vacancy+residual`. Every reported area is a directed rational
interval obtained from saved integer inner/outer rectangles. Decimal
values below are only displays. The new vertical x-event interval-merging
sweeps use direct interval merging rather than the old vertical segment tree.
The initial v1 doc called this horizontal; that orientation label was wrong.
The mathematics and all region values are unchanged. Method correction and
a v2 metadata-only replay preserve the original bundle and its hashes.

Full inner minus residual outer gives the deletion lower bound; full
outer minus residual inner gives its upper bound. We do not subtract
two approximate displays. Exact clipped C totals are checked against
the inherited complete deletion intervals; clipping only tightens
possible outer-error bands outside the exact deletion container C.

## Four controls

| input | deletion in J | deletion in side strips | deletion in four corners | total deletion in C |
|---|---:|---:|---:|---:|
| S25 only | 1.0816 | .2779936521 | .0144 (exact9/625) | 1.3739936521 |
| R71 positive33 | 1.0815020306 | .3327177757 | .0277612836 | 1.4419810900 |
| R72 aligned38 | 1.0816 | .4221583272 | .0514298908 | 1.5551882181 |
| R72 contracted unique38 | 1.0816 | .4221583268 | .0514298903 | 1.5551882171 |

The corner-only sufficient route proposed by root is

    corner missing >= |C|-7/2 = 427/1250,

or equivalently

    corner deletion <= 631/1250 = .5048,

since total corner area is529/625. It passes with substantial **exact
rational** margin in all four controls. This is not a general theorem.

For the contracted input, corner missing is approximately.7949701097:
full-halo vacancy contributes.3426814214 and residual halo.4522886883.
The credit above the proposed required gap is about.4533701097.
The original aligned input's **full** corner halo is approximately
.5037185822, close to.5048, but its actual corner deletion is only
.0514298908 after residual subtraction. The diagnostic consequently
supports keeping actual residual credit in any proposed payment.

The contracted input has C-minus-J deletion approximately.4735882171,
vacancy1.2907276235, residual.9956841593, and missing2.2864117829.
The R71 control alone has a small residual incursion into J,
approximately.0000979694; the other three controls have zero such
incursion within the certified endpoints.

## Exact unit-window caps

For a closed unit square with origin (x,y), all finite capture changes
occur at x or y equal to a source coordinate or that coordinate minus1.
Every such point and every open interval between consecutive events
is represented. Exterior origin intervals capture no source in that
axis and cannot add a nonempty capture pattern. Thus these are complete
finite exact unit-window mass diagnostics, including boundary ties.

The report groups windows by captured internal k by l grid shape and
by left/right/bottom/top edge or corner families. Each group saves the
maximum outside mass, its exact capacity `1-kl/25`, the minimum slack,
a representative closed unit square, and its full source capture mask.
All exact masses satisfy the unit cap. The numbers of distinct
saturating full capture masks are1,2,4,1 for the four controls, matching
their inherited maximum-block counts.

The aligned input has saturated right-edge windows with k=3,l=5,
outside mass2/5, and right-corner windows with k=3,l=4, outside mass13/25.
After the exact outside contraction, their slacks become respectively
`1/250000000` and `13/2500000000`. Only the all-S capture mask remains
saturated. The full per-k,l and family tables remain in the JSON;
they are not compressed to a claim that all edge capacities saturate.

## Reproduction and boundary of evidence

From the topic root, using the existing Python environment and data
temporary/cache directories:

    python computation/round073/regional_credit_diagnostic_v2.py --archive /absolute/path/to/integrated/topic --out computation/results/round073-regional/four-control-regions-v2.json

The archive must contain the four exact inputs under inputs/r072 and
their complete geometry under results/round072-actual. Every source is
bound by SHA in the output and component manifest. The first diagnostic process
took.233063 seconds internally. A v2 orientation-description replay used the
same inputs and every semantic output field is exactly equal. Both processes
exited0 and were checked absent.
The saved execution receipt records actual versions, PID and command.

No new source input, MIP, general corner bound, 3.5 budget, or C2 bound
is claimed. These four fixed controls guide the next edge/corner
inequality; they do not validate it for arbitrary actual outside measures.
Root is independently running additional structured controls. This
agent is waiting for a specific payment inequality before further probes.
