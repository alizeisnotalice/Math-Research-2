# R74 one preliminary Borel-corner clipping diagnostic

Actual date: 2026-11-16, Asia/Shanghai. This is only the one minimal
pre-diagnostic requested before the formal R74 topic is fixed. The R73
frozen v2 bundle is unchanged. No classifications, subdivisions, MIP,
source search or MATLAB process was run.

For each inherited R72 Borel family member n=5,17, put

    e=1/n^2, w=1-e,
    C=[-1/2,w+1/2]^2, J=[-e,1]^2.

The four corners are the products of the left and right components
of C-axis minus J-axis. Each corner has side1/2-e. The candidate
corner-only deletion threshold is

    total corner area - (|C|-7/2) = 1/2+3e^2.

The Borel outside measure has density t=n/(n+1)<1, so its residual
strict halo is empty by the inherited proof. Hence full-halo area in
these corners equals actual deletion area for this Borel family.
This interpretation does not apply to finite atomic approximations.

The saved adaptive certificates partition C into classified and
unresolved leaves. For each leaf, this diagnostic computes its exact
integer rectangle intersection with each corner. State1 intersection
area enters both lower and upper sums; state0 enters the upper sum
only; state-1 enters neither. It verifies old global state counts and
area sums, plus that the three clipped state areas sum to each exact
corner area. It does not replay the old leaf classifiers or partition
proof. Source summaries, gzip certificates, and independent-review
records are all SHA-pinned.

| n | corner | lower | upper |
|---|---|---:|---:|
| 5 | left-bottom | .04161516552 | .04218033266 |
| 5 | left-top | .04132683182 | .04246678314 |
| 5 | right-bottom | .04132683182 | .04246678314 |
| 5 | right-top | .04132683182 | .04246678314 |
| 17 | left-bottom | .03252975811 | .03476512448 |
| 17 | left-top | .03252975811 | .03476512448 |
| 17 | right-bottom | .03252975811 | .03476512448 |
| 17 | right-top | .03143104293 | .03606400854 |

The unequal enclosure widths reflect the saved adaptive refinement,
not a claim that the symmetric actual corner areas differ. No symmetry
tightening has been applied to these raw clipped sums.

The exact total-corner bounds are

    n5: [868198179/5242880000, 1778182333/10485760000]
        approximately [.165595661,.169580682],
        threshold 631/1250 = .5048;

    n17:[664668019/5151653888, 12292410251/87578116096]
        approximately [.129020317,.140359382],
        threshold 83527/167042 approximately .500035919.

Both upper bounds are strictly below their thresholds by exact
rational comparison. These two Borel controls therefore do **not**
refute the corner-only candidate. They do not prove that candidate
for arbitrary Borel or finite atomic outside measures.

The n5 lower bound does refute the stronger Borel-only proposal that
each corner receives at most one quarter of its internal corner mass.
Its total proposed budget is4/25=.16, while the saved lower bound exceeds
it by the exact amount29337379/5242880000. Each individual corner lower
bound also exceeds1/25. A separate endpoint-comparison script verifies
these inequalities without rereading or reclassifying any leaf. There
is no finite-atomic extrapolation.

Inherited review scope remains unchanged: for n5, the partition and
all classified leaf predicates were independently replayed; for n17,
the whole partition and area ledger were checked, but only16 classified
leaves of each sign were independently replayed. The rest of n17's
classification remains inherited production evidence under its reviewed
algorithm, not a new full independent replay.

The new script uses Python3.13.18 and standard-library integer/Fraction
arithmetic. It reads 198139+306433 existing leaves, took1.137696 seconds
internally, exited0 and its PID51172 was verified absent. The receipt
and source pins record the actual execution. Reproduce from topic root:

    python computation/round074/clip_saved_borel_corners.py --archive /absolute/path/to/integrated/topic --out computation/results/round074-borel-corners/saved-n5-n17-corners.json

Only r074 files were written. No commit or push was performed.
