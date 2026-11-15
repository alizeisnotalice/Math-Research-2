# R68 N128 repeated-capacity family

One actual input was evaluated: N=128, 257 positive sources. The producer
enumerated all 25,153 tight bbox closures, of which 15,205 are strictly eligible,
and computed exact directed union areas:

| Region | Upper display |
| --- | ---: |
| W | 1.0116576104946515 |
| Wprime = true deletion | 0.9846617144651398 |
| W overpayment | 0.02699589602951174 |
| Old halo H(S) | 0.9383075929213543 |
| Full halo | 1.23466171446514 |

L=1. These finite results are far below 3.5, and no larger N was run. N32 was
not rerun. Independently checked class carriers give the all-family bound
W < 399/128+7sqrt(2)/8-25sqrt(3)/32 < 389/128 < 3.5. The core's
[finite-carrier and exact-deletion proof](../../proofs/r068-two-row-family-finite-carriers-and-exact-deletion.md),
sections1–3, supplies the containment argument and mathematical scope behind
the producer's scalar comparison. The scalar assertion in code is not used
as a replacement for that proof. No new general C2 bound follows.

## Complete special source catalog

Old sources are p_i^+-=(7i/(8N),+-7/16), i=0,...,N-1, each of mass1/(2N).
The external source q=(39/32,0) has mass1/4; S consists of all old sources.

Without q, a tight closed bbox contains one contiguous interval [i,j] in the
lower row, upper row or both rows. With q and old sources, its horizontal left
edge is some x_k and its right edge is q_x; it contains a suffix [k,N-1] in
the lower, upper or both rows. The only remaining closure is q alone. Every
listed bbox is tight and closed, and these cases exhaust all coordinate bboxes.
Thus the exact count is 3N(N+1)/2+3N+1, with no 2^(2N+1) enumeration.

`catalog.json.gz` stores every descriptor and its eligibility test. Kinds are
0/1/2 for internal lower/upper/both-row intervals, 3/4/5 for the corresponding
external suffixes, and 6 for the external singleton. Source coordinates in the
catalog are integer multiples of1/(32N), masses are integer multiples of1/(4N).
The radius-floor column, present only on eligible rows, is in the scale
D=(32N)2^bits. The catalog is deterministic gzip with mtime0.

## Exact area method and state definition

Mass units u give the exact root-floor test

    k² <= 64 N u 2^(2bits) < (k+1)².

All non-singleton mass radii are enclosed by (k-1)/D and (k+1)/D. The residual
source radius1/4 is rational and used exactly. At N128 with bits64, the actual
coordinate mesh is2^-76. The run has139 distinct mass radii.

Internal rings have b=0. Every crossing label has b=1/4: W subtracts the
receiver at this mass on the full-T bbox; Wprime subtracts the actual external
square Q. Rectangle differences are decomposed and then aggregated as a true
union. The program uses an x sweep with an exact integer y-coverage segment
tree, not the earlier strip-merging routine. There are O(N²) rectangles/events
and O(N) distinct y levels, giving O(N² logN) union work and O(N²) memory.

The old halo lies strictly left of Q: its right edge is at most53/56, while
Q starts at31/32, a gap5/224. The producer additionally checks separation even
at its outer directed endpoints. Consequently Wprime equals true deletion for
this family. Its independently assembled directed union agrees at both levels
with full-halo area minus1/4. W's overlap with Q is separately measured and
strictly positive; it gives the actual difference between the two states here.

## Execution and reproduction

Python3.13.18, standard-library integer/Fraction arithmetic. The mathematical
child PID13783 exited0 and was confirmed absent. Wall1.35069 seconds; sampled
peak RSS52,256 KiB. Sampling is not a guarantee of the exact peak. The receipt
retains the absolute actual command, source hash and stdout/stderr.

From the repository root:

```sh
python3 -B planar-centered-square-weak11/computation/round068/run_n128_control.py \
  --producer planar-centered-square-weak11/computation/round068/two_row_erosion_family.py \
  --out-dir /Volumes/data/CodexTemp/T2/r068-replay/N128
```

This is self-contained: there are no runtime imports from old workspaces or
R64/R66 kernels. The minimal explicit257-source input is archived at
`computation/inputs/r068-two-row-N128.json`. Previous family sources, locations,
hashes and the distinction between inherited census and this new area run are
recorded in `family-source-index.json`.
