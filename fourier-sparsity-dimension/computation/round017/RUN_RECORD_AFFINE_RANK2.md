# R17 finite characteristic-two affine-rank classification audit

## Scope

This is a complete, exact enumeration only of alternating 4-by-4 matrices over
F2 and their affine subspaces.  It is finite evidence for the candidate proof
route in `pending-proof-route.md`, not a proof of its arbitrary-dimensional
classification claim.

Forms use the six coordinate order `(01,02,03,12,13,23)`.  The script first
generates every subspace of F2^6 by closure, verifying the Gaussian-binomial
counts `(1,63,651,1395,651,63,1)`.  It then enumerates every coset of every
subspace, does GF(2) Gaussian elimination on each 4-by-4 alternating matrix,
and accepts precisely cosets whose forms all have rank two.  For every
accepted coset it tests common support intersection (the star condition) and
the dimension of the sum of support planes (the fixed-core condition).

## Exact result

The full affine-subspace coverage counts by affine dimension were
`(64,2016,10416,11160,2604,126,1)`, totaling 26,387.  The accepted
constant-rank-two counts were:

| affine dimension | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| constant rank two | 35 | 595 | 840 | 105 | 0 | 0 | 0 |
| with rank-four difference | 0 | 280 | 630 | 105 | 0 | 0 | 0 |

Thus the largest observed affine dimension with a rank-four difference is
three.  Of the 560 families whose difference subspace has rank at most two,
455 have a common star factor, 455 have summed support dimension at most three,
350 have both properties, and zero have neither.  The underlying full form
rank histogram is rank 0: 1, rank 2: 35, rank 4: 28.

Two exact higher-dimensional constructions are included in the result:

- In F2^8, `e0^e7 + span(e0^e1,...,e0^e6)` is an affine dimension-six star
  family of 64 rank-two forms.
- The line `{e0^e1,e2^e3}`, embedded in four coordinates of F2^8, consists of
  rank-two forms and has rank-four difference `e0^e1+e2^e3`.

## Run evidence

- Remote host/work directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r017/affine-rank2-cpu`.
- Source SHA-256:
  `416d7f9fa61dd49ccd4ced57aaec5919473138d756b26c0b1fa56a110f3eac07`.
  Result SHA-256:
  `65d2c5312537422b351ac20cfb1aa67aceb706a77cec767091eb3961a5b36c6e`.
- Command: `/usr/bin/python3 ./enumerate_affine_rank2_4x4.py`; Python 3.9.6.
  Final driver exit code was zero, and the enumeration reported 0.293137 s.
  A post-run process query found no owned process.
- An initial implementation exited one before producing output because its
  assertion incorrectly imposed the <=3 fixed-core predicate on rank-four
  difference families.  The correct finite test applies the candidate's
  separate `dimension <=3` rule to that branch, and applies star-or-core<=3
  only to rank<=2-difference families.  The repaired source was rehashed and
  rerun once to the recorded final result.
- Host samples adjacent to the first short run reported 12.68% and 13.8% idle;
  a later postcheck reported 4.95% idle, but these are aggregate host readings,
  not attributable to the 0.293-second single-process enumeration.  No claim
  of sustained A28 CPU saturation is made.

The compact JSON output is
`/Volumes/data/CodexTemp/T2/a28/r017-enumeration-result.json`.
