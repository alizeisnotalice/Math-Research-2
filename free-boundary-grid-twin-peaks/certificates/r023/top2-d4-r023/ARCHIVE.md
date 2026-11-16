# r023 exact birth-rank-two D4 accounting

`top2_d4.cpp` implements the D4-orbit independent-set integer calculation of
the birth-rank-two (`tau=2`) top-two-root counts, Euclidean moment and
Manhattan histogram for n=3,…,7. `summary-n3-n7.json` is the authoritative
normalized view; it uses `Q=D/(n^2 B)`. The raw per-n JSON retains an old
*derived display* error: its `tau2_metrics.Q_den` and `K2_compare.Q_den`
omit n². Its integer counts `Btop2,D,F,G` and histograms are unaffected.

The independent full-subset boundary-DP implementation and its n=3,…,5
outputs live in `../fixed-state-completion-root-r023/top2-reference/`.
`top2-independent-receipt.json` records equality of every integer field and
histogram. All numbers here are finite exact results; the tau=2 law is not the
uniform law conditional on exactly two peaks.
