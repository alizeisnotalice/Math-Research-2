# Exact tau=2 top-two root D4 computation

`top2_d4.cpp` is the D4-orbit independent-set integer computation. It uses
12 worker threads and computes each state recurrence once. The pair weight is
calculated from \(B=V\setminus N[I]\) at that representative, rather than
rerunning a recurrence for each root pair.

This archive is deliberately flat: `n3.json` through `n7.json` and
`summary-n3-n7.json` are adjacent to this README. The independently written
full-subset boundary-DP source and its n=3,…,5 results are at
`../fixed-state-completion-root-r023/top2-reference/`.

From this directory, run the portable finite recheck:

```sh
python3 -B verify_top2_archived.py
```

It recompiles the independent subset-DP source in a temporary directory,
checks its n=3,…,5 output against the retained D4 integers and histograms,
checks \(F-2G=4D\), and recomputes the normalized summary. It does not rerun
the n=6/7 enumeration. To generate a separate summary without altering the
archive:

```sh
python3 -B summarize_exact.py n3.json n4.json n5.json n6.json n7.json \
  --out /Volumes/data/CodexTemp/T2/top2-d4-r023-summary.json
```

The raw n-specific JSON retains a historical *display* defect in
`tau2_metrics.Q_den` and `K2_compare.Q_den`: \(n^2\) was omitted from those
two derived denominators. `summary-n3-n7.json` is the authoritative view and
uses \(Q=D/(n^2B)\). The raw integer fields \(B,D,F,G\) and all histograms
are unaffected. `MATH.md` gives the normalization and distinguishes the
birth-rank-two law from the law conditional on exactly two peaks.
