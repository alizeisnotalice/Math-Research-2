# r022 tau=2 lift calibration

This package uses the root-audited executable
`/Volumes/data/CodexTemp/T2/a022-tau2-lift-root-audit-r022/tau2_lift`
(binary SHA-256 `2a56637121f9a9d24a6f5771613ec00b176e719735e905ac20e1eabf779eae81`;
source SHA-256 `e268035742b1c8eb7b0e8c366e3b8e393054352ed0755cd523770e1d8d7d9c3a`).
The source was previously replay-audited at n=32; no source or binary was
modified here.  Each run used integer stratified resampling, tau=2 only,
16 MH moves, and a distinct domain.

| n | P | R | domain | elapsed wall seconds | strict aggregation |
|---:|---:|---:|---|---:|---|
| 3 | 4096 | 64 | `9c24010000000000` | 2.8931 | pass |
| 5 | 8192 | 64 | `9c24020000000000` | 9.5489 | pass |
| 7 | 16384 | 64 | `9c24030000000000` | 27.3592 | pass |

Each raw file has one metadata row, exactly 64 tau=2 replicates, and one
`run_end` row; each monitor stream ends with return code zero.  The strict
aggregator verifies sidecar source/raw hashes, expected metadata, seed formula,
no tau=2 dead trajectory, complete replicate set, group capacity and integer
range bounds, exact per-group lift identity, scaled and complete-lift fields,
and log-normalizer handling.  Its JSON outputs are in `reports/`.

The numerical comparison in `reports/exact-comparison.json` uses:

- n=3: `a022-second-birth-survival-r022/second-birth-survival-n3.json`;
- n=5 and n=7: `a022-root-moments-exact-r018/n5-root-moments.json` and
  `n7-d4-root-moments.json`.

It compares the lift full normalizer with exact `A2/(n^2)!` and the three
ratios `Q=E[d^2]/n^2`, `V=F/(8A2)`, `C=G/(4A2)`.  The absolute differences in
Q were (1.27090000514e-4) for n=3, (4.69e-4) for n=5,
and (9.47e-5) for n=7 Q.  The n=3 exact normalizer is `0.3962962963` versus
estimate `0.3961801683`; n=5 is `0.0004553649562` versus
`0.0004566374127`; n=7 is `2.2535264803e-10` versus `2.2602400963e-10`.

This is a finite simulation calibration only.  The reported delta standard
errors are independent-repeat diagnostics, not mixing guarantees, confidence
intervals, or a proof of the target asymptotic statement.
