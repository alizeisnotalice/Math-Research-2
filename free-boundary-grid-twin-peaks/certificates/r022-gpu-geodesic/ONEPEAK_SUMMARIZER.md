# Strict one-peak geodesic summarizer

`summarize_onepeak_geodesic.py` is an offline aggregation step; it does not
start a GPU job or modify the frozen Swift source.  Invoke it with a complete
JSONL, an explicit sidecar, and the actual source file:

```sh
python3 summarize_onepeak_geodesic.py \
  --input RUN.jsonl --sidecar RUN.sidecar.json \
  --source geodesic_prefix_r022.swift --out RUN.summary.json
```

The sidecar has exactly three objects: SHA-256 of `--source`, every run
metadata parameter to be fixed, and the explicit deterministic seed formula.
`onepeak_summary_sidecar_example.json` is an actual n=8 P=257 example.  The
program rejects a SHA mismatch, any metadata mismatch, anything other than one
metadata line followed by exactly R uniquely indexed replicates and one final
`run_end`, a missing replicate, an incorrect seed, or a tau other than N+1.

For every non-dead replicate it checks all group dimensions, the actual final
group capacity (including the final partial group), per-group and total
monotonicity in d, d=1=d=2=P, and `scaled=sum(groups)/P`.  It also cross-checks
the established inverse group fields against their old mean/scaled/Zmean
fields.  A non-dead record may have `Z=0` only as a floating underflow: its
finite `logZ` remains mandatory and is used in aggregation.

All cross-replicate weights are computed at Decimal precision 80 as
`exp(logZ-max(logZ))`.  For R complete replicates with normalized weights
`alpha`, the independent-replicate diagnostic uses
`Var(p)=R/(R-1) sum [alpha*(p-pbar)]^2`; for
`qbar=p[d+1]/p[d]` it uses the direct residual formula
`R/(R-1) sum [alpha*(p_next-qbar*p)]^2/pbar^2`.  Dead onepeak records are
rejected: a connected-grid one-birth prefix always has an available boundary.
For each d the output records nonzero replicate count, event-weight ESS, and largest event-weight share.  If no non-dead replicate
has the event, `observed=false` and estimate/ESS/SE are null: the output does
not call that a true probability zero.  `qbar=p[d+1]/p[d]` and its delta-method
SE are provided only where the event estimates exist.  SEs assume independent
replicates and are explicitly diagnostic, not a mixing certificate.

Run the negative and underflow tests with:

```sh
python3 tests/test_summarize_onepeak_geodesic.py
```

Those tests cover a missing final record, a missing replicate, over-capacity
group count, wrong seed, and a non-dead `Z=0` finite-logZ underflow replicate.
