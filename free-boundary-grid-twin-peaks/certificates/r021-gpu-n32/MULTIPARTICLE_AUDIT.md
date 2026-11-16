# P=17 full-trace audit

The frozen n<=32 implementation has source SHA-256
`fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629`
and binary SHA-256
`cd208dca6b7424d3835018e258ff0ef654516b99b000a6484342b953cb21dcde`.

The following four complete stratified traces were actually run with P=17:

| trace | L | logZ | dead | result |
| --- | ---: | ---: | --- | --- |
| n=24, tau=2 | 4 | -666.9867092937521 | false | full root trace audit; inactive inverse fields/groups |
| n=24, tau=577 | 4 | -749.008140223285 | false | full root and inverse-fibre trace audits |
| n=32, tau=2 | 4 | -1261.5869515001514 | false | full root trace audit; inactive inverse fields/groups |
| n=32, tau=1025 | 4 | -1571.0044377428585 | false | full root and inverse-fibre trace audits |

Every trace has one row for every rank, real ordered prefixes, and all 32
occupancy/boundary words.  The independent graph reconstruction observed
nonzero words beyond word 8 for all four cases.  Thus these are not n<=16
states padded to 1024 bits.

Both P=17 one-peak traces have a finite `logZ` but `Z=0.0`, because the
representable Double underflows; `dead=false` remains correct.  In the n=32
P=1 one-peak trace, the same condition occurs with nonzero inverse statistics:
`logZ=-1624.3282781663117`, `Z=Zmean_m2=Zmean_m3=0.0`, while
`mean_m2=scaled_mean_m2=1.0` and
`mean_m3=scaled_mean_m3=1018.0`.  The inverse verifier checks exactly this
log-domain representation rule, so an external aggregator must recover the
weight from `logZ` rather than treating a non-dead zero `Z` as a dead run.

`run_negative_trace_tests.py` altered disposable copies of that n=32 trace.
The full root verifier rejected a changed high occupancy word and a truncated
prefix; the inverse verifier rejected a changed terminal m3 group.  The
actual return codes are recorded in `negative-tests.json`.

The P=17 traces are verification artifacts only.  They are not probability
estimates and do not validate a high-particle n=32 production batch.
