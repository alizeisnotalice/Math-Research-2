# n24 local GPU batch recovery, r025

The r024 launcher created a detached driver (PID 13117), `/usr/bin/time`
wrapper (PID 13119), and GPU process (PID 13121).  The compute process was
not interrupted:

* `raw/n24-p131072-r64-l16.jsonl` has 66 lines: one metadata line, 64
  unique replicate indices `0,...,63`, and one final
  `{"kind":"run_end","wall_seconds":2930.665172100067}` line.
* Every replicate has `dead=false`; deterministic seed values are unique.
* The time wrapper wrote exit code `0`, 2930.73 real seconds, zero swaps and
  no signals.  Its maximum resident set size was 1,225,441,280 bytes and its
  peak memory footprint was 1,526,679,136 bytes.
* The telemetry final line is an `exit` record for PID 13119 with
  `returncode: 0`.  A post-run `ps` query found none of PIDs 13117, 13119 or
  13121.

The r024 driver did fail *after* the GPU computation and before it wrote its
sidecar, reports and receipt.  Its `validate_raw()` method required the
printed floating field `F_minus_2G_minus_4D` to be exactly zero.  Fifteen
records instead contain residuals of magnitude roughly `1e-216` to `1e-218`;
their integer per-group residual lists are all exactly zero.  Thus this is
ordinary final floating arithmetic roundoff in a printed summary, not a
failure of the integer identity or of the compute kernel.

The r025 package preserves the raw bytes and runs the frozen aggregator
(`cb466dca...`), which explicitly validates the same printed identity with a
scale-aware tolerance while independently checking every integer group
identity.  It accepted all 64 records and produced the three threshold
reports recorded in `receipts/n24-p131072-r64-l16.completion.json`.

No raw output was overwritten and no GPU computation was restarted.
