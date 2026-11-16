# Completed n=16 million-particle batch

This is a minimal, self-contained archive of the **completed** n=16,
P=1,048,576, 16-replicate, tau=2,...,256 stratified SMC batch.  It is a
finite numerical experiment, not a proof about the limiting conditional law.
It deliberately does not include the earlier n=16 pilot data, an executable
binary, rank traces, or a claim of exclusive GPU use.

## Contents and identities

- `bitset_smc_r017_checked.swift` is the exact production host source used by
  this batch: SHA-256 `65aefc323b20995564df7a65656b4bfffbc1570487f8be132fafa916225c2f67`.
  The embedded Metal bytes (from immediately after `let metal = """` through
  the byte immediately before the closing delimiter) hash to
  `6c59a560e34e746633fa5e6c5ba886c759786d9da2df4d313af9a32258a32181`.
- `raw-full-n16-p1048576-r16-d7500.jsonl` is the complete raw run:
  SHA-256 `49ef2fc3ba1b63cf2eacb07ee267685554e37eeac595ded1b0d350bd0af5024a`.
- `telemetry-full-n16-p1048576-r16-d7500.jsonl` is the complete target-RSS and
  global-memory-pressure monitor record (not a GPU-utilisation log): SHA-256
  `9a645e3b84a76463ae5160af17a48abf1dc7beb6a34620078c4f9b2f36fd35c9`.
- `full-n16-p1048576-r16-d7500-by-replicate.json` is the original Decimal-80
  aggregate; `summary-recomputed.json` is a fresh local recomputation from
  this archive's raw JSONL.  They agree after their deliberately different
  `input` path values are ignored.

The accepted run command was `n=16 P=1048576 R=16 lo=2 hi=256`, domain
`7500000000000000`, method `strat`.  Its raw rectangle has 4,080 replicate
records, one metadata record and one final `run_end`.  The deterministic
seed rule is checked in `verify_completed_million.py`.

## Reproduce the archival checks

Run from this directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 summarize_full_by_rep.py \
  raw-full-n16-p1048576-r16-d7500.jsonl > summary-recomputed.json
PYTHONDONTWRITEBYTECODE=1 python3 verify_completed_million.py \
  > completed-million-verified.json
```

These commands are light checks; they neither compile nor run the Metal
program.  The verifier checks the three immutable-file hashes, exact metadata,
the complete tau-by-replicate rectangle, every source-specified seed, the 32
dead rows (all 16 repeats at each of tau=255 and tau=256), final `run_end`,
and the monitor's natural `returncode: 0` event for PID 6859.

`root-audit/audit_million.py raw-full-n16-p1048576-r16-d7500.jsonl` is an
additional direct Decimal audit of the same raw file (source SHA-256
`a67c8cb522c213d33024526aa6be0b043b4471d660a1fae9df9ee5c8d782f28a`).
Its replayed `root-summary.json` supplies the fixed-strata and paired-repeat
delta diagnostics, repeat-weight ESS, and birth-rank tail.  It is descriptive
finite-run evidence only, not a confidence or convergence certificate.

The raw final wall time is 15,432.929917931557 seconds.  `monitor-summary.json`
records its 15,105 one-second samples: target RSS 5,968--1,094,944 KiB and
global free memory 23--74%.  The monitor covers only the target process RSS
and system memory pressure.  Device access was
shared with other work, so neither this telemetry nor this archive certifies
exclusive GPU use or GPU-wide utilisation.  A fresh process listing after
archiving found no `bitset_smc` or `run_with_memory_monitor` process from this
batch; the historical monitor provides the stronger PID-specific exit receipt.

## Numerical summary (finite experiment only)

The fresh Decimal-80 aggregate reports

- mean of per-replicate `sum_tau Z`:
  `2.0815567336456518498295129676100556235880143905069916003881023063207074472283856E-79`;
- corresponding Z-weighted `U_over_Z`:
  `0.30148899518092508752329714191200089962543691139461263285025095449327875422778032`;
- repeat log-normalizer ESS: `10.716576158232591360618977575981272534260508326787927642716146543938793524216484` out of 16.

These are output summaries for this particular SMC implementation and finite
run.  They are not an unbiasedness, convergence, mixing, or asymptotic
certificate.
