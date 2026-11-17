# Packaging replay receipt

Actual replay in an isolated scratch directory used CPython 3.15.0.

* The generic fibre, seven-vertex exhaustive, selected nine-vertex, and
  counterexample scripts in `analytic-r45/` were executed. Their saved JSON
  outputs matched the replayed outputs; the seven-point exhaustive output has
  no `status` field but records zero counterexamples over 5184 hosts.
* `independent-prefix-dp/verify_r45_prefix_dp.py` was rerun and passed.
* Both path-parameterized n=6 and n=7-centre CRT auditors passed. They use
  integer arithmetic, recursive Lucas/Pratt validation, direct CRT, and
  Garner CRT; no production enumerator was rerun.
* `short-gap-consequences/audit_short_gap_consequences.py` was then rerun
  against those two newly produced public auditor JSON files. It reports
  lifted `L>=6` fractions 0.3256225227465938 (n=6 global) and
  0.639629933728095 (n=7 centre only); all n=6 t=4,5,6 QE8 slacks are
  nonnegative and all centre-only slacks are negative.
* The initial packaging interpreter was CPython 3.15.0 without SymPy, so it
  could not run the symbolic phase. The coordinating agent subsequently ran
  the unchanged script successfully with Python 3.13.18 and SymPy 1.14.0:
  exit 0, elapsed 0.999 s. The path-normalized receipt and stdout are in
  `spatial-r45/replay-receipt.json` and `spatial-r45/replay-stdout.txt`.

The n=7 centre run requested 12 threads with one process at a time. The four canonical time receipts total
872.64 CPU seconds (user+system) and 269.96 wall seconds, or about 3.23248
core-equivalents averaged across these runs. This is process CPU time, not
whole-host utilization and not evidence of full twelve-core use.
