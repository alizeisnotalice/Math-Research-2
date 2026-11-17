# r034 N=9 bipartite-prefix certificate

This package preserves the exact finite certificate for
\(R_{27}=h_2-5h_4+5h_6-2h_7\) on the four fixed nine-vertex bipartite
hosts \(1+8,2+7,3+6,4+5\).  The shipped strict ledger has complete,
non-overlapping coverage of

\[
2^8+2^{14}+2^{18}+2^{20}=1,327,360
\]

edge subsets: 326 completed checkpoints split as 1, 4, 64, and 257 by left
host size.  Each reports exact integer `max_abs_R27 = 0`.

`proof/NINE_PREFIX_LIFTING.md` is a separate analytic argument: conditional
on this finite local statement, its prefix--completion decomposition lifts it
to every connected bipartite graph with \(N\ge9\).  The finite machine run
does not itself establish an asymptotic or peak-distance result.

The main finite scan is one subset-DP program over all four N=9 host families.
The direct `9!` enumeration was used only on three gate graphs (`P9`,
`K23-tail`, `C6-tail`) to compare with that DP.  It is not a second independent
full-N=9 scanner.

## Contents

- `source/`: portable source. `bipartite_prefix_dp.cpp` is byte-identical to
  the certified source; `local_cpu_run.py` has only its private cache location
  replaced by the required `CLANG_MODULE_CACHE_PATH` environment parameter.
- `evidence/results/`: all 326 completed raw checkpoints, unchanged.
- `evidence/strict-aggregate.public.json`: the strict interval aggregation,
  with only its local output path redacted.
- `evidence/execution-exits.public.json`: one path-redacted exit record for
  every counted completed checkpoint (325 production plus the counted
  benchmark prefix).
- `evidence/gates/`: raw gate results. `evidence/run-receipts.public.json`
  records normalized compile, gate, aggregation, completion, and resource
  receipts without a binary or log flood.
- `scripts/audit_n9_certificate.py`: recomputes all checkpoint hashes and
  coverage and checks the source, gates, and exit ledger without rescanning.

The certificate is frozen: use a pre-created output directory outside this
tree for any optional auditor JSON output, as shown in `REPLAY.md`.

`PROVENANCE.md` gives original hashes and redaction details.  No executable,
large trace, compiler cache, or private absolute path is shipped.
