# Provenance and I/O-only public adaptations

The completed source run is identified here by hashes rather than by a private
local path.  Its certified C++ source SHA-256 is
`be18e9f8c593648c4e2121305f400db5199ed2fcdde46c4522ddbae27b15b5c4`.
The historical compiled binary SHA-256 was
`a0a80708e1bc2cc732bc94e8f3bc824a2cf019adf187103b5227a1cc3e1cf18f`;
no binary is distributed, so public users cannot be expected to reproduce that
hash across compilers or platforms.

| Original item | Original SHA-256 | Public status |
| --- | --- | --- |
| `bipartite_prefix_dp.cpp` | `be18e9f8c593648c4e2121305f400db5199ed2fcdde46c4522ddbae27b15b5c4` | Byte-identical in `source/`. |
| `resume_aggregate.py` | `8d60d1d1f196a202a76221ff963a5eeb3396d1a827f6c0d64265e20ffd151edc` | Byte-identical in `source/`. |
| `make_shard_plan.py` | `8be79bc530a8a1f38ff4dbf3273e575e126d668ca47026ef89eee6bec34cd310` | Byte-identical in `source/`. |
| `remote_run_template.sh` | `632ff1a2d6df437f1970e38f1e4016560cc99549c2abf36aaff3670f997aa1b4` | Byte-identical in `source/`. |
| `local_cpu_run.py` | `a3c88ff649be0fa59d90af8efebbf049a555db3b8d8281351523e96094d70ed2` | Public SHA `16bb813270dcec33483e241a02925fb1dec204cab43cad9b4931720d404731f0`: one cache-path line replaced by required `CLANG_MODULE_CACHE_PATH`, plus its explanatory comment; no algorithm, parameter, or scheduling field changed. |
| analytic lifting manuscript | `f21f8996341d757f7884efb1ca69ea01903e8e255ed52e2e48a00eb8865812e7` | Public SHA `a236112e57419423b0ddcb8fd50009bcb90040d6780b245b60e0e49725eafbf2`: only the private run-directory reference was changed to `evidence/`. |
| strict aggregate | `6d4fc2f453ef8ef17ec14a58ebe808936516c8d04cf5b408fb73d59679fd3a5e` | `evidence/strict-aggregate.public.json` retains status, source hash, rule, and every coverage origin; only local output/plan paths are removed. |
| production outcomes | `e95ca401e906b4a935c51bbccdcabfbed7a89d983e93b97dc0e5b0c74e29e243` | `evidence/execution-exits.public.json` retains all 325 production exit/task/checkpoint/output-hash/timing records while removing paths and PIDs. |

The counted benchmark adds the 326th completed checkpoint.  Its original
receipt SHA-256 is `82efd875636a8ad054f4a09ea960d9182b7a0976b7a127300da89094a9a115a6`.
Original gate-result hashes are checked by the portable auditor.  Original
receipt hashes and an uncounted isolated resource-profile hash are retained in
`evidence/run-receipts.public.json`.

No source claim depends on a private path.  Redactions are I/O-only: local
paths, cache paths, PIDs, and nonessential binary/log artifacts; the C++ source,
raw completed checkpoints, source SHA, interval endpoints, zero values, gate
results, exit codes, and material timings remain available for checking.
