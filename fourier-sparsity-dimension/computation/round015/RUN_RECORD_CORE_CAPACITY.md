# R15: exact core-capacity audit of completed R14 records

Scope: this is a join against the already completed R14 `tasks.jsonl` and
`results-rerun.jsonl` files.  It does **not** rerun the 29,085 lifted T DPs,
and it does not run any H=C+wB capacity computation.

## Method

For each distinct supplied core basis `B_basis`, the driver enumerates its
F2 span on its original at-most-seven-point domain.  It applies the exact
support-partition recurrence used by the capacity DP, with integer containment
counts, verifies the recovered partition by an independent sum of
`subcode_dimension + 1`, then joins its result to every matching completed T
row.  The checked comparison is

`kappa(T) <= base_capacity (= kappa(C)) + kappa(B)`.

## Run evidence

- Remote host/directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r015/core-capacity-audit`.
- Actual foreground command:
  `/usr/bin/python3 ./audit_core_capacity.py --tasks
  /Users/zhengzhihao/CodexRemoteMath/a28-selector/r014/tasks.jsonl --results
  /Users/zhengzhihao/CodexRemoteMath/a28-selector/r014/results-rerun.jsonl
  --out-dir ./out --workers 15`.
- First attempt exited 1 before results because macOS Python 3.9 does not
  offer `int.bit_count`; the source was changed to the exact equivalent
  `bin(value).count("1")`, locally smoke-tested, copied, and the subsequent
  foreground driver returned exit 0.
- Final source SHA-256:
  `8f4314a78bf934131866f52cc7b2803c9b406e4fa867b7b3d5a227f4b379ac45`.
  Input SHA-256 values: tasks
  `d88aaa304d30126a0a1338fa06a10d0a00026abbd68cdca9dd403dba9fd6494c`;
  completed results
  `5221efeba0c8a4435a88c8bf388d42e31f23a29eb1d5e6e17677fcb5208382e6`.
- Version: remote Python 3.9.6 on macOS arm64.  The 15-worker calculation
  took 1.874993125 seconds, so it did not sustain a 60-second load sample.
  A post-run process check found no audit worker.

## Exact finite result

All 29,085 completed T records and 1,106 distinct B codes were covered:

| `kappa(C)+kappa(B)-kappa(T)` | count |
| ---: | ---: |
| 0 | 26,565 |
| 1 | 1,680 |
| 2 | 840 |

There are zero violations of the tested inequality and 2,520 strict cases.
Consequently equality is frequent but is not universal in this finite R14
family.  `minimal-counterexample.json` records `exists: false`; its only
meaning is no counterexample within this supplied finite input, not a general
theorem.

The raw per-task comparison stays remote at
`out/per-task-comparison.jsonl` (3.2 MiB; SHA-256
`127683790d67177ed1cdd98af8458b041737acd674a1c36613ebe824a1c26ebd`).
Compact results were retrieved to
`/Volumes/data/CodexTemp/T2/a28/r015/summary.json` (SHA-256
`fb2604d2f3ba51efd123413fa15560d84be850f4fe298339d87d3355ea3a12ae`) and
`minimal-counterexample.json` (SHA-256
`1967ff1e84d6ff1ae24424b3482351dce35472eb28e28510b2d3ea0ac852f76f`).

This is an exact finite calculation supporting/ruling out only the displayed
comparison on R14's task family.  It does not establish the proposed general
non-unital capacity identity.
