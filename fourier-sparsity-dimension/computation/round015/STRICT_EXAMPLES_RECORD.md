# R15 strict examples extracted from completed R14 data

No T DP was rerun.  The extractor joined the already-generated remote R14 task
and result files to R15's remote `per-task-comparison.jsonl`, then chose the
least row for each positive gap in lexicographic `(lift_dimension,id)` order.
It wrote the full C, B, E and T bases, exact capacities, and the original T
optimal partition to
`/Volumes/data/CodexTemp/T2/a28/r015/strict-examples.json`.

- Extractor source SHA-256:
  `cf2f002630047d68da5bbc3a4e5ef36754e3333d25203522a2e3d8f5398cdfda`.
- Bound input SHA-256 values: tasks
  `d88aaa304d30126a0a1338fa06a10d0a00026abbd68cdca9dd403dba9fd6494c`;
  completed results
  `5221efeba0c8a4435a88c8bf388d42e31f23a29eb1d5e6e17677fcb5208382e6`;
  R15 comparison `127683790d67177ed1cdd98af8458b041737acd674a1c36613ebe824a1c26ebd`.
- Output SHA-256:
  `827489628cd34856a65c543df0d25f951e90e8b7c89e0bc04c4c85c768899af4`.
- Remote command returned 0; postcheck found no extractor or audit process.

The selected strict rows are:

| gap | R14 ID | lift D | kappa(C) | kappa(B) | kappa(T) | RHS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1509 | 8 | 6 | 4 | 9 | 10 |
| 2 | 10530 | 10 | 8 | 6 | 12 | 14 |

They are finite exact witnesses that this upper bound need not be attained on
every supplied R14 instance.  They do not challenge the proved upper bound or
assert a general sharpness classification.
