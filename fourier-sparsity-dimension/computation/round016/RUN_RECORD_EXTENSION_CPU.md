# R16 six-variable extension: CPU integer-FWHT verification

Scope: each of the 896 supplied four-variable relative masks was used for
both blocks, with both values of `a`, giving exactly
`2 * 896^2 = 1,605,632` phases.  This is a finite exact verification of this
specified construction, not a proof beyond the supplied seed family.

## Construction and arithmetic

For lexicographic `(a,i,j)`, the C++ verifier constructs

`g = h_i | (h_i << 16) | (h_j << 32) | (h_j << 48)`,

and for `a=1` xors `0xffff0000ffff0000`.  It xors this with the truth-table
mask of `x0*x2 + x1*x3 + x4*x5`, converts values to exact `+1/-1` integers,
then applies a 64-entry integer FWHT.  Each phase is marked bad if any
coefficient's magnitude differs from 8.  An atomic counter assigns index
chunks to exactly 15 C++ threads; processed-count equality supplies the
coverage check.

## Input and execution evidence

- Remote host/directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r016/extension-cpu`.
- JSON input SHA-256: `c94010a703187ddeccc15f94d3390e7777324cd2d99742a1476872a0b61f3df0`.
  The extractor checked 896 distinct masks in `[0,65535]`, including
  zero, and emitted `seed_masks.txt`, SHA-256
  `16c5b23d9b53ae0d533dd9296735970dd46534b2c04ef3eb16be2a9b280bd745`.
- Extractor source SHA-256:
  `faa9ad6b880f6a1e81da55acfef80db428f1fc5143f244965b168ed574d94647`.
  CPU verifier C++ source SHA-256:
  `6c341da3a336588ffb4e076281b91c89e28ad887fb361bb2cc1057684cf43ed0`.
  Runner source SHA-256:
  `4a94cd436a7953cab5a22c2a0bce999393a87033e145c27626a8c94713548d62`.
- Compiler: Apple clang 21.0.0; Python 3.9.6; binary SHA-256
  `ff1ff173b7a21257484be07534a717a8f0f9c6dd5af1d0e792f982f66dc693e3`.
- One SSH banner exchange timed out before the run command was issued.  A
  process/directory check found no runner or CPU verifier; the next foreground
  launch was therefore the only actual run.  Its child and driver both exited
  0.  It did not interact with the separate GPU task.

## Result

| exact check | value |
| --- | ---: |
| expected / processed phases | 1,605,632 / 1,605,632 |
| bad phases | 0 |
| bad Walsh terms | 0 |
| maximum `abs(abs(W)-8)` | 0 |
| C++ elapsed time | 0.166073 s |
| runner wall time | 1.194532 s |

The only recorded CPU samples are start/end because the owned workload finished
in under the 0.5-second sampling cadence.  They therefore do not evidence
sustained full-machine CPU use.  The post-run query found no owned verifier or
runner process.

Compact outputs are
`/Volumes/data/CodexTemp/T2/a28/r016/verify-result.json` (SHA-256
`262be721ed15df0a5a2a09f192efb8542f7d011efcf1d466c5cb53f2f2c80e60`) and
`run-receipt.json` (SHA-256
`1bb6e3549dc9fb35c5607ccbd88b841c6cabcf295b7bb79bf3017fb8849c2900`).
