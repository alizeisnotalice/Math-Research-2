# r014 corrected GPU SMC calibration

This is numerical calibration of a fixed-prefix completion kernel.  It is
neither an estimate of the uniform labeling law conditioned on two peaks nor
a proof about the original twin-peaks problem.

## Correction and audit

The r013 propagation kernel copied a parent RNG state after resampling.  In
`corrected_smc.swift`, every child instead receives a deterministic hash of
the parent state, child index, and stage before the uniform-boundary draw.
Thus children no longer simply inherit the identical parent RNG state; this is not a proof of collision-free streams.  The streams
are deterministic and the change is not an assertion of probabilistic
independence.

The inherited CPU audit had also hard-coded `n*n-4` transitions.  It now uses
`n*n-initialCount`, which is necessary because extended-L has `2*n-1` initial
vertices.  For each of n=7,8,9, a 64-particle CPU/GPU audit compared the
normalizing log and event fraction within 1e-12, and minimum surviving root exactly; all three
passed before the production repetitions began.

## Production run

`production-foreground-262144x512/smc.jsonl` contains 512 repetitions at
262144 particles for each shape.  Aggregation is

`sum(exp(logZ)*p_hat) / sum(exp(logZ))`.

| shape | scaled estimate | exact finite value | difference | Z-ESS / 512 | maximum Z share |
| --- | ---: | ---: | ---: | ---: | ---: |
| n7 | 3.4295308891 | 3.4361735254 | -0.0066426363 | 511.93 | 0.002032 |
| n8 | 4.0640195531 | 4.0552255355 | 0.0087940176 | 511.69 | 0.002079 |
| n9 | 4.6349840682 | 4.6562852063 | -0.0213011381 | 510.44 | 0.002312 |

The run used a managed foreground session and exited normally.  Summed Metal
kernel time was 16.822, 64.140, and 87.932 seconds respectively; summed
per-stage submit/wait time was 33.447, 90.756, and 135.384 seconds.  These are
not complete end-to-end wall times.  Minimum final genealogical-root counts
were 43512, 22601, and 10928.

The production JSON SHA-256 is
`b883004d0128ac368430e964138299db506f49aa1010d57db72eaafc1acd4ba5`.

## AGX sample

The separate 32-repetition, 262144-particle mini run recorded 43 unprivileged
`AGXAccelerator` samples.  Device utilization was mean 77.02%, maximum 93%,
with 41 samples at least 70% and 3 at least 90%.  Renderer utilization was
mean 55.91%, maximum 76%.  Two samples were zero near work boundaries.  These
are system-wide instantaneous measurements, not exclusive or full-duration
utilization claims.  The raw record is `agx-mini.txt`.

## Reproduction

Compile with Swift 6.4.2 and Metal runtime, using the data-volume module
cache:

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules \
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
corrected_smc_262k512.swift -o corrected_smc_262k512
./corrected_smc_262k512 NEW_OUTPUT_DIRECTORY
```

The source records the fixed inputs, particle count, repetitions, and seeds.
It checks every Metal command-buffer completion.  Do not overwrite the saved
output when reproducing.

归档说明：SOURCE-RUN-HASHES.txt保留原交接/运行哈希及历史路径，主Agent已核对原文件。归档对说明文字作适用范围修正后，当前文件以课题manifest-r010.json为准；不把历史绝对路径当可移植复现入口。
