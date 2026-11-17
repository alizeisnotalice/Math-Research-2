# R18 quadratic graph-cases native Metal receipt

Remote task directory: `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r018/graph-metal`.
The foreground command was `./run.sh >run.stdout 2>run.stderr`; it completed
under `set -euo pipefail`.  A subsequent process query found no `verifier` or
`swiftc` task, so the Metal process, its queue, and its buffers had exited.

## Bound artifacts

- Source SHA-256: `900a71abd506c3a1e6d10d748083c0534cf9a7dfb54c4c212827a1ca90a00cf4`
- Input SHA-256: `52d992aa77f2075575fee144069602babc54266e03b7ee3deec48ed9ba9daa81`
- Binary SHA-256: `db11f1134ccbdb5100e925a1594f59a9fbf23c2e610b83ac976f1839c5c7ce80`
- Result SHA-256: `72a3edfc1bf50b563b22fd1ffc8143c1d7eb6e57bf8cd86d2f3c7ee5924fd919`

Environment metadata in `environment.txt` records Apple Swift 6.4.2 on macOS
27.3 arm64; the report records device `Apple M5 Pro`.

## Actual check

The repaired source uses truth-table length `N=2^n`, not the variable count
`n`, for hex limbs, rank validation, CPU FWHT, sign/output arrays, and GPU
threadgroups.  It checked 512 phases at `N=64`, `q=8` and 16,384 phases at
`N=256`, `q=16`.  Every GPU coefficient had the requested magnitude and CPU
int64 FWHT agreed with every GPU signed coefficient: exactly
`512*64 + 16384*256 = 4,227,072` comparisons.

The original signal 133 was caused by confusing the variables `n=6,8` with
truth-table lengths `64,256`; it occurred before any valid result.  The final
source above was rebuilt and rerun successfully.

The 38.518124 ms report total is host-side command-buffer submission, wait,
and shared-buffer readback time only.  It excludes CPU phase generation and
CPU comparison; it is not a GPU timestamp, kernel-only duration, utilization,
or speedup claim.  Maximum GPU output buffer was 2 MiB.  The reference array
used roughly 4 MiB at the largest block, but this is not a measured peak RSS.
