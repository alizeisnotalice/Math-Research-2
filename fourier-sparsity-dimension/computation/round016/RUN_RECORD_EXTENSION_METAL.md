# R16 nonlinear extension-family native Metal receipt

Remote directory: `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r016/extension-metal`.
The foreground command was `./run_extension_metal.sh >run.stdout 2>run.stderr`.
It returned successfully under `set -euo pipefail`; `run.stdout` records a
passed 1,605,632-record run and `run.stderr` only contains the Swift driver
version emitted while capturing environment metadata.  A completion query
returned no `extension_metal_verifier` or `swiftc` process, so the task's Metal
context and buffers were released.

Input SHA-256: `c94010a703187ddeccc15f94d3390e7777324cd2d99742a1476872a0b61f3df0`.
Source SHA-256: `d75365e2401d76b6468f77b8728e5f0ecabf0bcf59bece338201ad6e0729a6f9`.
The manifest additionally binds the native executable and report.

The host generated, without assuming linear-code closure, every ordered
`(a,h0,h1)` with `a in {0,1}` and `h0,h1` among all 896 seed masks.  Thus all
`2*896^2 = 1,605,632` requested phases were covered in 784 chunks of at most
2,048 phases.  The first full chunk (2,048 phases x 64 coefficients) agreed
entrywise with CPU int64 FWHT before the sweep.  All GPU output coefficients
then had absolute value 8; the report records every completed block, zero bad
entries, and the source-defined rolling checksum `8c8eb544649d0383`.  Its
offset is `1469598103934665603`, so it is neither the conventional FNV-64
offset nor a cryptographic or independent full-coefficient certificate.  The
CPU/GPU entrywise claim is limited to the first 2,048 phases; a separate CPU
process is responsible for any full-family independent amplitude check.

Apple M5 Pro Metal was used with unified memory and 32-wide execution.  The
512 KiB output-buffer cap was respected.  `total_host_ms=1367.080665` is the
sum of host command-buffer submit/wait/shared-readback intervals, not a GPU
timestamp, pure kernel duration, speedup, or utilization measurement.
