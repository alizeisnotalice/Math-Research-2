# R17 nonparallel-case native Metal receipt

Remote path: `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r017/nonparallel-metal`.
The foreground `./run_nonparallel_metal.sh >run.stdout 2>run.stderr` completed
under `set -euo pipefail`; `run.stdout` says `status=passed phases=3072` and
`run.stderr` only contains the Swift-driver version recorded during environment
capture.  A completion query found no verifier or `swiftc` process, so this
run's Metal context and buffers were released.

The report/manifest bind input SHA-256 `f6496c64...`, source SHA-256
`8809f77b...`, binary SHA-256 `a6bb4f1e...`, and report SHA-256 `7f7a7f25...`.
It used Apple M5 Pro Metal, unified memory, a 32-wide execution width, and a
maximum 2,048-phase/524,288-byte output block.

For both specified cases, the verifier checked GF(2) basis rank, and that the
given blocks are a disjoint partition of the 64 truth-table points.  It then
checked all 3,072 phase masks: all GPU unnormalized Walsh entries had magnitude 8,
and an independent CPU int64 FWHT compared every one of all 64 coefficients
for every phase with the GPU output.  This block test is **only** truth-point
coverage/disjointness; it is not a local Walsh-support or branch check.  Such
branch analysis remains the separate Python computation recorded in the input.

The 4.390501 ms total is host command-buffer submit/wait/readback time, not a
GPU timestamp, pure kernel duration, performance claim, or utilization value.
