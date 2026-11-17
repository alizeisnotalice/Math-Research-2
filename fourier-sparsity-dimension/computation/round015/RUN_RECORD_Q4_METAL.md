# R15 Q=4 native Metal run receipt

Remote directory: `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r015/q4-metal`.
Command: `./run_q4_metal.sh >run.stdout 2>run.stderr`.

The foreground command completed successfully under `set -euo pipefail`.
`run.stdout` records `status=passed rows=3 gpu_ms=69.623661 bad=0 cpu=98304`;
`run.stderr` contains only the Swift driver version emitted while the script
captured environment metadata.  A completion query found no matching
`q4_metal_walsh_verifier` or `swiftc` process, so this run's Metal process and
its device/queue/buffers had exited.

Environment is specific to this R15 run: Swift 6.4.2 on macOS 27.3 arm64;
`PYTORCH_ENABLE_MPS_FALLBACK=0`.  The accompanying manifest binds source
`26fba2c7...`, input `cec2b5f1...`, binary `fee108ff...`, and report
`c0314d3c...`.  It must not be replaced with the R14 receipt.

The report records Apple M5 Pro Metal, three rank-15 codes, full CPU int64
coefficientwise cross-checking of all 98,304 phases, no Walsh failures, and
stabilizer count 2 for every code.  Its reported time is the host-side
command-buffer submission/wait/readback interval, not a GPU timestamp or
kernel-only measurement.
