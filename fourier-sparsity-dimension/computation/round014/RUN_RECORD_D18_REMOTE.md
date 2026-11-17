# R14 remote D18 exact-capacity scan

Scope: only `remote_d18_free_pair` and `remote_d18_mixed` from the supplied
eight-point candidate list.  `local_d19` was deliberately not dispatched.

## Remote execution

- Host/work directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r014/d18-parallel`.
- Command used to launch the guarded foreground-child driver:
  `nohup /usr/bin/python3 ./run_d18_remote.py > driver.log 2>&1 &`.
  The driver invokes `capacity_dp_parallel 15 D <D basis masks>` one case at
  a time, waits for each child, imposes a 900-second limit, and stores that
  child's `returncode`, stdout, stderr, elapsed time and CPU samples in its
  JSON receipt.
- Tool versions: Apple clang 21.0.0 (clang-2100.4.1.1), Python 3.9.6,
  arm64 macOS.
- SHA-256: parallel source
  `a8217f49d07d16991c8914933c2a680c17d47519f67b47f3e8b6669860598467`;
  candidate input
  `f139a59424542739cefb04115809e8c4984cd29d3f85df80ad263140b1754de7`;
  driver `567bd5fe6e87c804cb3cf4b9057fb727fde61886742aed358a8114a78aca6b6e`;
  remote summary `fd29943229f0793a50c53c745ed6fd208f18e7c0b93c9a98f99c7d1a328aeba5`.

## Results

The parallel code first reproduced historic serial tests: D=7 had kappa 8
(return code 0, 2.706 s), and D=13 had kappa 16 (return code 0, 2.376 s).
The authorized D=18 results were:

| candidate | result | wall seconds | upper |
| --- | ---: | ---: | ---: |
| `remote_d18_free_pair` | kappa = 24, rc 0 | 19.777 | 24 |
| `remote_d18_mixed` | kappa = 24, rc 0 | 17.270 | 24 |

Thus these exact instances attain, but do not exceed, their supplied upper
bound.  Both runs were shorter than the prescribed 60-second sampling cadence.
Start samples saw the DP process at 129.3% and 146.0% CPU respectively; this
is not evidence of sustained full-machine utilization.  A post-run remote
process query found no matching driver or DP process.  The GPU/MPS job was not
modified.

## Independent local witness checks

The compact remote summary was retrieved to
`/Volumes/data/CodexTemp/T2/a28/r014-d18/remote-summary.json`.  On the local
machine, `verify_eight_point_case.py` checked all partition cells using exact
F2 ranks and independently decomposed the witnesses: free-pair
`8+8+4+4=24`, mixed `12+4+4+4=24`.  Both its serial and parallel checks on
each component returned zero.  Receipts are
`free-pair-local-receipt.json` (SHA-256
`914f69144d63df88c492885efbcbf8b144c72f19196577491642d78882a365bc`) and
`mixed-local-receipt.json` (SHA-256
`b3a00bc9304d5e57f7f33c8c9c05a4887f32a03c4456052a75d3730b3046905b`).

This is finite exact evidence for the two supplied D18 codes, not a proof for
the general selector/addback assertion.
