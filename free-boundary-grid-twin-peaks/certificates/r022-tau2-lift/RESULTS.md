# r022 validation receipt

## Scope

Only terminal-kernel implementation and small deterministic checks were run.
No production batch was launched.  `tau2_lift.swift` starts from frozen r021
source SHA-256 `fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629`.
The original `metal` and `inverseMetal` string bytes were independently
extracted and compared equal before/after (SHA-256 respectively
`0f93453297b587362fc5fb100568b931ed49dfc928f3de20dedb7448cf15e574`
and `5c75e5b87bf60981b6ac2e9441c20775739381f3e009395dbf3e5a360c6e29ac`).
The added source SHA-256 is
`e268035742b1c8eb7b0e8c366e3b8e393054352ed0755cd523770e1d8d7d9c3a`.

## Actual tests

All following use `tau=2`, stratified mode and zero MH moves, with full trace:

| n | P | domain | audit result |
|---:|---:|---|---|
| 3 | 17 | `9000000000000000` | root/rank/group audit passes; same-seed r021 common fields match |
| 3 | 257 | `9000000000000100` | audit and old-field comparison pass |
| 8 | 257 | `9000000000000800` | audit and old-field comparison pass |
| 32 | 17 | `9000000000002000` | audit and old-field comparison pass at 1024-vertex stride |

`source/verify_tau2_lift_trace.py` independently reconstructs each traced
particle's roots, the least one-based low-root-neighbour rank, L and all four
group sums.  Every reported group passed `LF - 2 LG - 4 LD = 0` exactly.
The same-seed comparison excludes only `gpu_seconds` and `wall_seconds`.
`n3-p17-tau3-inactive.jsonl` verifies every `tau2_lift_*` field is null for an
inactive layer.  No terminal invariant guard fired.

For an actual two-repeat integration check, `n3-p17-r2-new.jsonl` was passed
read-only to the independent strict aggregator under
`/Volumes/data/CodexTemp/T2/a022-tau2-lift-aggregate-r022/`.  It verified
metadata, full replicate/seed coverage, source/raw hashes from the sidecar,
group capacity, integer bounds, exact identities and dead/underflow handling.
Its result is `n3-p17-r2-strict-summary.json`: `Q=0.5610191682068427`,
`V=1.57859709188937`, `C=-0.945989165041426`, with a two-repeat delta SE for
Q of `0.04358`.  This is a tiny finite numerical diagnostic, not a confidence
statement or asymptotic result.

All r022 process names were absent after test completion.  Other agents'
independent GPU jobs may continue; this package starts none.
