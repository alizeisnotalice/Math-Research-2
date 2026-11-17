# Incomplete n=7 attempt

This is a minimal failure/partial-evidence record. `n7-v24-p0.json` is the
single completed first-prime checkpoint. The second pass exited nonzero and
did not write its checkpoint; no later pass was started. The saved resource
logs are evidence of those two process outcomes only. They do not diagnose
the failure and cannot be combined into a CRT value.

`local-n7-exit-check.json` records that the six listed local driver/child
PIDs were absent after cleanup.  It deliberately records the p1 failure cause
as unknown; it does not identify an OOM, jetsam, or software cause.
