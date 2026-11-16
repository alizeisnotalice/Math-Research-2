# r012 GPU parent-selection baseline

`tau_smc_integer_gpu_parent.swift` preserves the previous multinomial
resampling random streams: each child uses the existing `childkey` domain 1
for the integer target and domain 2 for mutation.  The host still constructs
the exact group CDF from `groupSums`, but it no longer loops over children or
uses a host `Set`: the GPU binary-searches the CDF and scans only its selected
256-particle group.  `min_distinct_parents` is intentionally null: a one-step
parent count is not an initial-ancestor (genealogical) count.

Completed numerical baseline (not a proof): with exactly the prior n=8
parameters and seed schedule, the old and new weighted ratio agree bit-for-bit
in printed output, `0.3369305512539116`.  The new end-to-end wall time was
109.805 seconds versus 212.987 seconds in the prior timed run; hardware GPU
seconds changed from 66.176 to 68.171.  No AGX time series was collected for
this run, so it makes no utilization claim.

The proposed independent future mode is integer stratified resampling.  Given
integer weights summing to `T`, child `j` independently samples
`t_j in {0,...,T-1}` and selects CDF target
`floor((j*T+t_j)/P)`.  This must be maintained separately from the audited
multinomial baseline.  In particular, it is not the shared-offset systematic
scheme.
