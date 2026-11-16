# r023 tau=2 terminal L-weighted root-distance histogram

`tau2_distance_terminal_r023.swift` freezes r022 tau2-lift source
`e268035742b1c8eb7b0e8c366e3b8e393054352ed0755cd523770e1d8d7d9c3a`
and adds only an independent terminal Metal library.  It does not modify the
SMC, RNG, resampling, MH, state, inverse-fibre, or existing tau2-lift kernels.

For each tau=2 terminal particle, it independently recovers
`L = j-2`, where `j` is the earliest one-based rank of a neighbour of the
second (high) root after the two roots.  It bins L at the *original* Manhattan
distance `d(root1,root2)=0,...,2n-2`, and accumulates `L^2`.  New JSON fields:

- `tau2_distance_L_hist_groups`: 256-particle group-by-distance integer bins;
- `tau2_distance_L2_groups`: group sums of L²;
- `tau2_distance_scaled_mean_L_hist`, `tau2_distance_scaled_mean_L2`: each raw
  total divided by P.

The existing `tau2_lift_L_groups` supplies the denominator for direct
L-weighted distance tails: sum any desired high-distance bins and divide by
the sum of all histogram bins.  The second raw moment follows directly from
`L2_groups`.

## Invariants and actual checks

At every group, the host enforces: histogram sum equals the independently
computed r022 lift-L sum; bins d=0 and d=1 vanish; and
`sum L² <= group_particle_count*(N-3)^2`.  The kernel independently rejects
malformed terminal roots/prefixes through `tau2_distance_bad`.

`oracle/n3_tau2_distance_oracle.py` exhaustively reads all 9! n=3 words and
finds 81,984 tau=2 words, with L-histogram `[0,0,63552,64512,15744]`, total
L 143,808 and total L² 358,080.  This is exact finite enumeration, not an SMC
estimate.

`oracle/verify_tau2_distance_trace.py` recomputes L, original root distance,
all group bins and L² from every terminal prefix/state.  It passed actual GPU
traces at n=3 P=17, n=8 P=257 (two groups), and n=32 P=17.  A separate r022
binary reproduced all old non-timing replicate fields exactly at the same
seeds for those three runs; see `results/*old-fields-receipt.json`.

## Reproduce a small trace

```sh
TASK=/Volumes/data/CodexTemp/T2/a022-tau2-distance-terminal-r023
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  "$TASK/tau2_distance_terminal_r023.swift" -o "$TASK/bin/r023"
"$TASK/bin/r023" 8 257 1 2 2 "$TASK/results/replay.jsonl" a4b5c6d7 strat 4 trace
python3 "$TASK/oracle/verify_tau2_distance_trace.py" \
  "$TASK/results/replay.jsonl.trace.json" "$TASK/results/replay.jsonl"
python3 "$TASK/oracle/n3_tau2_distance_oracle.py"
```
