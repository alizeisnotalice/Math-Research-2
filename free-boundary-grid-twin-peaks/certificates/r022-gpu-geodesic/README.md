# r022 terminal ordered-geodesic-prefix kernel

`geodesic_prefix_r022.swift` is an isolated extension of the fc45 n<=32
terminal implementation.  It neither changes the SMC Metal library nor its
RNG, resampling, Metropolis moves, state buffers, or inverse-fibre kernels.
At the terminal step only, `geodesicPrefixGroups` reads the existing ordered
prefix and, for every particle and every `d=1,...,2n-1`, sets an indicator iff
every consecutive pair in positions `1,...,d` is adjacent in the grid and the
two endpoints have Manhattan distance `d-1`.  Both clauses are checked; an
endpoint-only test is not used.

Each GPU thread sums these 0/1 indicators over one consecutive 256-particle
group.  The JSONL replicate fields are:

- `geodesic_prefix_group_counts`: group-by-d exact integers;
- `scaled_geodesic_prefix_p`: their total divided by the final allowed total.

For a non-dead terminal run every final allowed count is one (the only missing
vertex is adjacent to the occupied grid), hence this denominator is `P` and
the latter field is the particle average.  Using the terminal total also keeps
the normalization explicit alongside the established terminal statistics.

## Actual checks

- `oracle/n3_onepeak_oracle.py` independently enumerates all `9!` ordered
  words and finds exactly 29,568 one-birth words.  Its geodesic counts for
  `d=1,...,5` are `[29568,29568,14784,4000,448]`.
- `oracle/verify_geodesic_trace.py` recomputes every particle's indicator
  directly from the saved prefix and checks all 256-group totals and scaled
  fields.  It passed actual GPU traces for n=3 P=17 and P=257, n=8 P=17 and
  P=257, and n=32 P=17.
- The same verifier also passed a two-peak terminal trace at n=8, P=257,
  tau=2.  Its d=2 count was zero, as required because the two prescribed
  roots are non-adjacent.
- `oracle/compare_legacy_replicate.py` proved every pre-r022 replicate field
  (excluding unavoidable timing) exactly unchanged against fc45 for n=3 and
  n=8 at the same seeds.

The results are verification evidence for this computation, not a proof about
the asymptotic conditional law.

## Reproduce

```sh
TASK=/Volumes/data/CodexTemp/T2/a022-gpu-geodesic-prefix-r022
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  "$TASK/geodesic_prefix_r022.swift" -o "$TASK/bin/geodesic_prefix_r022"
"$TASK/bin/geodesic_prefix_r022" 8 257 1 65 65 \
  "$TASK/results/replay.jsonl" a4b5c6d7 strat 4 trace
python3 "$TASK/oracle/verify_geodesic_trace.py" \
  "$TASK/results/replay.jsonl.trace.json" "$TASK/results/replay.jsonl"
python3 "$TASK/oracle/n3_onepeak_oracle.py"
```
