# Canonical inverse-fibre checks

The mathematical statement and its scope are in
`EXACT_INVERSE_FIBRE_DP.md`; the label-gap refinement and deletion-layer
constraints are separate, strict inequalities, not asymptotic estimates.

For fixed output and fixed candidate endpoint, the state graph has at most
`N * 2^(d-1)` states.  Each state considers at most `d` internal-path
insertions plus its next retained vertex, so the implemented rolling DP costs
`O(N*d*2^(d-1))` per candidate.  A square-grid output has at most four
candidates at the last path vertex.  The Python checker implements this same
state graph with dictionaries; it is consequently an exact small-size audit,
not a claim of fast large-n enumeration.

Reproduce the finite checks from this directory:

```sh
python3 exact_inverse_fibre_n3.py > /Volumes/data/CodexTemp/T2/a022-archive-prepare-r020/exact-replay.stdout
python3 gap_audit_n3.py > /Volumes/data/CodexTemp/T2/a022-archive-prepare-r020/gap-replay.stdout
```

The first command exhausts all `9!` words, matches every eligible one-peak
output and distance `d=2,3,4` against its exact inverse fibre, and writes
`exact-inverse-fibre-n3.json`.  The second exhausts the same `9!` words and
checks every `d,g` instance of the low-peak label-gap inequality.

The verified totals are `A1=29568`, `A2=143808`, and inverse-fibre totals
`63552, 64512, 15744` at distances `2,3,4` respectively.
