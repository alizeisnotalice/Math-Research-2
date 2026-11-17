# Round 005 GPU sign-code search: corrected receipt

Date: 2026-11-17 (Asia/Shanghai).  These files were produced in an independent worktree and integrated into
the round-5 archive.  It is not a theorem or an upper-bound claim.

## Corrected finite search

`mm_affine_flat_gpu_search.py` screened every one of the 11,160 affine
three-flats of `F_2^6` for each of 16 seeded six-variable
Maiorana--McFarland spectra.  The MPS float32 Walsh matrix product was checked
entrywise against an int16 CPU product.  Every span extension then used a
cached, exact Python-integer FWHT for the actual XOR mask; it did **not** ask
whether that XOR was another flat.

With 96 random orders and beam width 32 per seed, all 16 seeds recovered
`dim V=15` (eight dimensions modulo the translation/constant code `C`), with
affine support rank 14.  Every displayed basis was rechecked on all 256 of its
`W` slices by integer FWHT at tolerance zero.  No seed produced `dim V>=16`.
This is only finite evidence for this generator family, not a bound on it or
on general Boolean spectra.

GPU screen time including upload/download averaged 8.31 ms per seed; the CPU
integer matrix screen averaged 30.88 ms.  MPS was therefore useful for the
batched screen; strict certification remained integer CPU arithmetic.

## Required positive and negative controls

For every MM seed, the eight three-flats `{v=a}` passed all 256 exact slices;
their span has dimension 8 and joins `C` to dimension 11 (four new quotient
dimensions).  `triangular-positive-baseline.json` independently reruns the
known support-64 triangular construction at sign-code dimension 15 and affine
rank 14.

The three-mask example `u0v1,u1v2,u2v0` accepts all proper subsets and rejects
their triple combination: its exact inverse-Walsh absolute values are `0` and
`16`.  This validates the full-span predicate.

The earlier 192-seed output in `results.json` is explicitly rejected debugging
output; see `REJECTED-192-SEED-DEBUG.md`.  It imposed the false flat-closure
condition and supports no negative conclusion.

## Reproduction

For the original 16-seed corrected-predicate run (limited old beam depth):

```sh
PYTORCH_ENABLE_MPS_FALLBACK=0 python3 mm_affine_flat_gpu_search_prebeam_9e26ee2b.py --seeds 16 --orders 96 --beam-width 32 --rng-seed 20261117 --out /Volumes/data/CodexTemp/T2/a28/r005/replay-gpu16.json
```

For the corrected bounded beam continuation:

```sh
PYTORCH_ENABLE_MPS_FALLBACK=0 python3 mm_affine_flat_gpu_search.py --seeds 4 --orders 96 --beam-width 32 --beam-prefix 6 --beam-target 9 --rng-seed 20261117 --out /Volumes/data/CodexTemp/T2/a28/r005/replay-gpu4.json
python3 verify_corrected_results.py corrected-results.json --out /Volumes/data/CodexTemp/T2/a28/r005/replay-gpu16-verify.json
```

Both use the same RNG root; run 0 repeats the same MM seed. See BEAM_CONTINUATION.md for the precise bounded search coverage. Rejected raw results.json is intentionally not archived; REJECTED-192-SEED-DEBUG.md records the implementation failure.
