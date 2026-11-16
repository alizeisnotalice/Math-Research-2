# Completed CPU root-moment batch: n=11

This directory is a self-contained record of the completed finite CPU SMC
batch for the original two-peak law.  Parameters are `n=11`, `P=16384`,
`R=64`, `L=16`, 12 CPU threads, and seed base `2686592768`; see the metadata
JSON for definitions.  It is not an exact count, a mixing proof, or a
confidence certification.

`covsmc.cpp` is byte-identical to the audited r019 source (SHA-256
`6f1031104b964649ced764cf58b5f6f738063987b985333e10d7eaad411c82de`).
The raw CSV has SHA-256
`2072a396d5e8c895f8076c5fc5bf3a377e3d0f5b67969d0f416589732f757a9f`.

Recompute the Decimal-80, zero-inclusive ratio-of-stratum-means summary:

```sh
python3 summarize_covariance.py n11-p16384-r64-l16.csv > recomputed-summary.json
cmp recomputed-summary.json n11-summary.json
```

The saved runtime record reports 3532.64 seconds wall time, 12-thread CPU
use, 97,533,952 bytes maximum RSS, and zero swaps.  The summary gives
`Q=0.317702275832959056...` with a stratified delta diagnostic
`0.002075136169805317...`; this diagnostic does not establish mixing or a
probabilistic coverage guarantee.
