# Public replay of the R74 diagonal equality

The independently executed historical checker `check_reverse_credit_diagonal.py` read the original local 1200-row raw output and its original receipt paths. Its frozen source, result and receipt preserve that actual execution, but that checker is not the public replay entry point. The regenerable 1200-row output is not required for the mathematical equality certificate and is not included in the public archive.

Use the following public-only replay from the topic directory:

```sh
python3 -B computation/round074/check_reverse_credit_portable.py --summary computation/results/round074-root/reverse-credit-summary.json --output /Volumes/data/CodexTemp/T2/r074-portable-equality.json
```

The portable runner was adapted from the independent checker (source SHA256 `947376361b4d6bc4cb208abed75fee1b0a3ae2558cd8b4d8512483f8c9b74680`). It reconstructs seed742090 from the compact public summary, checks the complete tied center prefix and all three nonempty subsets, then checks the exact clipped singleton quarter areas. It does not read the original receipt, geometry path or 1200-row raw output. This is a portability check of the same strict equality4, not a second independent mathematical proof or a replay of the full search.

Root actually ran this version with Python3.13.18, PID5062, exit0; PID absence was verified, with wall time about0.043s. Frozen public result and receipt accompany it. Existing source/output hashes in the old independent package that refer to the omitted raw output are historical local provenance, not unresolved public-file dependencies.

No general reverse-credit inequality, finite-sample exhaustive optimum or change to C₂ follows.
