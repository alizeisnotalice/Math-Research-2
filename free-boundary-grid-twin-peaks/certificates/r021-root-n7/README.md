# Exact n=7 one-peak root vector

`root_counts_d4.cpp` computes all 49 specified-root one-peak counts by a D4-compressed independent-set inclusion-exclusion recurrence using exact U256 arithmetic. The completed n=7 calculation is finite and exact. It does not prove the displayed transport inequalities for all n.

The local verifier uses `reference/n3-roots.json` and `reference/n4-roots.json` for entrywise calibration, and `reference/n7-d4-root-moments.json` for independent r018 values of A1, A2, and the squared Euclidean two-peak total. Thus all inputs are local and relative.

```sh
python3 verify_root_counts.py
```

This reruns 49 root-sum checks, n3/n4 vector comparisons, D4 symmetries, 42 strict inward coordinate inequalities at n=7, and the two finite transport-ratio checks. The completed n=7 resource record gives 129.52 wall seconds, 2,585,559,040-byte peak RSS, zero swaps, and normal exit. The compiled executable is intentionally omitted.
