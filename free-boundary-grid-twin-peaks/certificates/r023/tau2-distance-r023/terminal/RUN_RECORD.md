# Actual r023 check record

Source SHA-256: `01b77bca1735f07293eae5ffeea0326810fdc5ea785c7a9e99d9030cba9b0702`.
Frozen source input SHA-256:
`e268035742b1c8eb7b0e8c366e3b8e393054352ed0755cd523770e1d8d7d9c3a`.
Compiled binary SHA-256:
`06d78d3b5c1071521e9cce369f420e394c22448a91c6b7e2372886603f509977`.

On Apple M6 Metal, macOS 27.2, Swift 6.4.2, actual terminal traces were run
at n=3 P=17 tau=2 L=16, n=8 P=257 tau=2 L=4, and n=32 P=17 tau=2 L=4.
All were independently replayed from their saved terminal prefixes/states;
the n=8 run exercises a 256-particle and a one-particle group.  The n=32 trace
is 159 MB and remains only in this temporary calculation directory.

The exact n=3 9! oracle, trace receipts, and same-seed r022 comparisons are
stored in `results/`.  All commands completed; no r023 GPU process remains.
