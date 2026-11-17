# Provenance

Certified C++ source SHA-256: `e1bc6ec7753f60f65f61da7b332e481c8cbe4613b1a31f6147935da081a52595`.
The source completion-audit SHA-256 is
`f7099d5b90b1cb69e8881408653acf6da7de05ac69015e4a37e5afad924441c3`.
It performed only scalar h comparisons before this independent CRT audit.
Public r033 global h audit and its Pratt certificate are copied as references.
No binary, remote absolute path, large stdout, or production launcher is included.

`references/r018-n7-distance-histogram.json` (SHA-256
`19b0a3672ad8fcd07c94c35388889ccb7a626e05b1d62470548b2c0e841cd01d`)
supplies the independently stored full Euclidean integer. The separate
`evidence/n7-spatial-residuals-audit.json` records the exact subtraction and
explicitly records the unavailable all-L Manhattan numerators.

The owning r018 producer is shipped at
`references/r018-owner-source/n7_d4_peak_dp.cpp`, SHA-256
`24ed7b5f67b38b327d4be051670962fbc473c39a7767bedee7266f6adf782c25`.
Line 114 defines `dist` as Manhattan distance; lines 162--163 serialize its
histogram separately from `squared_euclidean_total`. The residual auditor
asserts `sum(distance_counts)=A2`.

`spatial-small-gate/` preserves a relative source/input/raw tree: 60 newly
run n3/n5 residues, the independent joint-table inputs, and the dual-CRT gate
source. Its receipt reports 15 points, 60 files, and all 25 h/E/far fields per
point checked exactly. It is a small implementation gate, not n7 production.
