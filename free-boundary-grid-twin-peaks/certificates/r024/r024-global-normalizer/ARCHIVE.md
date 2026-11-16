# r024 global top-two normalizer audit

For a connected graph, let \(A_2\) count all two-birth orders and \(B\)
count those in which the first two vertices are births. The files here retain
exact, auxiliary results only; none is a free-boundary square-grid result.

`double-binary-b-upper.md` is the final all-height proof (SHA-256
`c20a6b13a336e58ab37b1301f2eaba837e871531dfd34f29a4bac1fe93c17421`).
For the degree-three tree \(D_h\) made by joining two height-\(h\) complete
binary roots, it proves

\[
B(D_h)\le720(h+1)A_1(D_h),\qquad
\frac{A_2(D_h)}{A_1(D_h)}\ge\frac{N_h}{48}\quad(h\ge2),
\]

and hence \(A_2(D_h)/B(D_h)\ge N_h/[34560(h+1)]\to\infty\).  Thus no
bound \(A_2\le C(\Delta)B\) can hold for every finite connected graph even
when \(\Delta=3\).  The proof includes the cut-edge identity, the ancestral
hook correction \(G\), root-ratio bounds, and the selected-subtree-size sum.
It depends on the general weighted-promotion theorem, retained verbatim as
`weighted-promotion.md` (SHA-256 `26c13189…`).  `verify_constants.py` and
`constants-audit.json` are finite rational sanity checks only, not a
substitute for that all-height argument.

1. `path-cycle-normalizer.md` proves \(A_2/B<2\) for every nondegenerate
path and cycle by rank-1 insertion. This is a degree-two calibration only.
2. `leaf-bridge-identity.md` proves the one-leaf insertion identities and
their average over all attachment vertices. It exposes a second \(L\)-moment
but does not give a bounded-degree or grid estimate.
3. `exact_kernel.py` exhaustively scans all connected labelled graphs on six
vertices; `exhaustive-n6-all.json` records the maximum \(9/5\).

`wolfram_tree_audit.wls` and its retained output independently evaluate the
same integer rank-order recurrence on a 38-vertex degree-three tree and give
\(A_2/B=263569067611/74252195533\approx3.54965\). It is an independent
implementation check of the recurrence, not a proof for all trees. The
separate direct-integer N=30 counterexample lives in `../r024-tree-integer/`.
Neither tree result is a free-boundary grid counterexample.

`tree-bridge-counterexample.md` is a scope summary updated after the
all-height proof. `normalizer-audit.md` is deliberately excluded: its former
"no counterexample found" status was superseded by the theorem above.
