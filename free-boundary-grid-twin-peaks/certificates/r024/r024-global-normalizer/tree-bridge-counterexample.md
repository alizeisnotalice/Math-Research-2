# Binary-tree bridge: a finite obstruction and the actual asymptotic target

Let \(D_h\) be two copies of the complete rooted binary tree of height
\(h\), joined by an edge between their roots.  It is a tree of maximum
degree three and has \(N_h=2(2^{h+1}-1)\) vertices.

## A strict counterexample to the coefficient-three tree guess

For \(h=3\), \(N_3=30\), exact integer rank-interleaving/marked-forest
counting gives

\[
 A_2(D_3)=4570935802581417984000000,
\]
\[
 B(D_3)=1382099471396359372800000,
\]
and hence
\[
 A_2(D_3)/B(D_3)=838907595/253657849\approx3.30724082975>3.
\]

Therefore the proposed statement \(A_2(T)\le\Delta(T)B(T)\) is false
even for a degree-three tree.  This is a finite, exact counterexample to
that auxiliary inequality, not a counterexample to the free-grid problem.

The historical independent tree-message audit is retained in the sibling
`r024-tree-root/` archive package; its source-version scope is documented
there.

## Independent height-(2,4) check

For the asymmetric tree made from height-two and height-four complete
binary blocks joined at their roots (38 vertices), an independent Wolfram
Language rational-integral implementation gives

\[
 \frac{A_2}{B}=\frac{263569067611}{74252195533}
 \approx3.5496467911694497423.
\]

The construction is explicit: vertices 0--6 and 7--37 are the two binary
blocks, with the usual parent-child edges and bridge \(0\! -\!7\).
The Wolfram run used exact rational polynomial integration, took 9.66
seconds, had maximum RSS 93.5 MB, and reported zero swap.  Its source and
stdout are `wolfram_tree_audit.wls` and `wolfram_tree_audit.stdout` in this
directory.

## What remains open

For equal blocks, exact values through \(h=7\) are

\[
1.4211,\ 2.2849,\ 3.3072,\ 4.7980,\ 7.1636,\ 11.0520,\ 17.5825.
\]

They motivated the all-height proof now recorded in
`double-binary-b-upper.md`.  It combines the exact tree cut-edge formula,
hook-length root comparisons, and the previously proved general weighted
promotion inequality to give

\[
 {A_2(D_h)\over B(D_h)}
 \ge {2^{h+2}-2\over34560(h+1)}\longrightarrow\infty.
\]

Thus no graph-uniform normalizer inequality \(A_2\le C(\Delta)B\) can
hold even at \(\Delta=3\).  This result remains a comparison-class result:
it does not prove or disprove the free square-grid twin-peaks statement.
