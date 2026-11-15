# Global top-two exact-peak identity

Let `C(I)=N! P(I is contained in the peak set)` for every independent set
`I` in a finite graph.  For an ordered pair `(u,v)` of distinct vertices,
fix `u` and `v` as the largest and second-largest labels, respectively.  Put
`S={u,v}`.  Every possible remaining peak lies in
`Z_S=V\closed(S)`, because every vertex in `closed(S)\S` has a higher
neighbour in `S`.

For an independent `I` contained in `Z_S`, conditioning on this fixed top-two
prefix leaves the relative order in every comparison defining the event
`I is contained in peaks` uniform.  Thus inclusion-exclusion over possible
remaining peaks gives the exact number of completions with no peak outside S:

\[
h(S)=\frac{(N-2)!}{N!}
 \sum_{I\subseteq Z_S\atop I\text{ independent}}(-1)^{|I|}C(I).
\]

If `u,v` are adjacent, then the second-highest vertex is not a peak and these
completions have exactly one peak.  Conversely, a one-peak labelling has its
two highest vertices adjacent, since otherwise its second-highest vertex is a
second peak.  If they are nonadjacent, both are peaks, so these completions
are exactly the two-peak labellings whose two peaks are the global top two.

Reordering the sum by `I`, with `Z=V\closed(I)`, gives

\[
A_1=\sum_I(-1)^{|I|}C(I)\frac{2e(Z)}{N(N-1)},
\]
\[
A_{2,\mathrm{top2\ nonadj}}
=\sum_I(-1)^{|I|}C(I)
\frac{|Z|(|Z|-1)-2e(Z)}{N(N-1)}.
\]

The second number is generally smaller than `A2`: a two-peak labelling can
have adjacent global top two and a lower second peak.

## Exhaustive checks

`audit_top2_identity.py` checks both formulae on all `9!` labellings of the
3x3 grid and all `4!` labellings of a four-vertex path.  On the 3x3 grid it
returns `A1=29568` and `A_{2,top2 nonadj}=81984`, while `A2=143808`; on the
path it returns `A1=8`, `A_{2,top2 nonadj}=12`, while `A2=16`.

This identity costs only the closed-neighbourhood size and the edge count of
`V\closed(I)` per independent-set representative.  It can therefore be
added to an existing n=7 D4 traversal without retaining another table.  It
does not, by itself, decide the original distance problem.
