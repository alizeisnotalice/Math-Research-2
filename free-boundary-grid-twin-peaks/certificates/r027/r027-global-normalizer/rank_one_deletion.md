# Rank-one deletion: a localized exact normalizer identity

Let \(G=(V,E)\) be a finite connected simple graph.  A rank word is written
in decreasing-label order, so its first vertex is the global maximum and a
peak is a birth.  Write \(A_2(G)\) for the number of two-peak words and
\(B(G)\) for the number of such words whose first two vertices are both
peaks (equivalently, the top-two sector used in the first-merge identity).

For \(v\in V\), set \(H_v=G-v\).  If \(\pi\) is a rank word of \(H_v\),
let \(P_v(\pi)\) be its peak set and let \(r_v(\pi)\) be its first vertex.
Put
\[
 \mathcal C_v=\{\pi:\; |P_v(\pi)\setminus N_G(v)|=1\}.
\tag{1}
\]
Then the following two identities are exact:
\[
 A_2(G)=\sum_{v\in V}|\mathcal C_v|,\qquad
 B(G)=\sum_{v\in V}|\{\pi\in\mathcal C_v:r_v(\pi)\notin N_G(v)\}|.
\tag{2}
\]

## Proof

Delete the first vertex \(v\) of a two-peak word.  It was a peak.  In the
remaining word, precisely the peaks adjacent to \(v\) can be newly exposed:
every vertex outside \(N_G(v)\) has exactly the same earlier-neighbour
relations as before.  Conversely, prepend \(v\) to any word of \(H_v\): it
is a new peak and it suppresses exactly the peaks of that word lying in
\(N_G(v)\).  Thus the resulting word has exactly two peaks if and only if
(1) holds.  These operations are inverse after recording \(v\), proving the
first equality.

The second vertex of the prepended word is \(r_v(\pi)\).  It is a peak in
the new graph exactly when it is not adjacent to \(v\); if so, the two first
vertices are the two peaks.  This proves the second equality.

For a square grid, \(\deg v\le4\), so every \(\pi\in\mathcal C_v\) has at
most five peaks:
\[
 |P_v(\pi)|=1+|P_v(\pi)\cap N_G(v)|\le5.\tag{3}
\]
Consequently the desired uniform normalizer estimate is exactly the local,
punctured-grid root-transport inequality
\[
 \sum_v |\mathcal C_v|
 \le C\sum_v |\{\pi\in\mathcal C_v:r_v(\pi)\notin N_G(v)\}|.\tag{4}
\]
Equivalently, under the uniform measure on pairs \((v,\pi)\) with
\(\pi\in\mathcal C_v\),
\[
 {B(G)\over A_2(G)}=
 \Pr\{r_v(\pi)\notin N_G(v)\}.\tag{5}
\]
This formulation is sharper than merely bounding all \(\le5\)-peak words of
\(G-v\): it keeps the requirement that all but one peak lie in the four-site
neighbourhood of the puncture.

## Exact stable-promotion fibres

Let \(\mathcal D_v\) be the good subfamily on the right of (2).  In a word
\(\zeta\in\mathcal D_v\), its root \(p\) is the unique peak outside
\(N(v)\).  Let \(J(\zeta)\) be the one-based position of the first
neighbour of \(p\) in \(\zeta\).  There is an exact, disjoint fibre formula
\[
 \boxed{\quad |\mathcal C_v|=
 \sum_{\zeta\in\mathcal D_v}(J(\zeta)-1),\qquad
 A_2(G)=\sum_v\sum_{\zeta\in\mathcal D_v}(J(\zeta)-1).\quad}\tag{6}
\]
Thus \(A_2/B\) is exactly the \(\mathcal D_v\)-average of \(J-1\).

For the proof, take \(\pi\in\mathcal C_v\), let \(p\) be its unique
outside peak, remove \(p\) from its old position, and put it first.  This
stable promotion preserves the whole peak set: vertices not adjacent to
\(p\) retain their relative neighbour order, and every neighbour of \(p\)
was already lower than \(p\).  It gives a word \(\zeta\in\mathcal D_v\).
Conversely, remove the first point \(p\) from \(\zeta\) and reinsert it at
any position \(t=1,\ldots,J(\zeta)-1\).  It stays above every neighbour.
When \(t>1\), the new first point is the old second point.  It was already a
peak (it is nonadjacent to \(p\)) and therefore lies in \(N(v)\); no outside
peak is newly created.  These inverse positions are exhaustive because a
peak must precede its first neighbour.  This proves (6).

The fibre is not pointwise bounded by the hole degree.  For every \(n\ge4\),
take the opposite corners \(v=(0,0)\), \(p=(n-1,n-1)\).  Delete \(v,p\) and
the two neighbours of \(p\).  The remaining induced graph is connected and
contains \(a=(1,0)\): its \((n-1)\times(n-1)\) lower-left rectangle minus a
corner is connected, and each remaining boundary-strip point attaches to it.
Start at \(a\), use any connected ordering of that graph, and finally append
the two neighbours of \(p\).  Prepending \(p\) gives
\(\zeta\in\mathcal D_v\) with
\[
 J(\zeta)=n^2-2,\qquad J(\zeta)-1=n^2-3.\tag{7}
\]
This is sharp for a corner root: a word on \(G-v\) has \(N-1\) positions
and the two neighbours of \(p\) must occupy positions at most \(N-1\), so
always \(J\le N-2\).  Thus the example realizes the maximum possible fibre
length rather than merely a long special case.
So a bounded average in (6), if true, must use global weights; it cannot come
from degree or one-hole local geometry alone.

## What it does and does not solve

Equation (4) would give \(A_2/B\le C\), hence a bounded mean first-merge
Rao--Blackwell weight.  It is not proved here.  In particular, a pointwise
comparison inside a fixed \(\mathcal C_v\) is false as a method: the unique
peak outside \(N(v)\) need not be the root, and lifting it to the root changes
which nearby vertices are suppressed.  The identities isolate the missing
task to a local root-transport statement on a one-vertex-punctured grid.

There is already a literal \(3\times3\) obstruction to the most tempting
swap.  Number the grid row by row.  Take \(v=0\), and on \(G-v\) take the
decreasing word
\[
 (1,2,3,4,6,8,5,7).
\]
Its peaks are \(\{1,3,8\}\); the only one outside \(N(0)=\{1,3\}\) is
\(8\), while the root is the bad neighbour \(1\).  Swapping the root and
the outside peak gives
\[
 (8,2,3,4,6,1,5,7),
\]
whose peaks are \(\{2,3,8\}\), so it has two peaks outside \(N(0)\) and
has left \(\mathcal C_0\).  Thus any proof of (4) must repair the newly
exposed neighbours of the demoted root; simple peak/root swapping is not a
valid injection.

`verify_rank_one_deletion_n3.py` independently exhausts the two small grids.
For \(3\times3\), it gives \(A_2=143808\), \(B=81984\), and verifies (2)
word by word.  The local success fractions according to the degree of the
deleted maximum are respectively \(193/312,200/361,17/38\) for degrees
\(2,3,4\).  These finite fractions are checks, not a uniform lower bound.
