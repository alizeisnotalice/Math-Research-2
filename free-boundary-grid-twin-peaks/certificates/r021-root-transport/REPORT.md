# r021: root transport is an exact reversible edge chain

Let \(G=(V,E)\) be a finite connected graph.  For every nonempty connected
set \(S\), let \(h(S)\) be the number of legal connected completion words
after the vertices of \(S\) have already appeared.  Thus \(h(V)=1\), and

\[
h(S)=\sum_{w\in\partial S}h(S\cup\{w\}). \tag{1}
\]

Write \(a_x=h(\{x\})\), the number of one-peak words with root \(x\).
For an edge \(xy\), put \(c_{xy}=h(\{x,y\})=c_{yx}\).  Specializing (1) at
a singleton gives the exact root-edge identity

\[
a_x=\sum_{y\sim x}c_{xy}. \tag{2}
\]

Consequently

\[
P(x,y)=c_{xy}/a_x\quad (x\sim y) \tag{3}
\]

is a reversible Markov kernel with stationary measure proportional to
\(a_x\), because \(a_xP(x,y)=c_{xy}=a_yP(y,x)\).  This is not an auxiliary
tree law: it is the exact first-two-vertices law under the uniformly chosen
one-peak word.

There is also a useful exact reverse-time recursion.  If
\(a_x(H)\) denotes the number of connected words of a connected graph
\(H\) beginning at \(x\), then splitting by the last vertex gives

\[
 a_x(H)=\sum_{\substack{v\in V(H)\setminus\{x\}\\H-v\ \mathrm{connected}}}
 a_x(H-v). \tag{4}
\]

Indeed, deleting the last letter of a connected word leaves a connected
induced graph, and conversely a connected word of \(H-v\) may be followed
by \(v\).  Thus \(a_x\) also counts deletion orders which keep the
unremoved induced graph connected and end at \(x\).  A square grid has no
cut vertex, so its first recursion step sums over every \(v\ne x\).
This does **not** immediately induct root unimodality: after one deletion
the graph has a hole and later deletions can create cut vertices.  Any
induction must therefore state and prove a comparison for this larger class
of punctured connected subgraphs, not merely for rectangles.

This is a real obstruction to a naive heredity proof of the sharper
square-only estimate below.  In \(P_5\square P_5\), the four vertices

\[
(1,0),(2,0),(2,1),(3,1)
\]

induce a \(P_4\).  Its four root counts are \(1,3,3,1\), by the direct
interleaving count for a path.  Thus along the first horizontal edge the
ratio is \(3\), whereas the square candidate at \(n=5,i=1\) would have
right side \(9/4\).  Every connected induced subgraph can occur as a
remaining set in some connected deletion order (add vertices back in a
connected breadth-first order), so this is not an artificial state.
The file audit_deletion_state_obstruction.py exhaustively checks all
\(4!\) words.  It does not refute the estimate for the *full* square, but
it rules out an induction which assumes that estimate for every deletion
state.

## What an elementary root move really proves

For an edge \(xy\), take an \(x\)-rooted word and move the occurrence of
\(y\) to the first position.  The resulting word begins \(y,x\).  Given an
output beginning \(y,x\), recovering a preimage amounts to inserting \(y\)
back into one of at most \(N-1\) positions after \(x\); only some such
insertions are legal.  Hence the general sharp elementary comparison is

\[
\frac1{N-1}a_y\le a_x\le (N-1)a_y. \tag{5}
\]

It cannot be improved uniformly over bounded-degree graphs: on the path
\(P_N\), \(a_i=\binom{N-1}{i-1}\), so the endpoint/next-vertex ratio is
\(N-1\).  Thus a constant-factor adjacent-root transport requires actual
two-dimensional information, not merely bounded degree or local swapping.

## Exact square evidence and its precise limitation

The existing independent-set DP was recompiled and actually run for
\(P_n\square P_n\), \(2\le n\le6\); its exact root vectors are in
`n2-roots.json` through `n6-roots.json`.  The largest adjacent ratio occurs
on a corner--edge-neighbour and is

\[
1,\quad {11\over4},\quad {19681\over3560},\quad
{6786834357579\over719889959846},\quad
{68708495735539617742323275\over4793987217870606674220112}
\]

for \(n=2,3,4,5,6\), respectively; numerically these are
\(1,2.75,5.52837,9.42760,14.333\).  They grow on the observed scale of
\(n^2\), not a fixed constant.  This is finite exact evidence only.  It
does not prove an asymptotic lower bound and therefore does not refute a
square-specific \(\exp[O(n)]\) corner-root comparison.

At a corner \(o\) with the two symmetric inward neighbours, (2) gives
\(c_{o,u}=a_o/2\).  Therefore

\[
P(u,o)={a_o\over2a_u}. \tag{6}
\]

The observed large adjacent ratio is equivalently a small *inward* edge
transition.  Any viable proof of a sharp square root-transport estimate must
therefore prove a quantitative conductance bound for (3), rather than posit
an injective local label swap.  A possible scale compatible with these finite
data is a radial comparison of order \(N/k^2\) at distance \(k\) from a
corner: its product over \(k\le n\) has logarithm \(O(n)\), unlike applying
(4) independently at every step.  No such conductance estimate is proved
here.

This identifies a concrete first missing inequality: a lower bound on
\(c_{xy}/a_x\) for an outward/inward grid edge with its distance-to-corner
dependence.  The reversible identity (2) alone supplies no such bound.

## Coordinate unimodality audit

The exact vectors above are coordinate-nondecreasing from a corner toward
the central row or column for every one of the 6, 8, 20, and 24 oriented
edges available at \(n=3,4,5,6\), respectively.  This is compatible with,
but does not prove, square root-coordinate unimodality.  The implication
needed for a corner-root surface bound is stronger than a good staircase
comparison: one must also show that a centre point maximizes \(a_x\), or
otherwise control the total mass away from a chosen staircase.  D4 symmetry
only gives equality on D4 orbits and gives neither conclusion.

Two obvious local maps do not establish unimodality.  In the full 3-by-3
enumeration, interchange of the positions of a corner root 0 and its
neighbour 1 is a legal 1-rooted word for 672 of 1344 corner-rooted words and
fails for the other 672.  A concrete failure is

\[
(0,3,1,2,4,5,6,7,8)\longmapsto(1,3,0,2,4,5,6,7,8),
\]

whose second entry has no earlier neighbour.  The always-legal stable move
of 1 to the front has an 8-to-1 fibre on 60 outputs in the same enumeration.
Thus neither raw position transposition nor the stable-move construction is
an injective proof of coordinate monotonicity.

The next genuine statement to prove or disprove is therefore a **global**
conductance comparison, for example a distance-sensitive lower bound for
the inward transition in (5), together with a coordinate-unimodality or a
mass-domination theorem.  The present exact checks cannot supply either.

## A sufficient global transport lemma (and what it would buy)

Index the square by \(0\leq i,j<n\), with the selected corner at \((0,0)\).
The following is a precisely sufficient, but presently unproved, replacement
for coordinate unimodality.  Suppose a fixed \(C\geq1\), independent of
\(n\), satisfies on **every** horizontal grid edge

\[
 {a_{i+1,j}\over a_{i,j}}\leq {Cn^2\over(i+1)^2}
\qquad(0\leq i<n-1,\ 0\leq j<n). \tag{7}
\]

The transposed statement then holds on vertical edges.  Multiplying first
horizontally and then vertically gives, for every vertex,

\[
 {a_{i,j}\over a_{0,0}}
\leq { (Cn^2)^{i+j}\over(i!j!)^2}. \tag{8}
\]

Since \(n^k/k!\leq e^n\) for \(0\leq k\leq n\), (8) is at most
\((C^2e^4)^n\).  Thus, without any separate unimodality assertion,

\[
{A_1\over a_{0,0}}\leq n^2(C^2e^4)^n. \tag{9}
\]

This would establish the desired corner-root lower bound
\(a_{0,0}/A_1\geq\exp[-O(n)]\).  It also pinpoints why a comparison only
along one corner-to-centre staircase is insufficient: it says nothing about
the other \(n^2-O(n)\) root masses unless a separate mass-domination theorem
is supplied.

As an exact finite audit only, `verify_root_vectors.py` checks (7) with
\(C=1\) for all horizontal edges for \(2\leq n\leq6\).  The least observed
constants are respectively

\[
 {1\over4},\ {11\over36},\ {19681\over56960},\
 {6786834357579\over17997248996150},\
 {68708495735539617742323275\over172583539843341840271924032}.
\]

They are approximately \(0.25,0.3056,0.3455,0.3771,0.3981\), and each
maximum is at a corner--neighbour edge.  This confirms consistency of (7)
in the available exact range; it is not a proof for arbitrary \(n\), nor
does the reversible identity (2) presently imply it.

There is a sharper square-specific candidate.  For
\(i<\lfloor(n-1)/2\rfloor\), suppose
one could prove

\[
 {a_{i+1,j}\over a_{i,j}}
 \leq\left({n-1-i\over i+1}\right)^2. \tag{10}
\]

For even \(n\), the omitted central edge has equal endpoint counts by the
left--right reflection and its right-hand side in (10) would be exactly one.
Applying (10), this central equality when needed, and D4 symmetry afterwards
would give, for all \(i,j\),

\[
 {a_{i,j}\over a_{0,0}}
 \leq {n-1\choose i}^2{n-1\choose j}^2. \tag{11}
\]

Summing the two one-dimensional binomial-square bounds gives the especially
clean sufficient estimate

\[
 {A_1\over a_{0,0}}
 \leq\left(\sum_{i=0}^{n-1}{n-1\choose i}^2\right)^2
 ={2n-2\choose n-1}^2\leq16^{n-1}. \tag{12}
\]

The step from (10) to (11) uses symmetry only to reverse direction after
the centre; it does not falsely infer unimodality from D4.  Exact integer
data for \(2\leq n\leq6\) find no counterexample to (10) on its stated
strict-half range.  The largest observed ratio of its left-hand side to its
right-hand side is respectively

\[
0,\quad {11\over16},\quad {19681\over32040},\quad
{51539717273924\over61081509218211},\quad
{1640030468451435164898660033\over2129963367801728150012021525}.
\]

For \(n=2\) the strict-half range is empty; the single central edge was
separately checked to be an equality.  For \(n=3,4,5,6\), the displayed
ratios are \(0.6875,0.6143,0.8436,0.7700\ldots\); the central equality was
also checked in every even case.  This is still finite evidence, not a
theorem.  Unlike (7), it is deliberately stated only for squares;
the thin-rectangle interface obstruction gives no reason to expect an
unrestricted rectangular analogue.
