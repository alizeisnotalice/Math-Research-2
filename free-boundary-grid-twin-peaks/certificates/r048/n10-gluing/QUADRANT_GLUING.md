# Exact one-peak quadrant gluing and overlap of the two axes

Divide the free-boundary \(2m\times2m\) grid into its four \(m\times m\)
quadrants. A quadrant word lists its vertices in decreasing label order. A
word has exactly one peak in the induced quadrant if and only if every
vertex after the first is adjacent to an earlier vertex in that quadrant.
Consequently its occupied prefixes are connected, and the first vertex is
its unique peak. Throughout the construction all four restrictions have
this property.

In the vertical pairing, pair the upper and lower blocks in each of the
two column halves. Select one winner in each pair. In a reference quadrant,
measure \(x\) inward from the outside edge of the full grid and \(y\) inward
from its top or bottom outside edge. The winner root must have \(x<m-1\).
The loser root must lie on the shared horizontal seam, and its corresponding
neighbor in the winner must precede it in the whole word. Then every
nonroot has a larger same-quadrant neighbor and each loser peak has a larger
neighbor across its pair seam. The winner root has no neighbors in the
other column half. Any of its neighbors outside its own quadrant are in the
assigned losing quadrant, and all losing labels are below its root because
they are below the winning connector. Hence exactly the two winner roots
are global peaks, even if a winner root lies on its own pair seam. The fully
safe subclass further requires \(y<m-1\) for both winner roots. The
horizontal pairing is the transpose of this definition.

## Exact ranks and unrestricted interleaving counts

Let \(M=m^2\), let \(R_u\) be the number of quadrant words rooted at \(u\), and
let \(W_u(c,j)\) count such words with boundary connector \(c\) at rank
\(j\) in the winner. After the relative orders of a winning and a losing
quadrant are fixed, exactly \(\binom{2M-j}{M}\) interleavings place the
connector before the first losing label. In such a pair word the first
\(j\) labels must come from the winner, and every subsequent interleaving
is allowed. Thus the pair count for a fixed winner root \(u\) is
\[
 A_u=\sum_{c,j} W_u(c,j)R_{\ell(c)}\binom{2M-j}{M},
\]
where \(\ell(c)\) is the boundary point in the loser directly across from
\(c\). Two such pair words admit arbitrary \(\binom{4M}{2M}\) interleavings.
This counts every valid word in the stated family, including all allowed
ranks of the connector. The full word determines the orientation's winner
blocks from its two global roots, the unique local loser roots, the local
orders, and both interleavings. Hence a single orientation has no duplicate
words from the choices involved in the construction.

Compute \(W_u(c,j)\) by subset DP. For a set \(S\) containing \(u\), let
\(g_u(S)\) be the number of connected growth prefixes on \(S\) starting at
\(u\). For any connected nonempty set \(T\), let \(h(T)\) count the orders of
the remaining vertices that continue connected growth. The forward update
adds an available frontier point, and the backward recurrence sums over
the same choices with terminal value one for the entire quadrant. Then
\[
 W_u(c,j) = \sum_{\substack{|S|=j-1,\ c\notin S,\ c\in N(S)}}
                  g_u(S)h(S\cup\{c\}),
\]
for \(c\ne u\), with \(W_u(u,1)=R_u\) if the root is itself a connector.
Every connected suffix count is at most the word count of any root in the
prefix, by attaching one fixed valid order of its prefix. The C++ arrays
store 64-bit counts, calculate the suffix additions at 128-bit width and
use a saturation marker on unused masks. It checks that no single-root or
used connected completion count is saturated, and that each rank histogram
sums to the correct rooted count. The four observed saturated suffix masks
at \(m=5\) do not enter any counted connected prefix.

## Exact overlap of the axes

A word in both orientations has one winner root in each column half and
one in each row half, so the winning blocks are diagonal. Each root avoids
both central seams by the two orientation rules. Each local losing root
must lie on both of its relevant seam boundaries, so it must be its block's
central corner. For a winning block the neighbors across its two seams to
those losing roots are the same vertex, its central corner. The two
orientations are therefore satisfied exactly when the two winner corners
both precede both losing roots, for these fixed local root assignments.

Let \(i,j\) be the ranks of those two winning corners. If the first losing
label occurs after exactly \(a\) labels from the first winner and \(b\) from
the second, then \(a\ge i,b\ge j\). There are \(\binom{a+b}{a}\) orders of
that prefix, two choices of the first losing block and a freely interleaved
suffix. Thus its exact interleaving coefficient is
\[
 \sum_{a=i}^M\sum_{b=j}^M
 2\binom{a+b}{a}
 \binom{4M-a-b-1}{M-a,M-b,M-1,M}.
\]
A two-dimensional suffix sum obtains all coefficients. Multiply by the two
rank frequencies of the winning corners and by the squared rooted count of
a losing corner. This gives the intersection for each diagonal unordered
peak pair. On that pair, add the two orientation counts and subtract the
intersection. Different unordered global peak pairs are disjoint.

## Integer evaluation for \(m=5\)

The full selector permits a winner on its own pair seam, but excludes the
other seam. It contains the earlier aligned-local-root construction, for
which the winning connector is the winner's root itself and \(j=1\). For
that subcase the merger separately reproduces the earlier exact numerator
by using \(2\binom{2M-1}{M}=\binom{2M}{M}\). This earlier numerator is not
added to the enlarged count.

Let \(U\) denote the separate exact upper bound on \(A_2(10)\) from the
protected-peak `5555555555` count. The finite lower bounds from this
construction, obtained by dividing the exact union numerators by \(U\), are

| Family | \(\rho_{10}>10/3\) | \(\rho_{10}>5\) | \(\rho_{10}>6\) |
|---|---:|---:|---:|
| Fully safe winners | \(5.50296\times10^{-11}\) | \(4.18987\times10^{-11}\) | \(3.27370\times10^{-11}\) |
| Winners may lie on their pair seam | \(2.31917\times10^{-10}\) | \(1.76375\times10^{-10}\) | \(1.15342\times10^{-10}\) |

The 5-threshold enlarged numerator is
`75870104005340810432537386962117955628138810535172888351553029757904140709909413534534308813098281211758536525070764036761600000000`.
The lower bound is about 249.09 times the earlier aligned-root value. The
union ledger records all 2,688 supported unordered peak pairs with exact
integer counts and cuts the true original Manhattan distance per pair.
No claim is made for an infinite sequence.

Ranks were checked against all permutations on the \(2\times2\) and
\(3\times3\) quadrant, and the overlap interleavings against direct position
enumeration for local word lengths one, two and three. The \(5\times5\)
rooted totals agree coordinatewise with the earlier independent one-peak
inclusion–exclusion ledger. The two m5 versions agree on all shared roots
and ranks. The enlarged m5 run took 1.65 seconds wall with peak RSS
404,504,576 bytes. Its two dense tables require 384 MiB of payload. The
local programs terminated after producing the output files.

## Replay

Use Apple `clang++ -std=c++20 -O3 -Wall -Wextra -Werror quadrant_ranks.cpp
-o quadrant_ranks`. In a fresh data-volume directory run the counter with
arguments `2`, `3`, and `5` to create `ranks-m2.json`, `ranks-m3.json`,
and `ranks-m5.json`. The saved `ranks-m5-safe-only.json` is the earlier
shared-rank gate. Run `python3 -B merge_quadrants.py <fresh-output-json>
<incumbent-json> <A2-upper-summary-json>`. The merger checks the exact input
hash tying the incumbent to the upper count. The files used here are the
round-48 `n10-quadrants.json` and the separate `cut-5555555555/summary.json`.

## An additional disjoint family with one winner suppressing both losers

Put the two fully safe winning roots in diagonal quadrants. Choose one
of those quadrants as the dominant winner. In each losing quadrant,
require its local root to lie on the seam facing that dominant quadrant
and to avoid its central corner. The root then has only one neighbor
outside its own block, the corresponding point of the dominant winner.
Require both such corresponding points to precede both losing roots.
All four subwords still have one peak, so both losing peaks are suppressed
and both safe winning peaks persist. The exact pair of global peaks again
specifies the winning blocks.

The other winning block cannot suppress either local loser peak. Therefore
none of these words has the horizontal or vertical one-to-one suppression
pattern already counted above. The two choices of dominant winner are
also disjoint, since the facing seam of either off-center loser root
specifies its only possible neighboring winner block. Both disjointness
claims hold pointwise, before counting and independently of interleavings.

Reflect the dominant block to the upper-left reference quadrant. Its
bottom connectors are \((a,m-1)\) and its right connectors are \((m-1,b)\),
with \(a,b<m-1\). Let \(W_u(a,b,j)\) count rooted words in which the later
of these two connectors is at rank \(j\). Use the same prefix and suffix
growth recurrence as above. When adding either connector to the frontier,
check whether the other connector already belongs to the prefix. For
that pair, multiply prefix count by suffix count and add to the later
rank bin. Every valid word contributes exactly once at that later rank.
Each pair of connectors has a histogram summing to \(R_u\), and transposing
root and both connectors gives the corresponding transposed histogram.

After deleting the second, nondominant winner's word, the three remaining
words must start with \(j\) dominant labels, followed by an arbitrary
interleaving of lengths \(M-j,M,M\). The exact factor for the three words
is \(\binom{3M-j}{M-j,M,M}\). Then interleave the nondominant winner's word
in \(\binom{4M}{M}\) ways. Thus with
\[
 D_u=\sum_{a,b,j}W_u(a,b,j)R_aR_b\binom{3M-j}{M-j,M,M},
\]
a safe diagonal root pair \(u,v\) contributes exactly
\(\binom{4M}{M}(D_uR_v+R_uD_v)\) additional words after reflecting both roots
into their quadrant coordinates. Here \(R_a,R_b\) denote the relevant losing
boundary-root counts, which equal the indicated row-zero rooted counts by
square symmetry. All words are retained with unrestricted placement of
the nondominant winner's labels and no duplicate among the choices above.

The (m=5) additional family contains 512 unordered global peak strata.
At threshold \(\rho_{10}>5\) it adds the exact count
`4852481754626542337493562946825994196093399592475050083010902637799170608342836810255386107197040339549905149297421582336000000000`,
giving a combined count of
`80722585759967352770030949908943949824232210127647938434563932395703311318252250344789694920295321551308441674368185619097600000000`.
The same protected-peak denominator therefore gives finite bounds
\[
\begin{aligned}
 \mathbb P_2(\rho_{10}>10/3)&\ge2.43197899584563298\times10^{-10},\\
 \mathbb P_2(\rho_{10}>5)&\ge1.87656241514237226\times10^{-10},\\
 \mathbb P_2(\rho_{10}>6)&\ge1.26335227102482102\times10^{-10}.
\end{aligned}
\]
All decimals here are rounded downward. For threshold five the combined
bound is about 265.02 times the original aligned-root bound. The axes-only
results remain in `gluing-summary.json` and the added and combined results
are in `dominant-summary.json`.

Build `double_connectors.cpp` by the same clang command, with output binary
`double_connectors`, and run inputs `2`, `3` and `5` into `double-m2.json`,
`double-m3.json` and `double-m5.json`. The (m=5) run took 1.27 wall seconds
at 404,570,112 bytes maximum RSS. Run `python3 -B add_dominant.py
<axis-union-json> <fresh-output-json>`. It independently checks every
pair-and-rank bin against all permutations at (m=2,3), and the three-word
interleaving factor by direct position enumeration at word lengths one
through four. All computation processes were checked to have exited.
