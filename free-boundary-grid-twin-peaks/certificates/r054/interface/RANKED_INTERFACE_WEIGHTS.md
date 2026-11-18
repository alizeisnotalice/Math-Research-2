# Certified interleaving weights from boundary ranks

Fix a partition of the free-boundary grid into blocks. A full descending
word is uniquely specified by its relative word in each block and a color
word interleaving those block words. Under the original uniform law the
relative words are independently uniform, and conditional on them the
color word is uniform among all interleavings. This identity preserves
the original normalizer even when conditioning on the final peak count.
The rank of a vertex within a block below always starts at one.

## Exact fixed-word kernel

A vertex that is not a local peak in its block cannot be a global peak.
When a local peak is revealed, it becomes a global peak exactly if no
cross-block neighbor has been revealed. Its future neighbors are all
lower. Thus for fixed relative words, use the lattice point of counts
already revealed in each block. A transition adds a surviving peak if
its next vertex was a local peak and the counts of the other blocks
are all below the local ranks of its cross-block neighbors. A dynamic
program on this lattice counts color paths by global peak set. It can
discard any state with more than two recorded peaks when the goal is
exactly two, since every recorded peak remains a peak. The probability
weight of a prescribed peak pair is its path count divided by the full
multinomial interleaving count.

For two rectangles joined by a straight interface, every seam endpoint
has one external neighbor. A seam pair with two local peaks contributes
exactly one full-grid peak, the earlier endpoint. A seam pair with one
local peak contributes a peak exactly when that endpoint precedes the
other. This explains both the count lower bound obtained from the joint
peak mask and the additional dependence on relative ranks.

## A constant interleaving weight from bilateral rank margins

Take (r\ge2\) blocks of equal size (M), with fixed relative words. Suppose
the entire set of desired surviving peaks consists of two protected
local peaks. For each other local peak of rank (i), designate at least
one cross-block neighbor of local rank (j\le i-\Delta\), with the same
integer margin (\Delta\) for every designated edge. The neighbor need
not be a local peak. Require along the color word that all pairs of block
counts differ by at most (\Delta\). At the transition revealing the
undesired peak, the count of the suppressor block is at least
(i-\Delta\ge j\), so the suppressor has already arrived. Both protected
peaks remain. This works for conflicting suppression directions provided
the rank margins hold individually.

A specified two-color projection of the uniform (r)-color shuffle is
uniform on the \(\binom{2M}{M}\) binary bridges. With (X=a-b) denoting
their count difference, reflecting after the first hit of height (d>0)
gives the exact identity
\[
 \Pr(\max X\ge d)=\frac{\binom{2M}{M+d}}{\binom{2M}{M}}.
\]
The reflected paths end at (2d), and symmetry gives the same count
for negative hits. Thus the desired complete interleaving fraction is
at least
\[
 \gamma_{r,M}(\Delta)=
 1-2\binom r2\frac{\binom{2M}{M+\Delta+1}}{\binom{2M}{M}}.
\]
If the displayed number is negative, the trivial lower bound zero is
used. For (d\le M), the ratio is bounded above by
\[
 \prod_{q=1}^{d}\frac{M-q+1}{M+q}
 \le \exp\left(-\frac{d^2}{M+d}\right),
\]
by applying \(\log(1-u)\le-u\) and summing (2q-1\). If (d>M)
the ratio is zero. In particular, for
\(\Delta\ge\lceil\sqrt{2M\log(4\binom r2)}\rceil\) this certifies
interleaving weight at least (1/2\). For a fixed number of macroscopic
grid regions the rank margin is of order the grid side length.

The same reflection count bounds the cost of an adverse mandatory edge.
If a particular unwanted local peak of rank (i) must receive its only
possible higher cross-neighbor of rank (j>i), then its elimination
requires a two-color excursion to height at least (j-i\) in the
neighbor's direction. The upper fraction on all such accepted
interleavings is at most \(\binom{2M}{M+j-i}/\binom{2M}{M}\).

## Polynomial interleaving loss with directed block order

Suppose all undesired local peaks can instead use designated cross
neighbors in earlier blocks of one fixed total block order, each with
neighbor rank (j\le i). Impose (a_1\ge a_2\ge\cdots\ge a_r\) at every
prefix. This guarantees the needed earlier arrival at the step for any
such local peak. Its exact fraction among all color words is
\[
 \theta_r(M)=(M!)^r\prod_{i=0}^{r-1}\frac{i!}{(M+i)!}.
\]
To count these paths, translate the coordinates by (r-i), so the weak
order on the (a_i) becomes a strict noncollision condition. Take
alternating counts of paths starting at every permutation of the shifted
initial vector and ending at the same shifted terminal vector. A
path with a collision cancels with the path obtained by swapping the
two coordinate histories before the first collision. The uncancelled
paths are exactly the ordered ones. Their count is
\((rM)!\det[1/(M-i+j)!]_{i,j=1}^r\), with the negative factorial entries
taken as zero. Writing each entry as a falling factorial in (M+r-i)
divided by \((M+r-i)!\) gives a Vandermonde determinant. Division by
the full multinomial proves the formula. For (r=4), the fraction is
\(12/[(M+1)^3(M+2)^2(M+3)]\). This permits systematic chain or other
acyclic suppression patterns at polynomial, rather than exponential,
interleaving cost under the local rank hypothesis.

## Same seam-mask data with exponentially different weights

For every even (n=2w\ge6), take the equal vertical halves and a fixed
row. Put the left protected target peak (p) in its outermost column
and a second left local peak (u) at that row's seam endpoint. Start
the left relative word with (u,p). They are nonadjacent and the rest
can grow adjacent to the revealed set, giving exactly those two local
peaks. On the right, start at the protected vertex (q) next to the seam
endpoint (v) across from (u). There are two allowed one-peak words:
reveal (v) immediately, or postpone it until the end while growing
in the connected region with (v) removed. Taking the chosen row at
the outer horizontal edge suffices to ensure this deletion is connected.
Both words have the same local peak set, with the same seam-peak mask.
The protected peaks (p,q) have original distance (w+1\).

In either example the exact original two-peak event at this fixed
word pair requires (v) earlier than (u). Since (u) has local rank
one, with (v) at rank (j) the complete interleaving fraction is
\(\binom{2M-j}{M}/\binom{2M}{M}\), where (M=nw\). For (j=2) it tends
to (1/4\), and for (j=M\) it is \(1/\binom{2M}{M}\). The central
binomial coefficient lies between (4^M/(2M+1)\) and (4^M\). Thus
identical local masks and identical protected roots can have
interleaving weights differing by \(\exp(\Theta(M))\). This rules out
using the seam mask alone to bound the fixed-word weight from below
at a subarea exponential scale. It leaves open an average comparison
inside a mask cell.

At size twelve, these two exact fractions are (71/286\) and
(1/1480212998448786189993816895482588794876100\). A separate size-eight
example has two unwanted seam peaks, both at local rank 30, with the
two corresponding opposite suppressors at rank 2. Both directions are
active. The protected roots have Manhattan distance 10. The width-12
corridor certifies the accepted peak pair and its whole exact interleaving
fraction is (3966718919789086/3966718919789157\).

## What this says under the original normalization

For a fixed partition write \(\mu\) for the product uniform law on
relative words, and let \(Z_2=A_2(n)/(n^2)!\). Let \(\mathcal L_{\Delta,\varepsilon}\)
be the word class satisfying the bilateral margin hypotheses, with its
two protected local peaks separated by more than \(\varepsilon n\).
Then the original conditional probability satisfies
\[
 \mathbb P_2(\rho_n>\varepsilon n)
 \ge
 \gamma_{r,M}(\Delta)\frac{\mu(\mathcal L_{\Delta,\varepsilon})}{Z_2}.
\]
The fixed-order version replaces the margin class by its ordered
rank class and replaces \(\gamma\) by \(\theta_r(M)\). Both formulas
use the true global normalizer. A fixed positive mass on an infinite
subsequence would follow from a matching lower estimate for the
respective relative-word class. The proof certifies the interleaving
part and supplies no such mass estimate. The same-mask example explains
why the exact joint-mask upper counts alone cannot give this estimate
pointwise.

## Replay

Run `python3 -B interface_weights.py <fresh-output-json>` on the data
volume. Six sampled pairs of relative words on the (4\times4) grid
were each tested over every interleaving, with 77,220 literal peak-set
checks against the lattice dynamic program. The sample only selects
word pairs, so each selected interleaving identity is tested exhaustively.
The reflection identity was checked on all bridges of side length one
through six. The ordered and corridor formulas were checked on small
multicolor words, including all 369,600 color words of length twelve
with four equal color classes. The example constructions and their
exact fractions were checked independently in the same run. Actual
wall time was 2.20 seconds and maximum RSS 27,672,576 bytes. The local
process exited and the exact tables and resource receipt are saved.
