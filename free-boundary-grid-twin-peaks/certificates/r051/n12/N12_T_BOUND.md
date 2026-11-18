# A size-twelve upper count from three protected rectangles

Divide the (12\times12) grid into its (12\times4) top strip and the two
(6\times8) rectangles below it. In each block protect the vertices whose
original graph neighbors are all in that block. The top block protects
(12\cdot3=36) vertices. Each bottom block protects (5\cdot7=35) vertices,
after excluding its top edge and the edge facing the other bottom block.
The total is 106, compared to 100 for the four equal quadrants. The peak
test at a protected vertex is the same in the induced block and the full
grid. All labels are the distinct ranks of the original grid.

Let (Z) count the protected peaks across the blocks. The event (K=2)
requires (Z\le2). When (Z=2), a peak at any unprotected vertex is an
additional third peak. The first two label positions give the same
necessary filters as in the round-47 method. The global maximum must be
protected. If the globally second largest vertex is unprotected, it must
be adjacent to the global maximum. A protected maximum and a vertex in
a different block cannot be adjacent.

A further filter applies when the first two global labels are protected
and belong to distinct blocks. If the third label is unprotected, it must
have an earlier neighbor to keep it from being an additional peak. In a
third block that is impossible because both earlier vertices are protected.
If it lies in either of the first two blocks, it is the second local label
of its own block and needs adjacency to that block's first label. The
other earlier label is protected in a different block and has no cross-edge.
This third-label filter removes only words with first two labels protected,
so its exclusions are disjoint from both earlier exclusions.

## Local histograms and exact interleaving

For a rectangular block (H), let (P) and (B) be its protected and
unprotected sets. Given exactly (i) local peaks among (P), write (F_i)
for the total number of relative words, (T_i) for the count with largest
label in (B), and (V_i) for the count with largest label in (P) and
second label in (B) nonadjacent to that largest. The three histograms
are the same quantities as in the earlier top-two method. They depend
only on the induced rectangular graph and (P).

For independent (S\subset P) let (W=S\cup N_H(S)), and let (C(S))
be the count of words in which every vertex of (S) peaks. The exact
moment recursion is (C(\varnothing)=|H|!\) and
\[
 C(S)=\frac{\sum_{v\in S}C(S\setminus\{v\})}{|W|}.
\]
The corresponding restricted moments for (T,V) use the factors
\[
 \frac{|B\setminus W|}{|H|},\qquad
 \frac{(|P\setminus W|+|W|)|B\setminus W|
             -e_H(P\setminus W,B\setminus W)}{|H|(|H|-1)},
\]
respectively, multiplying (C(S)). These factors come from placing the
top one or two labels as in round 47. For the second factor, the top
protected vertex either lies outside (W), requiring the top pair to be
chosen in the two outside sets and not adjacent, or lies in (S). In the
latter case the top unprotected nonneighbor is outside (W), and summing
the predecessor moments over (S) gives (|W|C(S)). Aggregate moments
by cardinality of (S) and invert by the binomial transform. The
implementation uses the ranked-layer method from round 49 and separate
masks for the two parts of a closed neighborhood. Its four certified
moduli have a product larger than (48!\), enough for every counted
block histogram entry. The fifth certificate node tree is also verified
by the driver but is not needed for these output counts.

For fixed local protected peak counts summing to at most two, interleave
all relative words by the usual multinomial. When those counts sum to two,
apply the first two exclusions with the exact top-two formula. Then subtract
the new disjoint third-label failures as follows. Denote by (i,j) the
ordered distinct block indices of the first two labels, and set
(S_i=F_i-T_i) at the prescribed local peak count. The case with third-label
block (i) contributes (V_i S_j), the case with block (j) contributes
(S_i V_j), and the case with another block (k) contributes (S_i S_j T_k).
Multiply each term by the unrestricted (F\) counts of the other blocks
and by the multinomial of the residual word after fixing these first three
block indices. Sum over permitted indices and local peak counts. The
subtractions count words failing the stated necessary conditions exactly.

## Measured size-twelve bounds

| Partition | Protected points | Top-two upper on \(A_2(12)/144!\) | Added third-label filter |
|---|---:|---:|---:|
| Four equal quadrants | 100 | \(3.986515204\times10^{-22}\) | \(2.508433872\times10^{-22}\) |
| Three-block T | 106 | \(1.330321645\times10^{-23}\) | \(8.818119201\times10^{-24}\) |

These upper decimals are rounded upward; exact integers and fractions
are available from the ledger. The T partition improves the final
four-quadrant denominator bound by a factor of 28.446. The two 48-point
T-side types process 6,879,979 and 3,729,091 protected independent sets,
with at most 22,962,896 and 13,017,336 bytes of live residue table
respectively per modulus. The final run has 32 completed and waited local
modular tasks, maximum actual worker RSS 36,716,544 bytes and 5.80 seconds
wall. The smaller scientific calibration included full relative-word
enumeration on two small rectangles, the exact round-49 size-eleven
quadrant histograms and top-two upper count, and all \(8!\) words in an
abstract three-block shuffle gate that exercises the third-label terms.
All these checks passed.

Two more ambitious decompositions were planned only. A straight
\(12\times6\) half has 60 protected points and requires about 352.676 GiB
of dense 64-bit residue tables per modulus. The T partition with a
height-five upper block has 108 protected vertices in total, but its
upper block alone requires about 2.733 GiB of live 64-bit tables per
modulus. These exceed the resource limit of this round.

## A numerator using the existing size-six completion counts

Partition the same grid instead into four \(6\times6\) quadrants and
require one local peak per block. Use the round-50 rank constraints
on the winning connectors, with all growth prefixes of rank at most
seven and the shortest-path rank beyond seven. Either paired winner
may be on its own paired seam provided it avoids the seam between
pairs. Each root precedes its connector, which precedes the entire
paired loser word, so its cross-edge neighbors are suppressed.
The exact interleaving count is \(\binom{72-j}{36}\) for a fixed local
connector rank \(j\), followed by \(\binom{144}{72}\) interleavings of
the two pair words. Compute the vertical and horizontal families
and subtract their overlap at each unordered full-grid peak pair.
The overlap requires diagonal winners and both winner central
corners to precede both losing central roots. The first-loser position
multinomial formula of round 50 evaluates it exactly.

An additional family has diagonal safe winners, both losing roots on
seams facing one dominant winner and both off their central corners.
The two corresponding points of that winner precede both losing roots.
This is disjoint from the two one-to-one pairings because neither losing
root has any cross-edge neighbor in the nondominant winning block. It is
also disjoint across the two dominant choices. The dominant connector
pair's later rank is restricted to seven, and its conditional completions
are already included in the single-connector prefix ledger, up to
transposition. No new independent-set queries were necessary. The two
pairs and dominant class are counted exactly under their stated prefix
restrictions. The final ledger resolves every one of 5,950 supported
unordered peak pairs in the original Manhattan metric.

The exact count for \(\rho_{12}>6\), divided by the T denominator upper,
gives
\[
 \mathbb P_2(\rho_{12}>6)\ge8.62465939651556539\times10^{-20}.
\]
The corresponding bounds for \(\rho_{12}>4\) and \(\rho_{12}>7\) are
\(1.11204963341321861\times10^{-19}\) and
\(5.99885379184909944\times10^{-20}\), with all lower decimals rounded
downward. At threshold six the numerator is 38.32 times the aligned-root
subfamily. Using the same numerator and the four-quadrant denominator
would instead give about \(3.0319\times10^{-21}\). At size eleven, the
previous conditional lower is about \(1.3928\times10^{-12}\), considerably
larger. The size-twelve T normalization controls only some necessary
conditions, and its upper count is not the true \(A_2(12)\). These finite
bounds do not determine the behavior as \(n\) increases.

## Replay

Compile `rect_protected.cpp` with Apple clang using `-std=c++20 -O3
-Wall -Wextra -Werror` to `rect_protected`. The driver uses Python 3.11+
and macOS `/usr/bin/time -l` for resource receipts. Supply
`pratt-61bit.json` in the same directory. Run the following with fresh
output paths:

```sh
python3 -B count_rectangles.py denominator-final --workers 2
python3 -B numerator_n12.py denominator-final/summary.json n12-conditional-final.json
```

The denominator calibration reads the round-49 four-block ledger. The
numerator reads the exact rooted counts for size six and the round-50
completion and construction utilities from the project archive. These
input hashes, exact regional histograms, resource receipts, exact
conditional fractions, and distance-resolved numerator counts are
recorded in the result files. All local computation processes were
checked to have exited.
