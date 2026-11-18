# Odd-quadrant gluing with exact short connected prefixes

Split the \(11\times11\) grid into four rectangles with heights 5 and 6
and widths 5 and 6. Every rectangle is required to have exactly one peak
in its induced labels. Its descending relative word therefore grows one
adjacent vertex at a time from that peak. The reference coordinates of
an upper or lower block count inward from the outside, so a block paired
vertically always has its connector on the reference bottom edge.

For the vertical pairing separately in each column half, require the winner
root to avoid the seam between the two column halves, and the loser root
to lie on the shared horizontal seam. The corresponding neighbor of the
losing root in the winning block must precede it globally. The winner root
precedes that connector or is the connector itself. Consequently it
precedes the entire losing word, so it is safe even if it lies on the
shared horizontal seam. It has no other cross-block edges. Every local
nonroot has a higher neighbor in its block, the two loser peaks are
suppressed, and precisely the two winning roots survive globally. The
horizontal direction uses the same rule after transposition.

## Count the allowed connector ranks without a full subset table

The \(5\times5\) block uses the exact full connector ranks and rooted counts
from round 48. Both 30-point rectangular orientations and the 36-point
square use the following rank criterion. A connector of local rank at
most seven is allowed. Beyond rank seven, it is allowed exactly at the
first possible rank, one more than the Manhattan distance from the root.
This criterion is invariant under reflections and transposition.

For local rank \(j\le7\), a connected prefix recursion counts all ordered
prefixes to rank \(j-1\), retaining only the set and its multiplicity in
that layer. If the next vertex is the desired connector, fix it next and
multiply the prefix multiplicity by the completion count from its occupied
set \(T\). At the minimum rank \(\mathrm{dist}(u,c)+1\), the prefix has just
enough vertices for a shortest path from the root \(u\) to the connector
\(c\). Any connected induced prefix of this size containing both endpoints
is precisely such a path. A shortest path has no chords, so its growth
order from \(u\) is the path order. Thus all minimum-rank frequencies beyond
seven come from enumerating the ordered shortest paths. Distinct initial
words are disjoint, even if their occupied sets agree.

For a fixed nonempty connected prefix set \(T\), put \(R=V\setminus T\) and
\(F=R\setminus N\(T\)\). Vertices of \(R\setminus F\) already have a higher
neighbor in \(T\). A remaining word is a valid completion exactly when
there are no peaks of its induced labels on (R) in (F). To count it,
sum by inclusion–exclusion over independent subsets \(I\subset F\), using
\[
 C_R(\varnothing)=|R|!,\quad
 C_R\(I\)=\frac{\sum_{v\in I}C_R(I\setminus\{v\})}{|I\cup N_R\(I\)|},\quad
 h\(T\)=\sum_I(-1)^{|I|}C_R\(I\).
\]
The moment recurrence follows by conditioning on the greatest remaining
label in the union of the relevant closed neighborhoods. For independent
\(I\), that label's vertex has no effect on the other peak conditions.
The completion count depends on \(T\) alone. Separate actual prefixes
can share a symmetry-reduced completion query, and transposed rectangles
use the same query table. The exact signed sums fit U192 since each
positive or negative sum is at most \(2^{36}36!\). The program checks
integer division and nonnegativity on every query. It enumerates only
independent sets in (F), with no \(2^{36}\) table.

If the two paired blocks have sizes \(M,L\), and the chosen connector has
local rank \(j\), the first \(j\) labels of the merged pair must come from
the winner. All subsequent labels can interleave arbitrarily. The exact
factor is therefore \(\binom{M+L-j}{L}\). Sum this against the winner's
rank count and the corresponding loser seam-root count. The two paired
halves of sizes 55 and 66 are interleaved in \(\binom{121}{55}\) ways.
Given a word in a fixed pairing, the true peak pair determines the winning
blocks, and the four relative words determine the losing local roots and
connector ranks. This is a one-to-one count.

## Exact overlap of the two pairings

In any word belonging to both axes, the winners lie in diagonal blocks
and avoid both central seams. The winning connectors and losing roots
are the central corners in their own blocks. Each winner corner must
precede each loser root. For a fixed winner corner-rank pair \(i,j\),
condition on the first losing label arriving after exactly \(a\ge i\) labels
of the first winner and \(b\ge j\) of the second. There are
\(\binom{a+b}{a}\) prefix interleavings. If the first loser belongs to the
first losing block, the suffix has lengths
\((M_1-a,M_2-b,L_1-1,L_2)\); if it belongs to the other, the lengths are
\((M_1-a,M_2-b,L_1,L_2-1)\). Sum their multinomial coefficients and then sum
over \(a,b\). The rank conditions on the central crosses agree between the
orientations, because their restrictions are invariant under transpose.
The implementation checks equality of those rank vectors for every
diagonal root pair. This overlap is subtracted within each unordered true
peak pair, and distinct unordered peak pairs are disjoint.

## One winner suppresses both losers

When the two winners are diagonal and both roots avoid the two seams,
choose one dominant winning block. Prescribe that the losing roots face
that block and avoid their own central corners. Each losing root can
have only one cross-block neighbor, which is in the dominant block.
Require both such neighbors to precede both losing roots globally. Then
all nonroots and both loser roots are suppressed, while both winner roots
remain. The other winning block cannot suppress either local loser root.
Both one-to-one pairings are therefore impossible. The two choices of
dominant block are disjoint because the seam occupied by either off-center
loser root determines its only neighboring winner block.

Delete the nondominant winner word temporarily. If \(j\) is the later rank
of the two dominant connectors, the three remaining words have an initial
\(j\) labels from the dominant block, followed by arbitrary interleavings of
lengths \((M-j,L_1,L_2)\). Restore the nondominant winner in any of its
(\binom{121}{M_{\rm other}}\) placements. For the 30- and 36-point blocks
restrict this later connector rank to at most seven. At the moment the
later connector is added, its completed set has already been queried in
the single-connector table, either directly or by transposition. The
program reuses those completion values and enumerates no new independent
sets. For the \(5\times5\) block use the complete double-connector ranks
from round 48. The additional family is counted separately by the true
peak pair and added to the axes union.

## Exact finite comparison

The local run certified 4,150 nontrivial completion contexts, with
37,598,697 independent sets processed and no context pruned. The maximum
worker RSS was 79,986,688 bytes, and the batch was waited in 4.11 seconds
wall. The \(3\times3\) short ranks were independently checked against all
\(9!\) words and the connected-suffix recurrence. A separate \(5\times5\)
check used 220 contexts to reproduce 450 bins in the older full-rank
ledger. Unequal-size interleaving equations were checked by exhaustive
small position lists. The added two-connector weights were checked
against round 48's independent \(3\times3\) double-rank table.

For the original Manhattan threshold \(\rho_{11}>11/2\), the exact
axes union numerator is
`50276944924755328697099177844258496482550658912593911356132675784937438089454318281967880441588008587020171862259818706326586964154101184654650909545922764732661760000`.
The disjoint dominant family brings the exact numerator to
`52351689802293466136392599474327417763846016283777900331091904619007483803562117035474957335425461404117113295762653347037690361691229230041429695613248835287285760000`.
Dividing by the actual round-49 dense integer upper \(A_2(11)\) bound gives
\[
 \mathbb{P}_2\(\rho_{11}>11/2\)\ge1.39276530389707903\cdot10^{-12}.
\]
This is about 67.43 times the previous aligned-root lower bound at the
same threshold. The corresponding lower bounds for \(\rho_{11}>11/3\) and
\(\rho_{11}>6\) are \(1.61560014245426908\cdot10^{-12}\) and
\(1.04333902403860655\cdot10^{-12}\). The displayed decimals are rounded
downward. These are finite bounds and do not settle the infinite-subsequence
question.

## Replay

Build with Apple `clang++ -std=c++20 -O3 -Wall -Wextra -Werror`, using the
project's data-volume module cache, and name the binary `prefix_completion`.
Python 3.11 or later with the standard library is sufficient. The replay
uses the archived rooted counts, full \(5\times5\) rank tables and protected
bound at the paths in `merge_n11.py`, and the existing local plan files.
For new output directories and files run:

```sh
python3 -B run_prefixes.py plan-small.json small-results --workers 1
python3 -B run_prefixes.py plan-n5-gate.json results-n5-gate --workers 1
python3 -B audit_n5.py
python3 -B run_prefixes.py plan-K7.json results-K7 --workers 3
python3 -B merge_n11.py plan-K7.json results-K7/summary.json n11-gluing-K7.json
python3 -B add_dominant.py n11-gluing-K7.json results-K7/summary.json n11-dominant-K7.json
```

`build_prefixes.py plan-K7.json --rect-cutoff 7 --square-cutoff 7` regenerates
the main plan. Its `--small-gate` mode regenerates the \(3\times3\) plan.
The result ledgers save exact integers and fractions, input hashes, peak
pairs, original distances, and the explicit inclusion–exclusion terms.
