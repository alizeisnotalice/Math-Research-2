# Disjoint extensions of the size-twelve construction

The previous numerator uses four (6\times6) blocks with one local peak
in each. Its two paired orientations use a one-to-one assignment of
losing roots to neighboring winning blocks. The previous extra family
uses diagonal fully safe winners, with both off-center losing roots
facing the same winner. Use the same local connector rules here. All
ranks through seven are allowed. Later ranks are allowed when equal to
the first possible rank, one more than the Manhattan distance to the
local root. The symmetry of a block recovers the inner-column roots'
rank vectors by reflecting the outer-column vectors and the connector
coordinate. No new completion queries enter this calculation.

## Roots on the seam between pairs

Consider two winning quadrants in the same row half, with one winner
in each column half. Each has a vertically paired loser. In a paired
word, the winner's root precedes every losing label because the connector
precedes the losing root. This handles a winner's possible cross-edge
to its own losing block. Now allow one or both winner roots to occupy
the central vertical seam. If both do, their rows must differ so that
the roots are nonadjacent. Each such seam root has one other possible
external neighbor, in the opposite winning block, at its matching row.

If that neighbor is at distance (d) from the opposite winner root
within their block, its local rank is at least \(\ell=d+1\). Its rank
inside the opposite paired word is at least as large. Therefore it
suffices to place the first label of the seam root's paired word before
the \(\ell\)-th label of the opposite paired word. For two paired words
of length (2M) this deletes exactly \(\binom{4M-\ell}{2M}\) interleavings.
With two seam roots, the two rejected sets are disjoint since their
first letters are in opposite blocks. All other labels interleave
freely. In particular a single rejection fraction is at most
\(2^{-\ell}\), by the product formula for the probability of an initial
run of labels of one color. The transposed family covers the other
orientation. The new true peak pair is excluded from the old paired
family, since its sole possible pairing requires both winners to avoid
this seam. The previous dominant-winner family needs diagonal winners.

## Two losses along a chain

Again take adjacent winning blocks, and choose one as the first
winner. It suppresses the local root in its neighboring losing block
along the paired seam. That losing block has a connector on the seam
to the second losing block and suppresses the second losing root. The
second loser root lies off its central corner, on the boundary facing
the first loser. Thus it has no neighbor in the other winner's block.
The first winner can be on its own paired seam, since its root precedes
all first-loser labels, but it avoids the seam between winners. The
other winner avoids both seams. These root conditions ensure that the
two winners are the exact global peaks and make this family disjoint
from the seam-root extension. The second loser cannot be suppressed
by its adjacent winner, so the original one-to-one family does not
apply. Reversing which loser is suppressed second requires the former
second loser's root to lie on its upper seam, which is excluded.
The other axis and the prior diagonal family have different winning
quadrants. Both chain choices therefore add without multiplicity.

Let (i) be the connector rank of the first winner (W), and (j\) the
connector rank of the intermediate loser (L_1). With the other winner
set aside, the three-word interleaving must satisfy
\(W_i<L_{1,1}\) and \(L_{1,j}<L_{2,1}\) in descending-order arrival time.
If exactly (a\ge i) first-winner labels arrive before (L_{1,j}\), then
the first (i) labels must all come from (W). The next (a-i+j-1)
labels interleave (a-i) first-winner labels and (j-1) first-loser
labels arbitrarily, and the next label is (L_{1,j}\). The suffix is free.
For fixed local words the exact triple coefficient is
\[
 \sum_{a=i}^M
  \binom{a-i+j-1}{j-1}\binom{3M-a-j}{M-a,M-j,M}.
\]
Multiply by \(\binom{4M}{M}\) to restore the other winner freely, then
sum against the two connector rank frequencies and the two unconstrained
rooted counts. The chain is counted once per word, since the four local
one-peak words, the chain direction and the interleaving are recoverable.

## A diagonal central-corner winner

Take diagonal winners, one with its peak at its own central corner and
the other with a fully safe root. Require the two loser roots to lie at
their central corners. The first winner's root is adjacent to both loser
roots and is protected from any competing neighbor if it precedes both.
The other winner is safe by its root location. With that winner's word
set aside, the three-word interleaving starts with the corner winner,
and the remaining labels interleave freely. This gives factor
\(\binom{3M-1}{M-1,M,M}\binom{4M}{M}\). The true peak pair violates both old
paired root restrictions and the old fully safe diagonal root rule.
This case is also disjoint from the two adjacent-winner additions.

## Integer results and gates

For \(\rho_{12}>6\), the combined numerator is
`4325025309755870589020775687294797524902748176662048864785034015629981080105822436094012025998508121031992881237908475909963776093673839529416538884645121444881731082007556791884855329621633633414545408000000`.
It exceeds the previous numerator by a factor of about 1.024601. Using
the round-51 three-block denominator upper gives the lower bounds
\[
\begin{aligned}
 \mathbb P_2(\rho_{12}>4)&\ge 1.14094564347528760\times10^{-19},\\
 \mathbb P_2(\rho_{12}>6)&\ge 8.83683432413816308\times10^{-20},\\
 \mathbb P_2(\rho_{12}>7)&\ge 6.14053439391517232\times10^{-20}.
\end{aligned}
\]
Decimals are rounded downward. The ledger retains the exact partial
counts, exact ratios and true peak pairs. This calculation does not use
any new size-twelve denominator computation. Three-color interleavings
were exhaustively verified for word lengths one through four and paired
interleavings for lengths two, four and six. The old one-axis count was
recomputed on 3,600 adjacent-block peak pairs. The actual run took 0.18
wall seconds and 55,476,224 bytes RSS and made zero independent-set calls.

## Entropy of local prefix restrictions

A fixed shortest path of length (d) from a local root to a connector can
be moved to the start of any rooted one-peak word in path order, keeping
the order of all other vertices. All other vertices retain an earlier
neighbor, either outside the path in the same relative order or in the
promoted path. The original positions of the (d) nonroot path vertices
encode the inverse, with at most \((M-1)_d\) possibilities. The same
argument with the path in reverse order changes the root to its other
endpoint, preserving this bound. Thus a fixed proposed root or path
prefix costs at most a surface-order \(\exp[O(m\log m)]\) factor per
\(m\times m\) block. The label interleaving part is smaller. For equal
block sizes (M=m^2\), selecting the first (j\) labels from one of the
two blocks costs the factor \((M)_j/(2M)_j\). For (m\ge4\) and
(j\le2m-1\), the inequality \(\log(1-s)\ge-2s\) for \(s\le1/2\)
gives a lower bound \(2^{-j}e^{-4}\) for this factor. These local
bounds separate the surface logarithm from the area-order count of
single-peak block words.

Run `python3 -B extend_n12.py <fresh-output-json>` for the exact finite
calculation. It reads the bound and the rooted and completion counts
from the round-51 and round-50 archives and records their hashes.
