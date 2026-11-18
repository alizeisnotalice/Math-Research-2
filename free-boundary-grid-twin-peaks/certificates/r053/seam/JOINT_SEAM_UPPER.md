# A joint seam-peak upper count in the three-rectangle partition

Use the height-five T partition of the (12\times12) grid. Its top block
has area 60, and its two lower blocks are (6\times7) rectangles. A
protected point has no full-grid neighbor outside its block. Let (Z)
count peaks on all protected points. These are genuine full-grid peaks.

On the common vertical seam of the two bottom blocks, mark each pair of
adjacent vertices at the same row except in the topmost row of those
blocks. The six remaining pairs have the following property. Each
endpoint has exactly one neighbor outside its own block, the opposite
endpoint. If both endpoints are local peaks in their own induced
blocks, the endpoint with larger label also peaks in the full grid.
Distinct pairs produce distinct global peaks outside the protected set.
Bit (j) in a bottom block's mask (L) records a local peak at its
marked vertex in row (j). Reflecting the right bottom block horizontally
identifies its local histogram with the left one's without reversing
row indices. Consequently the original (K=2) condition requires
\[
 Z+\lvert L_{\rm left}\cap L_{\rm right}\rvert\le2.
\]
The same necessary condition holds for any collection of disjoint
cross-block edges whose endpoints have no other cross-block neighbors.
It supplements the old highest-label exclusions rather than treating
conditional seam and protected peak counts as independent within a block.

## Exact mask-conditioned local counts

In one bottom block (H), write (P) for the protected vertices, (B=H\setminus P),
and (J\subset B) for the six marked endpoints. Enumerate independent
sets (S\subseteq P\cup J) in the induced graph on (H). This eligible
region is exactly a (6\times6) rectangle, while its closed neighborhoods
are taken in the original (6\times7) block. Let (W(S)=S\cup N_H(S)\),
and let (C(S)) count local words making all of (S) peaks. Conditioning
on the largest label in (W(S)) gives
\[
 C(\varnothing)=|H|!,\qquad
 C(S)=|W(S)|^{-1}\sum_{r\in S} C(S\setminus\{r\}).
\]
Write (p=|P\setminus W(S)|\), (b=|B\setminus W(S)|\), and let (e\) count
edges between these two outside sets. Define
\[
 Q_P(S)=\sum_{r\in S\cap P} C(S-r),\quad
 Q_B(S)=\sum_{t\in S\cap B} C(S-t),\quad
 Q_{PB}(S)=\sum_{r\in S\cap P} Q_B(S-r).
\]

For each joint peak moment, classify the local highest labels. Let (T\)
require the largest label unprotected, (V\) require largest protected
and second largest unprotected and nonadjacent to it, and (A\) require
largest protected and second largest unprotected and adjacent to it.
The moment numerators for these restrictions are
\[
\begin{aligned}
 |H|T(S)&=b C(S)+Q_B(S),\\
 |H|(|H|-1)V(S)&=(pb-e)C(S)+b Q_P(S)+p Q_B(S)+Q_{PB}(S),\\
 |H|(|H|-1)A(S)&=e C(S)+\sum_{r\in S\cap P} C(S-r)
          |(N_H(r)\cap B)\setminus W(S-r)|.
\end{aligned}
\]
For (T\), the top vertex is either outside (W(S)\) or a specified
unprotected vertex of (S\). For (V\), split according to which of the
highest and second highest vertices belong to (S\). When neither belongs,
both must lie outside (W(S)\), and the count excludes adjacent ordered
pairs. When only the highest belongs to (S\), the second vertex lies
outside (W(S)\) by the specified nonadjacency. When only the second
belongs to (S\), the highest also lies outside (W(S)\). When both
belong to (S\), independence gives nonadjacency and both are removed
from the moment condition. For (A\), neither can belong to (S\) and
be adjacent unless both lie outside (W(S)\), or only the highest
belongs to (S\). In that latter case the second vertex must be adjacent
to the highest and outside (W(S-r)\). The other two membership cases
are impossible. All factors arise from fixing one or two distinct top
positions, after which the remaining relative order is uniform.

Aggregate the four modular moments (C,T,V,A\) by \((|S\cap P|,S\cap J)\).
A local word with (i\) protected peaks and exact marked mask (L\)
contributes \(\binom{i}{k}\) for each prescribed mask contained in (L\)
at protected moment size (k\). Upward Möbius inversion over the mask
and then binomial inversion over (k\) recover the exact local histograms.
Carrying out both inversions modulo each prime before CRT keeps the
reconstruction range at (42!\); the intermediate moments can exceed it.
The fourth required histogram, (H_i(L)\), counts words with both highest
labels protected. It equals (F_i(L)-T_i(L)-V_i(L)-A_i(L)\). This yields
exact (H\) data on the two bottom sides. The formula does not assume
that swapping the first two labels preserves the seam mask.

The dense cardinality counter stores the two residues (C(S),Q_B(S)\)
per entry, so the mixed predecessor term (Q_{PB}\) uses the preceding
layer alone. The term in (A\) uses those unprotected neighbors covered
uniquely by an (r\in S\cap P\). Hence the two-table resource bound remains
linear in the two largest adjacent layer populations and does not require
an extra predecessor layer.

## Joint contraction with the earlier statistical filters

Start with the archived exact top-block (F,T,V\) histograms from round
52, together with its certified lower count for the number of words with
both highest labels protected. No new top-block computation is needed.
Within each bottom mask pair, multiply the appropriate local histograms
and count interleavings using the multinomial coefficient for the block
sizes 60, 42 and 42. Keep only local protected peak counts and masks with
(Z+\lvert L_{\rm left}\cap L_{\rm right}\rvert\le 2\).

Apply the old top-one through top-five exclusions inside these surviving
mask cells. For (Z=2\) the seam masks now must be disjoint, and the earlier
extra-peak arguments and prescribed color prefixes still apply verbatim.
For (Z<2\) the existing program makes no highest-label exclusions. The
old top-block (H\) lower value can use any fixed subset of its true
both-protected words of that size separately for each local protected
peak count. The subset is specified independently of the bottom masks.
The rank exclusions are disjoint by the first unprotected label, and
by its block and the counts in the preceding prefix within a rank.
Consequently subtraction of their counted subsets is valid in every
mask cell. Exact bottom (H\) counts also improve the older conservative
cuts without adding any extra fit assumption.

## Exact finite results

Summing over exact masks reproduces the three archived bottom histograms
(F,T,V\) entry by entry. The merger reproduces the original round-52
upper count when it uses the old (H\) lower inputs and no mask constraint.
Using exact bottom (H\) alone tightens that bound by about 1.00303. Adding
the necessary seam condition tightens it by another factor of about
2.24598, for an overall improvement of about 2.25280. There are 21
nonempty local masks and 441 paired mask cells.

The exact contracted upper, divided by (144!\), gives
\[
 A_2(12)/144! < 9.544161373306617932\times10^{-26}.
\]
The same round-52 disjoint constructed numerators then give
\[
\begin{aligned}
 \mathbb P_2(\rho_{12}>4)&\ge 1.05415177848781451\times10^{-17},\\
 \mathbb P_2(\rho_{12}>6)&\ge 8.16459984072341637\times10^{-18},\\
 \mathbb P_2(\rho_{12}>7)&\ge 5.67341247957661947\times10^{-18}.
\end{aligned}
\]
The upper decimal is rounded upward, and the lower decimals downward.
The count is a strict finite upper bound, not a full (A_2\) enumeration.
It leaves the infinite-subsequence fixed-positive-mass question open.

## Gates and resource receipts

All local words in (2\times3\) and (3\times3\) rectangles were enumerated
as a separate check for every mask-resolved histogram and the ranking
shortcut. A three-block test graph has one eligible matched seam pair.
All (8!\) whole words verified the necessary condition against global
peak sets and the exact color-based contraction. The earlier filters
left 21,960 words in that gate, the joint filter left 21,840, and the
true two-peak count was 9,680.

Three certified 61-bit moduli are sufficient for the 42-vertex bottoms.
Each modular worker actually processed 5,598,861 independent sets and
49,205,952 predecessor queries. The maximum live two-table payload was
38,252,064 bytes per modulus. All three workers completed and were
waited in 0.41 seconds wall with observed maximum per-worker RSS
62,111,744 bytes. The final merge and independent gates took 3.57
seconds wall. No new top-block table was allocated and no remote or
GPU work was invoked by this package.

On the data volume build with `clang++ -std=c++20 -O3 -Wall -Wextra
-Werror joint_seam.cpp -o joint_seam`, using the shared module cache.
Generate the local plan and use Python 3.11 or later to run the bounded
jobs and merge:

```sh
./joint_seam 6 7 0 --plan > n12-plan.json
python3 -B run_joint.py <new-batch-dir>
python3 -B merge_joint.py <batch-dir> <new-ledger.json> --prior <r052-t5-dir>
```

The prior directory contains `r052-denominator-refined.json` and
`combined.json`. The local folder contains `pratt-61bit.json`. The
batch driver uses macOS `/usr/bin/time -l` for RSS measurement. The
ledger records exact counts, fractions, calibration data, and hashes
of the original numerators and denominators.
