# Pairwise disjoint Lipschitz cut events with prescribed row parities

Let (G=P_n\square P_n) be the free-boundary grid, with coordinates
(0\le x,y<n). All (n^2) vertex labels are distinct. For
(h=(h_y)_{y=0}^{n-1}\in\{1,\ldots,n-1\}^n), put
\[
L_h=\{(x,y):x<h_y\},\qquad R_h=V\setminus L_h.
\]
Both sides are connected, since every occupied row meets the respective outer
edge. Write (d=\lceil n/4\rceil). Define (E_h) to require exactly one peak
in each induced side, in (x<d) on the left and (x\ge n-d) on the right.
Require the two roots to be safe, meaning all of their full-grid neighbors
are in their own side. Every nonroot has a higher neighbor inside its own
side, and the side roots remain full-grid peaks. In particular, (E_h\)
implies (K=2) and original Manhattan distance at least
(n-2d+1\ge(n-1)/2), hence greater than (n/3) for (n\ge4).

The disjointness argument also works with an alternative root rule: insist
that each side contains exactly one prescribed full-grid peak and that this
peak is the unique peak of the induced side. Cross-interface neighbors of
these roots are allowed. The argument uses the common assignment of the
full-grid roots to the two sides, supplied either by the outer root bands
in (E_h) or by prescribing a left root and a right root. If assignments
may switch between cuts, this argument does not give disjointness.

## Difference components

Suppose a word belongs to both cut events and has the same left and right
roots for both cuts. Put (U=L_h\cap L_k) and (W=R_h\cap R_k). In a connected
component (D) of (L_k\setminus L_h), let (v) have largest label. This
vertex is a nonroot in both (L_k) and (R_h). Inside either side, a higher
neighbor of (v) cannot lie in (D), and another component of the difference
set has no edge to (v). Thus its higher neighbors for (L_k) and (R_h)
respectively lie in (U) and (W). The same statement applies to components
of (L_h\setminus L_k). Every nonempty difference component must contain a
vertex with neighbors in both common sides. This is the round-42 bridge
condition, stated here for common full-grid roots.

If both (h) and (k) have consecutive row-width differences of absolute
value at most one, consider row (y) with (h_y<k_y) and a difference vertex
(v=(x,y)), where (h_y\le x\le k_y-1). A horizontal neighbor in (U) is
possible only at (x=h_y). For a vertical neighbor in (U), one needs
(x<h_{y\pm1}\le h_y+1), again implying (x=h_y). A horizontal neighbor in
(W) is possible only at (x=k_y-1). A vertical neighbor in (W) requires
(x\ge k_{y\pm1}\ge k_y-1), so also (x=k_y-1). Therefore a bridge vertex is
possible only if (k_y-h_y=1), in which case the only difference vertex in
that row is a bridge by its two horizontal neighbors. Interchanging (h,k)
gives exactly the same result in the other difference set.

Fix any row parity word (p\in\{0,1\}^n). If (h_y\equiv k_y\equiv p_y\pmod2\)
for all (y), every differing row has absolute width difference at least
two. For distinct cuts there is a nonempty difference component and no bridge
vertex in it, so the two events are disjoint. No distance estimate or counting
identity is needed for this disjointness step.

## An exponentially large subfamily with common safe bands

Take (h_0=\lfloor n/2\rfloor), allow each of the (n-1) steps to be (+1\) or
\(-1\), and keep the paths with all widths in \([d+1,n-d-1]\). All widths
have a fixed parity by row. The entire left and right outer bands are safe
for every retained cut, by their one-column clearance from each interface.
Let (r_n=\min\{\lfloor n/2\rfloor-d-1, n-d-1-\lfloor n/2\rfloor\}\).
For independent equiprobable signs, the probability of any partial sum
exceeding (r_n) in absolute value is at most
\[
 2(n-1)\exp\left(-\frac{(r_n+1)^2}{2(n-1)}\right)
\]
when (r_n\ge0). This follows from \(\cosh t\le e^{t^2/2}\), the
exponential Markov inequality in each direction, and the union bound over
partial sums. Since (r_n/n\to1/4), the family has
(2^{n-1}(1-O(n e^{-c n}))) cuts for some absolute (c>0). Each event is
nonempty: choose the two corner roots on opposite sides, grow each connected
side from its root, then shuffle the two descending orders preserving their
internal orders. Thus there are exponentially many disjoint nonempty far-peak
events. This does not bound their weight sum relative to (A_2(n)).

## Exact finite counts and implementation

For a side (H) and its set (A) of safe outer-band roots, define (C_H(S))
for independent (S\subseteq H) as the number of orders making every (v\in S)
a peak. By conditioning on the position of the largest label in the closed
neighborhood of (S),
\[
 C_H(\varnothing)=|H|!,\qquad
 C_H(S)=\frac{\sum_{v\in S} C_H(S\setminus\{v\})}{|S\cup N_H(S)|}.
\]
Indeed, the largest vertex of this neighborhood under simultaneous peak
constraints is in (S). Given a specified (v\in S) is the maximum of the
neighborhood, the relative order of its remaining vertices is uniform. Also,
independence means (v) is not in the closed neighborhood of (S\setminus\{v\}).
The recurrence follows and every quotient is integer. Inclusion-exclusion
over independent peak sets gives the restricted single-peak count
\[
 a_H(A)=\sum_{\varnothing\ne S\text{ independent}}(-1)^{|S|-1}|S\cap A|C_H(S).
\]
Cut event counts equal
\(\binom{n^2}{|L_h|} a_{L_h}(A_L)a_{R_h}(A_R)\), since both side orders can be
interleaved arbitrarily. The right side count is computed as a left side of
width sequence \(n-h\) by reflection. In a fixed parity family these
counts add exactly. The proof and code count the safe-root events above,
not the enlarged events with cross-interface root neighbors.

The implementation `cut_region.cpp` counts each side using 192-bit integer
arithmetic. The code caps the region size at 36. Each intermediate sum in
the recurrence is at most (36\cdot36!\), and the two signed inclusion-exclusion
sums are individually at most (36\cdot2^{36}\cdot36!<2^{180}\). The C++
program asserts no overflow, nonnegative subtraction, and zero division
remainders. Python uses arbitrary-precision integer arithmetic for the final
products, sums and rational probabilities.

`disjoint-ledger.json` contains the parity-of-active-row construction with
all inactive widths fixed at (m=\lfloor n/2\rfloor) and the active widths
chosen independently in \(\{m-1,m+1\}\). Phases separately use the zero-based
even and odd rows. For the even-active families the exact results are

| \(n\) | Events | Union count | Union probability given \(K=2\) |
|---|---:|---:|---:|
| 4 | 4 | 25170390960 | 0.03229875225995825 |
| 5 | 8 | 413190355925965577400 | 0.05849852182539256 |
| 6 | 8 | 12908537561486409541605351981360000 | 0.04241319329549939 |
| 7 | 16 | 3794533811016162859403422392783698307056558727862272 | 0.02768158115605677 |

The n=5 and n=7 odd-active probabilities are respectively
0.03087253239366084 and 0.01642330994566399. At even n the two phases agree by
row reversal. The ledger records the exact normalizer sources and SHA-256.
C++ rechecked every n=5 alternating-row term against the independent round-42
Lipschitz diagonal and all 128 old r076 n=7 central binary event counts, also
checking the two region counts, sizes, allowed root counts, independent set
counts and independence numbers for each old curve. The two n=7 families
with fixed inactive width (m) do not assert disjointness across phases.

`parity-ledger.json` records the additional families:

| \(n\) | Widths and row parities | Curves | Union count | Union probability given \(K=2\) |
|---|---|---:|---:|---:|
| 5 | All widths 1..4, maximum total over the 32 fixed parity classes | 8 | 623157272838614660640 | 0.08822514563319939 |
| 6 | All widths 1..5, alternating parity from either initial parity | 36 | 28255911773169508092723435994512000 | 0.09283959875916552 |
| 7 | Widths 2..5, alternating parity from either initial parity | 34 | 7692330776915440242599079443568689944840065826997760 | 0.05611647946375644 |

The n=5 totals use all 178 round-42 single-event entries (reconstructed from
the diagonal orbit representatives). The n=6 and n=7 totals are new exact
region computations. The two alternating parity variants are reported
separately; their common value follows by a reflection and their masses are
not added. The n=8 run uses the even-active fixed width 4 family with 16
profiles. It records an unconditional probability with denominator 64!,
without an exact (A_2(8)) normalizer.

## Replay

Use a fresh build directory on the data volume and Apple clang++ with
`-std=c++20 -O3 -Wall -Wextra -Werror`. Compile `cut_region.cpp` to a sibling
binary named `cut_region`. The other scripts require Python 3 standard library
only, with the original round-42 and r076 ledgers at the paths shown in the
scripts. From the authoring directory, run `python3 -B certify_disjoint.py
--workers 3`, `python3 -B certify_parity.py`, and `python3 -B run_n8.py
<fresh-output-directory> --workers 3`. The largest n=8 region profile
`54545454` has 36 vertices, 10,213,044 independent sets and exact rooted count
61909885147738447311261252268248. The local largest-input gate measured
789,987,328 bytes maximum RSS and 3.22 seconds on Apple clang 21.0.0. The n=8
driver independently checks all row-mask independent-set totals and records
per-child process receipts. No empirical or asymptotic conclusion is inferred
from these finite counts.

The completed MacBook n=8 ledger is at
`/Volumes/data/CodexTemp/T2/a022-n8-disjoint-remote-r082/results/summary.json`.
It reports 16 regional counts, 120 geometric non-intersection checks and
all 16 child processes waited. The union count is
15497852066351452649557374277633220714870471956684105383340393957380983808,
with unconditional probability about \(1.2213907137141383\times10^{-16}\).
