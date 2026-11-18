# An exact finite upper bound for the two-peak normalizer

Take a 1-Lipschitz row-width sequence \(h_y\in\{1,\ldots,n-1\}\) and its
left and right induced grid regions. Call a vertex protected if all its
original full-grid neighbors are in its own region. In the left side, the
protected vertices in row \(y\) are exactly \(x<h_y-1\). The excluded
rightmost vertex has a horizontal neighbor across the cut. Every other
vertex has both of its existing horizontal neighbors within the side and
also satisfies \(x<h_{y-1},h_{y+1}\) for existing adjacent rows by the
Lipschitz inequalities. The right side statement follows by reflection.
Each side has one excluded vertex in each row, and the combined protected
sets have \(n(n-2)\) vertices. At any protected vertex the local peak test
is the original full-grid peak test.

Write \(Z\) for the number of protected peaks on both sides. Global peak
count \(K=2\) implies \(Z\le2\). If \(Z=2\), all global peaks are protected.
The largest grid label must therefore be protected. If the second largest
label is unprotected, it must be adjacent to the largest label. Otherwise
it is a third global peak, since only the largest label can exceed it.
These conditions are deterministic necessary conditions for \(K=2\).
The method counts all words satisfying these conditions and uses that count
\(U\) as a finite upper bound for \(A_2(n)\).

## Side histograms from independent peak sets

For one side \(H\), let \(a=|H|\), let \(P\) be its protected vertices and
let \(B=H\setminus P\). A side word gives the descending relative order of
labels in \(H\). For each protected peak count \(i\), use three histograms:

- \(F_i\) counts all side words with exactly \(i\) protected peaks.
- \(T_i\) counts such words whose largest label is at an unprotected vertex.
- \(V_i\) counts such words with largest label at a protected vertex and
  second largest label at an unprotected vertex nonadjacent to the largest.

To compute these histograms, choose an independent set \(S\subset P\). Let
\(W=S\cup N_H(S)\) be its closed neighborhood in the **whole side**. Let
\(C(S)\) count side words making every vertex in \(S\) a peak. For independent
sets, conditioning on the highest label in \(W\) gives
\[
 C(\varnothing)=a!,\qquad |W|C(S)=\sum_{s\in S}C(S\setminus\{s\})\quad(S\ne\varnothing).
\]
The highest label in \(W\) belongs to \(S\) whenever all of \(S\) peaks.
Independence implies that its vertex lies outside the closed neighborhood
of the remaining peak set. This proves the integer recurrence under the
uniform relative order of labels in the side.

Among the words counted by \(C(S)\), having largest label at a given
unprotected vertex is possible only if that vertex lies outside \(W\).
Every such vertex is unconstrained by the peak tests, and fixing it as
largest divides the count by \(a\). Thus the restricted peak moment for
\(T\) is \(|B\setminus W|C(S)/a\).

For the \(V\) moment, write \(r\in P\) and \(t\in B\) for the positions of
the two largest labels. They must be nonadjacent. If \(r\notin S\), neither
can belong to \(W\). The number of eligible ordered pairs of vertices is
\[
 |P\setminus W|\,|B\setminus W|-e_H(P\setminus W,B\setminus W),
\]
where \(e_H\) counts the side edges between its arguments. Fixing each such
pair leaves exactly \(C(S)/(a(a-1))\) side words. If \(r\in S\), it must
peak automatically and the remaining constraints involve \(S\setminus\{r\}\).
The vertex \(t\) must be outside the neighborhood of these remaining peaks
and must be nonadjacent to \(r\), so it must be outside \(W\). For each fixed
\(r,t\) there are \(C(S\setminus\{r\})/(a(a-1))\) completions. Sum over
\(r\in S\) and apply the recurrence to obtain the single formula for the
\(V\) moment,
\[
 \frac{C(S)}{a(a-1)}\left[(|P\setminus W|+|W|)|B\setminus W|
                  -e_H(P\setminus W,B\setminus W)\right].
\]
For \(S=\varnothing\) the term involving \(|W|\) is zero, as expected.

For each of the three statistics, sum the corresponding moment over
independent sets of a fixed size \(k\). For any actual word with \(i\)
protected peaks, it contributes \(\binom{i}{k}\) to that aggregate, since
simultaneous peaks are independent. Binomial inversion therefore recovers
its exact histogram. The C++ two-layer modular counter enumerates only
independent sets of \(P\), retaining closed neighborhoods in \(H\), and
accumulates the three moments. The product of the three certified primes
is larger than \(41!\), which bounds each reconstructed histogram entry.

## Exact merge and calibration

Let \(a\) and \(b\) be the two region sizes and \(N=a+b\). Once the two
relative side words are fixed, there are \(\binom{N}{a}\) interleavings.
First count all words with \(Z\le2\) by convolving the \(F\) histograms
and multiplying by this interleaving count. For every \(i+j=2\), remove
these disjoint subcases:

1. The global maximum is unprotected. It lies in the left or right side,
   giving removed counts \(\binom{N-1}{a-1}T^L_i F^R_j\) and
   \(\binom{N-1}{a}F^L_i T^R_j\).
2. The global maximum is protected and the global second largest label
   lies in the other side at an unprotected vertex. The two cannot be
   adjacent because the maximum is protected. The removed count is
   \(\binom{N-2}{a-1}[(F^L_i-T^L_i)T^R_j+T^L_i(F^R_j-T^R_j)]\).
3. The global maximum is protected and the global second largest label is
   an unprotected nonneighbor in that same side. The removed counts are
   \(\binom{N-2}{a-2}V^L_iF^R_j\) and
   \(\binom{N-2}{a}F^L_iV^R_j\).

The number left by these exact subtractions equals \(U\). Complete
permutation enumeration on six small side regions verified all three
histograms. Full \(9!\) enumeration on the original \(3\times3\) grid verified
all three merged upper counts and the true \(A_2(3)\). The final third-stage
bound also lay above the exact archived \(A_2(5)\) and \(A_2(7)\).

## Size nine

Among four tested size-nine cuts, the strongest bound came from
`444445555` with right reflection `555554444`, and region sizes 40 and 41.
It has 63 protected vertices and requires only 799,989 and 1,202,391
independent sets on the two sides. The final exact value is
\[
 U_9=1340033358218795995643106830381546075076376640233944809266307521415217609541972784018289543651035170465382400,
\]
so \(A_2(9)/81!\le U_9/81!\approx2.31154774524987524\times10^{-13}\).
All 60 local modular computations for the calibrated batch completed and
were waited. Maximum per-worker RSS was 61,505,536 bytes and local wall
time for the complete batch was 12.06 seconds.

The independently computed disjoint root-stratified numerator from the
28-cut construction contains 1296 eligible root strata before the distance
filter in the supplied ledger. Using the exact integer \(U_9\) in the
denominator gives
\[
 \mathbb P_2(\rho_9>3)\ge1.0798401479712629\times10^{-8},\qquad
 \mathbb P_2(\rho_9>9/2)\ge1.0239552761846083\times10^{-8}.
\]
Each displayed decimal is rounded downward. The second bound retains only
root pairs at Manhattan distance at least five. These are finite lower
bounds. Exact \(A_2(9)\) remains unknown, and the bound does not control the
normalizer along an infinite subsequence.

## Replay

Build `protected_mod.cpp` with Apple `clang++ -std=c++20 -O3 -Wall -Wextra
-Werror` to a sibling binary named `protected_mod`. With the three-prime
certificate `pratt-61bit.json` and the reused `reference_crt.py`, run
`python3 -B count_upper.py <fresh-dir> --workers 3`. The output directory
contains per-prime JSONs, resource logs and the exact histogram summary.
`combine_family_bound.py` takes that summary and the separate exact
construction-family summary and writes exact fractional conditional bounds.
