# An all-size far-peak map by protected seam promotion

For (n\ge4\) divide each axis into lengths \(\lfloor n/2\rfloor\) and
\(\lceil n/2\rceil\). This gives four rectangles (H_i\) of sizes (M_i\).
For each block let (B_i\) consist of the two inward-facing edges.
Their union is a connected L shape in the block. It is exactly the set
of vertices with a full-grid neighbor outside the block. The total size
of the four seams is (4n-4\). There are ((n-2)^2\) protected vertices
outside their union, and their local and global peak tests agree.

Take any descending word with exactly two global peaks, and let (P\)
be the set of its protected local peaks, so \(|P|\le2\). In each block form
(S_i\) from (B_i\) and a deterministic shortest path from each peak in
(P\cap H_i\) to (B_i\). Each path contributes at most
\(\lceil n/2\rceil-1\) new vertices, and every (S_i\) is connected.
Consequently
\[
 \sum_i|S_i|\le 4n+2\lceil n/2\rceil-6.
\]

Fix a reference root (q_i\) at the horizontal paired seam on the
outermost left or right column of the block. Starting at (q_i\), write
(S_i\) in any deterministic connected order. Move this whole order
before the remaining vertices, while preserving the order on (H_i\setminus S_i\).
Every block peak of the original word is in (S_i\). Thus any vertex
outside (S_i\) originally had an earlier same-block neighbor. In the
new block word its earlier neighbor is either in the promoted prefix
or still earlier in the unchanged suffix. The new block word therefore
has exactly one peak, at (q_i\).

Keep the original interleaving of the four block indices. In the resulting
whole word, only the four (q_i\) can be peaks. Their cross-block neighbors
are exactly their aligned partners, since each (q_i\) is away from the
vertical seam between the two column halves. In each column half, exactly
the larger of its two aligned roots remains a peak. The two true peaks
are on the outermost left and right columns, hence have distance at
least (n-1\). This maps the entire original (K=2\) set to a subset of
the aligned-root far-peak family.

## Inverse capacity

Given an output and (P\), all the marked sets (S_i\) are determined.
For each (i\), encode the original slot within its local word of every
vertex of (S_i\), listing those vertices in a fixed order. There are at
most \((M_i)_{|S_i|}\) possible codes. Fill the remaining slots with the
remaining vertices in the order retained in the output's local suffix,
and keep the output's interleaving of block indices. This uniquely
recovers the original input. Thus the maximum fiber is at most the
explicit integer
\[
 D_n=\sum_{\substack{P\subseteq V\setminus\bigcup_i B_i\\ |P|\le2}}
             \prod_{i=1}^4 (M_i)_{|S_i(P)|}.
\]
In particular, setting (M_*=(\lceil n/2\rceil)^2\) gives the true
original-conditioned bound
\[
 \mathbb P_2(\rho_n\ge n-1)\ge D_n^{-1}
 \ge\left[
       \left(1+(n-2)^2+\binom{(n-2)^2}{2}\right)
       M_*^{4n+2\lceil n/2\rceil-6}
      \right]^{-1}.
\]
It follows that \(-n^{-2}\log\mathbb P_2(\rho_n\ge n-1)\to0\) for the full
sequence. The bound allows the probability to tend to zero and gives
no fixed positive mass on an infinite subsequence. It is a direct
quadrant counterpart to existing approaches at the same surface scale.

The gate `python3 -B promote_seams.py <output-json>` enumerates the slot
capacity for each test size and verifies the map and exact inverse on
2,000 generated genuine two-peak words across sizes 4, 5, 6, 7 and 12.
The fixed seed and witnessed second-peak rank range are recorded. The
sample is used only as a replay check. The probability inequality above
uses the deterministic fiber bound, not any estimate from these words.
