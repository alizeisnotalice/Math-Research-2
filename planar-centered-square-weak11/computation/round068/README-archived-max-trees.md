# R68: exact support-tree classification of archived global maxima

This is one read-only classification of the twelve pinned R63 deletion-batch
files, containing 72 inputs and 211 global maximum eligible masks (199 proper,
12 equal to the full input). It does not enumerate every component's maximum
blocks. It evaluates neither a halo nor a cut.

The R63 global maximum-mask lists and eligible-hypergraph component partitions
are inherited metadata, bound by SHA-256 and byte count in
`archived-max-tree-sources.json`. The program independently checks positive
input masses, coordinates, mask ranks, masses, strict eligibility and full bbox
closure. It also checks each saved component's bbox/mass ratio, its terminal
flag, the component partition, and which component contains the selected mask.
It does not regenerate the eligible-hypergraph partition or maximum directory.

For a fixed maximum mass L, every edge has cost P sqrt(L)-Q, where
P=|dx|+|dy| and Q=|dx dy| are nonnegative rationals. Kruskal comparisons reduce
the sign of a sqrt(L)-b to rational signs and, only when both a,b are positive,
the sign of L a^2-b^2. Ties are settled by source indices. No floating ordering
is used. For the resulting minimum spanning tree, payment is equivalent to

    L P_total^2 <= (5 L/2 + Q_total)^2.

Both sides before squaring are nonnegative. Each output row retains the exact
tree edges, P_total, Q_total, and the signed rational difference of these two
squares. A failed test proves only that this support-tree sufficient criterion
fails for this particular global maximum block.

The priority classes are mutually exclusive, in this order:

1. The **actual containing eligible component** satisfies d_component^2 <=
   35 m_component/33. This is not the whole-input terminal condition.
2. min(width_x,width_y)^2 <= 9 L/16.
3. Support rank <= 3.
4. Minimum support-tree cost <= 5 L/2.
5. None of these criteria applies.

The program also reports all overlapping flags separately. In order above,
the 211 rows split as **47, 120, 12, 26, 6**. For the 199 proper rows the
split is **35, 120, 12, 26, 6**. There are 192 tree-paid rows overall; this
overlapping total must not be described as 192 new payments. At input level,
59 have some global maximum covered by the first three criteria and 67 by
all four. Five inputs have all of their global maximum ties unclassified.
Those five inputs contribute six remaining masks, all in connected,
nonterminal inputs, with ranks 13--15 and bbox widths 3 by 3.

The actual invocation used Python 3.15.0. Script wall time was approximately
0.13236 seconds; the child exited with code zero and was absent after waiting.
The full command, source/output hashes, and actual process receipt are saved.
There was one classification run, no MATLAB, no halo/cut replay, and no new
input generation. This late author classification is separate from the frozen
independent R68 geometry review and is not claimed as an additional independent
reconstruction of R63's source hypergraph.

From the topic directory after integration, a portable replay is:

```sh
python3 computation/round068/classify_archived_maximum_trees.py \
  --archive-root . \
  --spec computation/round068/archived-max-tree-sources.json \
  --output /Volumes/data/CodexTemp/T2/r068-max-tree-replay.json
```

The twelve archived inputs must match their pins; no old workspace fallback is
implemented. Current general W-prime selection and the exact C2 problem remain
open.
