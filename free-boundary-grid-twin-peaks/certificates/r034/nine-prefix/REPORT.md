# N=9 exact finite certificate report

The certified program used integer subset DP for every ordered nonedge
\((r,p)\) of each fixed bipartite host and records the first complete rank of a
neighbor of \(p\).  The finite relation checked per graph is

\[
R_{27}=h_2-5h_4+5h_6-2h_7.
\]

The frozen evidence covers all `1,327,360` edge subsets of the fixed hosts
`1+8`, `2+7`, `3+6`, and `4+5`; it includes disconnected graphs and isolates.
Strict aggregation accepts completed intervals only.  The portable auditor
actually reconstructed the ranges from the 326 raw checkpoints, obtaining
`[0,256)`, `[0,16384)`, `[0,262144)`, and `[0,1048576)` respectively, with no
gap or overlap, all source hashes equal to the certified source, and all
`max_abs_R27 = 0`.  The 326th interval is the counted left-4 benchmark
`[0,1024)`; the other 325 have their own recorded successful exits.

The source compiled historically with `-O3 -std=c++20 -Wall -Wextra -Werror
-pedantic`.  Receipts record successful compile, two gates, production exits,
strict aggregation, and completion.  The reported resource measurement is
only an uncounted isolated left-4 `[0,1024)` reprofile: 0.14 s wall, 0.13 s
user, 0.00 s sys, 2,195,456-byte maximum resident set, and 1,556,888-byte peak
footprint.  It must not be read as a full-batch resource total.

The N=7 `R45` gate enumerated 5,184 fixed-host edge subsets.  The N=9 gate
compared DP to direct `9!` enumeration on exactly `P9`, `K23-tail`, and
`C6-tail`; all three agreed and reported zero.  This is a targeted
implementation cross-check only.  There was one full N=9 subset-DP scan, not
two independent N=9 full scans.

The finite result is an exact computer-assisted local statement.  The separate
analytic lifting argument in `proof/NINE_PREFIX_LIFTING.md` explains what it
implies for connected bipartite graphs with \(N\ge9\), conditional on the
finite certificate.  Neither component proves the original asymptotic problem
or a peak-distance conclusion.
