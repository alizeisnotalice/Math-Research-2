# R67 basis occupancy and four-relation graph audit

The addresses and data variables each have dimension three or four. The
generated eight- or sixteen-row Boolean selector tables use characters,
quadratic phases with two independent linear forms, arbitrary small
Boolean truth tables, and mixtures of these. Case IDs determine
the generator. Distinct IDs may give the same table.

Let \(S_y\) denote the actual support of the Walsh spectrum of slice \(y\),
and let \(B\) be the union of these slice supports. Independently, the
whole truth-table Walsh spectrum groups its points by their data label
\(\beta\). The nonempty groups give precisely \(B\). Write \(n_\beta\)
for the size of group \(\beta\), and \(q=\sum_{\beta\ne0}n_\beta\).
For every basis \(I\subset B\setminus\{0\}\) of the data-label span, the runs
check, with \(r=|I|\) and \(\mu_I=\max_y|I\cap S_y|\),
\[
r^2\le \mu_I\sum_{\beta\in I}n_\beta\le q\mu_I.
\]
Indeed the activity fraction \(p_\beta\) of the \(\beta\) coefficient
function satisfies \(n_\beta p_\beta\ge1\) by Cauchy–Schwarz and Parseval.
Averaging the count of active basis labels gives
\(\sum_{\beta\in I}p_\beta\le\mu_I\). Cauchy–Schwarz on these two lists gives
the stronger first inequality.

Graph \(A\) joins pairs in a four-point zero-sum relation within at least
one actual Boolean slice. Graph \(P\) joins pairs in a four-point zero-sum
relation within \(B\cup\{0\}\). Four distinct labels are required in each
relation. Both full graph clique counts use only the nonzero vertices.
Within every actual slice, the Boolean square identity gives an alternate
same-xor pair for each pair of active labels, so the slice itself induces a
clique in \(A\). The experiments verify \(\mu_I\le\omega(A[I])\le\omega(P[I])\)
for every basis. They also compare the two unrestricted clique counts.

| Statistic | Mac mini | MacBook |
| --- | ---: | ---: |
| Case-ID range | \([0,4000000)\) | \([4000000,7000000)\) |
| Individual bases checked | 1,243,933,928 | 932,853,185 |
| Cases with a strict projected-to-actual clique gain | 1,926,973 | 1,445,265 |
| Cases with some basis having \(\omega(A[I])>\mu_I\) | 1,088,951 | 816,677 |
| Maximum full-graph clique gap | 14 | 14 |
| Maximum basis-induced clique gap | 3 | 3 |
| Matching CPU and Metal checksum | 4c005faa3cd6b6bf | 50a7cef8231357dc |

The verified archive contains complete gap histograms for the full
graphs and all induced bases. The independent Python checker redoes
the 4,096-case pilot and its 1,264,729 basis comparisons using explicit
four-point relations and separate clique search. A control with a
different character in each of sixteen slices gives an edgeless actual
graph and a complete projected graph on fifteen nonzero labels. Its
coordinate basis has occupancy and actual clique size 1 and projected
clique size 4. A separate control runs through the six pairs of four
independent labels and uses a quadratic phase on each pair. The
coordinate basis has actual and projected clique size 4 but occupancy 2.

## Device receipts

| Measure | Mac mini | MacBook |
| --- | ---: | ---: |
| Project CPU workers | 12 | 15 |
| Project CPU task wall time | 8.743 s | 9.465 s |
| Project average logical-core equivalents | 6.98 | 5.80 |
| Project peak logical-core equivalents | 7.55 | 9.52 |
| Whole-machine CPU idle, valid full-worker samples | 0% in 36 samples | 0% in 61 samples |
| Metal device | Apple M6 | Apple M5 Pro |
| Completed Metal commands | 1,954 | 1,465 |
| Sum of Metal command durations | 8.500 s | 11.586 s |
| Metal workload wall time | 4.308 s | 7.461 s |

The whole-machine idle metrics have no sample when the system tick count
does not change. Median and peak whole-device GPU utilization were 100% on both
hosts and the baseline was 0% and 99% respectively. Device readings include
unrelated work. The command receipts establish project GPU work. Sums of
command durations include possible overlap from two queues per host.
The four CPU and GPU execution windows overlapped for 3.273 seconds
across 210 positive-length intersections, using the host clocks.
All four compute processes exited with status 0 and their PIDs were absent
upon a subsequent check.

## Reproduction

Operating systems were macOS 27.2 (26B123) and 27.3 (26C5035r) for the Mac
mini and MacBook respectively. Both hosts used Apple clang 21.0.0, Apple
Swift 6.4.2 and Python 3.9.6 for the supervisors. Python 3.15.0 ran the
independent verifier, which requires Python 3.10 or above.

~~~sh
clang -O3 -std=c11 -pthread clique_cpu.c -o clique_cpu
swiftc -O CliqueMetal.swift -o clique_metal
./clique_cpu 0 4096 12 0 pilot_cpu.json pilot_cpu.csv
./clique_metal pilot_gpu.json 0 4096 1024
python3 -B verify_cliques.py local macbook reproduced.json
cmp reproduced.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

To rerun the main batch, copy the same sources to both hosts and build
both executables. Execute `python3 -B run_host.py output 0 4000000` on
the Mac mini and `python3 -B run_host.py output 4000000 3000000` on
the MacBook. After the supervisors write their `ready.json` receipts in
fresh output directories, write the same future Unix time to both
`start_epoch.txt` gates. The actual shared start was 1794978530.691055.
