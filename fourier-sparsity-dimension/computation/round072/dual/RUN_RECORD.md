# R72 weighted bound on the complete Walsh rank

For each Boolean sign function \(F(t,z)\), with three or four bits in
each of \(t\) and \(z\), the full integer Walsh transform is grouped by
its \(z\)-frequency \(\beta\). Let \(n_\beta\) be the actual number of
points in a nonempty band, let \(c_\beta\) be the affine dimension of
its \(t\)-frequency set, and write \(m_\beta=2^{c_\beta}\),
\(M=\sum_\beta m_\beta\), and \(s=\sum_\beta n_\beta\). Both
implementations independently compute the actual active label set at
each \(t\) from the true \(z\)-slice, and compute the linear rank \(d\)
of the entire Walsh support. The generator uses character, two-form
quadratic, arbitrary Boolean and mixed true slices. Case IDs are
deterministic and can yield repeated models.

For four choices of positive weights, the two implementations check
\[
d^2\le 7\left(\max_t\sum_{\beta\ \mathrm{active\ at}\ t}w_\beta\right)
              \left(\sum_\beta m_\beta/w_\beta\right).
\]
They also check the stronger form in which the first factor uses only
one entire affine coset from each band's active set. For each band the
coset is selected by the least active address and the band's actual
affine direction space. The programs verify containment and its exact
size. The weighted factors use integer arithmetic, clearing denominators
with \(720720=\operatorname{lcm}(1,\ldots,16)\) and
\(840=\operatorname{lcm}(1,\ldots,8)\).

The choices of \(w_\beta\) are 1, intrinsic \(m_\beta/n_\beta\), intrinsic
times a seeded multiplier from 1 through 8, and the reciprocal of the
number of active addresses. With intrinsic weights the second factor
is exactly \(s\). The unweighted product is \(LM\), where \(L\) is
the actual maximum slice sparsity. The counts compare the full-activity
products before the shared factor 7.

| Weights | Strictly smaller than \(LM\) | Equal | Strictly larger than \(LM\) |
| --- | ---: | ---: | ---: |
| Intrinsic | 86,710 | 7,175,112 | 7,738,178 |
| Intrinsic with varying multiplier | 880 | 609 | 14,998,511 |
| Reciprocal of active address count | 22,623 | 148,699 | 14,828,678 |

The run covered 15,000,000 cases, including 7,835,093 with at least one
holed band. For intrinsic weights the selected cosets strictly improve
the full-activity bound in 7,156,252 cases. The exact aggregate and
positional checksum match CPU and Metal on both hosts.
The archive includes 8,192 independent Fraction audits and five
explicit witnesses. Case 959 reduces \(LM=1856\) to an intrinsic
weighted product of \(1800\). Case 1540 increases \(LM=104\) to \(752/3\).
Their respective ratios \(225/232\) and \(94/39\) are the extremes in
the independently checked pilot. For each witness a separate exact
dynamic program finds the minimum weighted normal cost
\(\sum_\beta w_\beta(1+d_\beta)^2/m_\beta\) among allocations of
independent subspaces spanning all band directions, and certifies it
against seven times the selected-coset overlap. The reported ratios
are finite-sample comparisons between the two prescribed weight choices.

## Runtime and resources

| Measure | Mac mini | MacBook |
| --- | ---: | ---: |
| Case-ID range | \([0,9000000)\) | \([9000000,15000000)\) |
| CPU logical workers | 12 | 15 |
| CPU task wall time | 8.591 s | 1.591 s |
| Project CPU mean logical-core equivalents | 4.23 | 14.22 |
| Project CPU peak logical-core equivalents | 7.56 | 14.74 |
| Whole-system CPU idle, 0% samples | 53 of 53 | 8 of 10 |
| Metal device | Apple M6 | Apple M5 Pro |
| Completed Metal commands | 1,099 | 733 |
| Summed command duration | 2.130 s | 1.605 s |
| Whole-device GPU utilization median / peak | 63% / 91% | 100% / 100% |
| Whole-device GPU utilization before launch | 7% | 99% |
| Matching CPU and Metal checksum | ee111972b9ea5900 | ba362e757370cbae |

CPU data use the process CPU clock and system tick counters. Intervals
without changed system counters supply only process measurements. The
whole-device GPU percentages can include unrelated work. Project GPU
work is established by Metal receipts. The two hosts' CPU and GPU work
overlapped across all four resources for 0.382 seconds over 101 intervals
using their host clocks. All four compute processes exited with status
0, and a subsequent PID check found none still active.

## Reproduction

The hosts used macOS 27.2 (26B123) and 27.3 (26C5035r), Apple clang
21.0.0, and Apple Swift 6.4.2. Python 3.9.6 supervised both hosts and
Python 3.15.0 ran the archive verifier. The verifier needs Python 3.10
or above.

~~~sh
clang -O3 -std=c11 -pthread weight_cpu.c -o weight_cpu
swiftc -O WeightMetal.swift -o weight_metal
./weight_cpu 0 8192 12 0 pilot_cpu.json pilot_cpu.csv
./weight_metal pilot_gpu.json 0 8192 4096
python3 -B verify_weights.py local macbook reproduced.json
cmp reproduced.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

To repeat both main ranges, compile on both hosts and run
`python3 -B run_host.py output 0 9000000` on the Mac mini and
`python3 -B run_host.py output 9000000 6000000` on the MacBook.
Use fresh output directories. Once both write `ready.json`,
write the same future Unix time into both `start_epoch.txt` gates.
The archived start was 1794980954.623257.
