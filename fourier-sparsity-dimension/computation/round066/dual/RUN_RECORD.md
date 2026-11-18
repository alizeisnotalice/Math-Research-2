# R63 extended dual-device audit

This run extends the earlier R63 full-core experiment to 140,000,000 new
case IDs. In addition to the full core, it tests \(L=\mathbb F_2^h\times\{0\}\)
with a proper spectral support that spans all of \(L\). Cases use \(h=5,6\),
three external coordinates and restriction dimensions \(1\le k<h\).
The addressed range starts at the earlier run's exclusive upper bound.
Each of the two generators uses deterministically seeded case IDs. A range
can contain repeated maps or functions.

For a holed core, the first selector branch has four negative points on
\(y\), namely \(0,e_1,e_2,e_3\). At a nonzero \(y\)-frequency its Walsh
coefficient is proportional to the sum of the characters of those four
points. That sum is zero exactly when the frequency has Hamming weight
two on the first three coordinates. The zero \(y\)-frequency belongs to
the support as well. The support in \(L\) has size \(5\cdot 2^h/8\) and
contains a basis of \(L\). An even-weight three-variable table in the
second branch keeps all external terms out of the top odd coefficient
layer, which occupies \(2^{h-1}\) points in this holed model. Half of the
extended cases use this holed model.

Both C and Metal compute the complete Walsh spectrum from the restricted
truth table and independently compare every integer coefficient to the
signed fold of the original Walsh coefficients. Both compute the external
rank using the nonzero bands of the restricted transform. They verify
pointwise that the rank loss is at most the number of missing basis bands,
which is at most the number of same-band collisions of their selected
representatives. The eligible number of pairs is bounded by
\(s-|S\cap L|-d_{\mathrm{ext}}\). Case 20000008 has \(|L|=32\),
\(|S\cap L|=20\), \(s=28\), and \(d_{\mathrm{ext}}=3\). In this example
the count of eligible pairs is 3 and \(s-|S\cap L|-d_{\mathrm{ext}}=5\).

| Statistic | Mac mini | MacBook |
| --- | ---: | ---: |
| Case IDs | \([20000000,110000000)\) | \([110000000,160000000)\) |
| Case count | 90,000,000 | 50,000,000 |
| Holed-core case count | 45,000,000 | 25,000,000 |
| Cases with a basis representative collision | 30,475,506 | 16,930,440 |
| Cases with actual external rank loss | 19,242,412 | 10,686,944 |
| Total external rank loss across cases | 37,100,655 | 20,608,216 |
| Number of same-band representative collision pairs | 95,811,526 | 53,236,369 |
| Matching CPU and Metal positional checksum | 3d37b652cf40e438 | 05d2cb64301033c8 |

A separate Python verifier reproduced the positional checksum over the
16,384-case pilot and recorded four explicit controls. It also enumerated
every nontrivial restriction subspace and all 28,740 affine cosets for two
fixed holed controls. Each nonzero core difference collided with probability
\((2^{h-k}-1)/(2^h-1)\) for uniformly selected \(k\)-dimensional subspaces.
For the fixed \(h=6,k=3\) example the expectation of the representative
pair count was \(2/3\) and the expected actual external rank loss computed
from the restricted truth tables was \(59/186\).

## Physical workload and attribution

| Measure | Mac mini | MacBook |
| --- | ---: | ---: |
| Logical CPU workers used | 12 | 15 |
| CPU work interval | 55.020 s | 34.995 s |
| Project CPU time | 285.722 s | 209.718 s |
| Project CPU average | 5.19 cores | 5.99 cores |
| Project CPU sampled median | 5.25 cores | 6.21 cores |
| Project CPU sampled peak | 9.93 cores | 10.29 cores |
| Whole-machine CPU idle, all valid full-worker samples | 0% in 262 samples | 0% in 189 samples |
| Metal device | Apple M6 | Apple M5 Pro |
| Metal completed command buffers | 1,374 | 763 |
| Metal summed command duration | 34.024 s | 178.738 s |
| Metal wall time | 22.783 s | 123.388 s |
| Whole-device GPU idle before launch | 100% | 0% |
| Whole-device GPU utilization during live GPU task, median and peak | 96%, 100% | 100%, 100% |

Whole-machine counters include other workloads. During intervals with no
update to the machine's CPU tick counters, idle samples are unavailable
and excluded from the idle denominator. The project CPU measurements use
that process's own CPU clock and include such intervals. Both hosts reported
positive CPU time for each worker. GPU participation is established by the
project's Metal command receipts. Duration summed over command buffers can
include overlaps on distinct queues. Using the two host clocks, CPU work
and actual Metal commands overlapped across all four resources for 18.654
seconds, in 556 disjoint intervals. Four compute processes exited with
status 0 and a subsequent process check confirmed they were absent.

## Reproduction

Both hosts used Apple clang 21.0.0 and Apple Swift 6.4.2. Python 3.9.6
supervised the compute processes and Python 3.15.0 ran the independent
checks. The verifier requires Python 3.10 or newer.

~~~sh
clang -O3 -std=c11 -pthread restriction_cpu.c -o restriction_cpu
swiftc -O RestrictionMetal.swift -o restriction_metal
./restriction_cpu 20000000 16384 12 0 pilot_cpu.json pilot_cpu.csv
./restriction_metal pilot_gpu.json 20000000 16384 4096
python3 -B verify_restrictions.py local macbook reproduced.json
cmp reproduced.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

To rerun both hosts, compile the same sources on each host, then run
`python3 -B run_host.py output 20000000 90000000` on the Mac mini and
`python3 -B run_host.py output 110000000 50000000` on the MacBook.
Use fresh output directories. Once both supervisors produce
`ready.json`, write an identical future Unix timestamp to the two
`start_epoch.txt` gates. The archived run used 1794978012.763810.
