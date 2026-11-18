# R63 restrictions of a full affine spectral core

The generated functions use \(h=5,6\) core coordinates \((t,y)\), three external
coordinates \((z_1,z_2,v)\), and output signs. For \(t=0\), use one negative
point on \(y\), independently of the external variables. For \(t=1\), use
\((-1)^v\) times a three-variable truth table evaluated on three parities
of \((y,z_1,z_2)\). In each case \(L=\mathbb F_2^h\times\{0\}\) is a full
spectral core. Its \(M=2^h\) coefficients are exactly the odd entries of
the integer Walsh transform divided by 16. The case generator selects
dimensions \(1\le k<h\), independent rows \(w_j\) in the core and an affine
shift \(b\). Its tested case IDs are deterministic, and different IDs can
generate the same restriction or function.

Write \(\pi_W(\alpha)=(\alpha\cdot w_j)_j\), where \(W\) is the span of the
rows, and let \(\beta\) denote an unchanged external frequency. For
\(F_{W,b}(u,x)=F(b+\sum_j u_jw_j,x)\), both implementations compare the
integer Walsh coefficient \(A_\mathrm{new}(\gamma,\beta)\) of the actual
restricted truth table against
\[
2^{h-k} A_\mathrm{new}(\gamma,\beta)
 = \sum_{\pi_W(\alpha)=\gamma}(-1)^{\alpha\cdot b}A_\mathrm{old}(\alpha,\beta).
\]
Each implementation computes the true external rank after restriction from
the surviving bands of this separately transformed truth table.

Choose a basis of the original external projection using increasing
\(\beta\), and choose the least \(\alpha\) in each selected band
as its representative. Distinct basis bands retain distinct \(\beta\).
A missing basis band requires a collision of its representative with some
other coefficient of that same band. Thus if \(d\) and \(d'\) are the
original and restricted external ranks, and \(C\) counts collisions of
representatives with individual same-band points, \(d-d'\) is at most
the number of missing basis bands, itself at most \(C\). For uniform
\(k\)-dimensional \(W\), every given nonzero core difference is in
\(W^\perp\) with probability \((2^{h-k}-1)/(2^h-1)\). With
\(D=\sum_{\mathrm{basis}\;\beta}(|S_\beta|-1)\), linearity of expectation
gives \(\mathbb E[d-d']\le \mathbb E[C]=D(2^{h-k}-1)/(2^h-1)\) and
\(D\le s-M-d\). The exact expectation experiment enumerates each \(W\)
in reduced row-echelon form and all its affine cosets.

## Results

The Mac mini verified case IDs \([0,8000000)\), and the MacBook verified
\([8000000,20000000)\). Both hosts ran C and Metal on their own full ranges.
Each method passed the individual checks and independently produced the
same aggregate and positional checksum on that host.

| Measure | Mac mini | MacBook |
| --- | ---: | ---: |
| Cases | 8,000,000 | 12,000,000 |
| Cases with a representative collision | 2,925,982 | 4,389,844 |
| Cases with actual external rank loss | 1,756,841 | 2,636,878 |
| Actual external rank loss summed over cases | 3,606,802 | 5,413,085 |
| Same-band representative collision pairs | 9,724,550 | 14,590,066 |
| Positional checksum | 4eb73ec95ef77cba | 43847353d0ee5c7c |

A separate Python implementation computed both Walsh transforms and the
basis diagnostics for the 8,192-case pilot and reproduced its common checksum.
For the two fixed control models it enumerated all 3,195 nontrivial linear
restriction subspaces and all 28,740 affine cosets of those subspaces.
For \(h=6,k=3\), the expected number of representative collision pairs is
exactly 1, whereas the expected actual external rank loss is \(104/279\).
The four controls include equal subspace directions and different shifts
that give full exterior rank or annihilate the entire exterior spectrum.
The deterministic main ranges are an experimental cross-check.
The exact expectation formulas were checked by full enumeration for the
two named models with full cores.

## Device records

The Mac mini used 12 CPU workers and Apple M6 Metal. The MacBook used 15
CPU workers and Apple M5 Pro Metal. All workers spent positive CPU time in
the enumeration. Project CPU peak was 7.48 and 7.97 logical core
equivalents, and the project CPU averages were 3.56 and 3.77 cores.
The devices finished 489 and 733 Metal command buffers. Their total
recorded GPU command durations were 6.540 and 35.025 seconds. The durations
can include concurrent commands on distinct queues. Using the observed
host clocks, 3.818 seconds of CPU work overlapped reported GPU commands
on both hosts, in 303 positive-length intervals. Whole-device CPU and GPU
utilization reached 100% on both machines. GPU utilization was already
100% before launch on both hosts. The command receipts track the Metal
work attributable to this experiment. All four experimental compute
processes and their supervisors exited.

## Reproduction

The Mac mini used macOS 27.2 (26B123); the MacBook used macOS 27.3
(26C5035r). Both used Apple clang 21.0.0 and Apple Swift 6.4.2. Python 3.9.6
ran the supervisors, and Python 3.15.0 ran the archive and expectation
checks. The verifier uses Python 3.10 or newer.

From this directory on Apple Metal hardware with suitable temporary and
module-cache paths:

~~~sh
clang -O3 -std=c11 -pthread restriction_cpu.c -o restriction_cpu
swiftc -O RestrictionMetal.swift -o restriction_metal
./restriction_cpu 0 8192 12 0 pilot_cpu.json pilot_cpu.csv
./restriction_metal pilot_gpu.json 0 8192 2048
python3 -B verify_restrictions.py local macbook recreated-summary.json
cmp recreated-summary.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

For a synchronized run, copy the sources to both hosts, build the CPU and
Metal binaries and run `python3 -B run_host.py output 0 8000000` on the
Mac mini and `python3 -B run_host.py output 8000000 12000000` on the
MacBook. Use fresh output directories. Once both supervisors write
`ready.json`, write a common future Unix timestamp to
`start_epoch.txt` in both output directories.
