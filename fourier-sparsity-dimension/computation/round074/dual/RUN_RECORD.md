# R74 exact Boolean multiway recursion

A Boolean sign function \(F(t,z)\) has \(k=1,2\) independent address
coordinates \(t\) and five or six data coordinates \(z\). Splitting into
\(B=2^k\) children gives the actual Boolean truth tables
\(F_b(z)=F(b,z)\). Generated children are characters, quadratic phases
on independent forms, point flips on an intersection of parities, or
arbitrary small Boolean truth tables lifted along independent linear
forms. Individual leaves have signed modulations. Case IDs can generate
the same model.

Let \(g_\beta(b)\) be the unnormalized Walsh coefficient of child \(b\)
at data label \(\beta\). Write \(a_\beta=|\operatorname{supp}g_\beta|\)
and \(m_\beta=|\operatorname{supp}\widehat g_\beta|\), where the latter
transform is over the addresses. The full integer parent coefficients
must satisfy
\[
A(\alpha,\beta)=\sum_b(-1)^{\alpha\cdot b}g_\beta(b)
\]
at every frequency, and this equation is checked against a full Walsh
transform of the original truth table. The finite uncertainty inequality
\(a_\beta m_\beta\ge B\) is tested for each nonempty band. Cauchy–Schwarz
gives \(\|g_\beta\|_\infty^2\le m_\beta\mathbb E_b|g_\beta(b)|^2\),
and the middle quantity is at most
\(m_\beta(a_\beta/B)\|g_\beta\|_\infty^2\). In an exclusive band,
\(a_\beta=1\), every address frequency is present and the displayed
formula specifies its sign.

Let \(C\) be the binary span of the labels that are repeated across
children. Each band outside \(C\) belongs to a unique child. If \(e_b\)
counts the labels of child \(b\) outside \(C\), the parent tail has
exactly \(q=B\sum_b e_b\) frequency points. Exclusive labels in \(C\)
can occur and are included with the parent core. If \(r_C=\dim C\),
\(r_b\) is the tail quotient rank of child \(b\), and \(r_T\) is the
true rank of the union of the child tail images, then the full data rank
is exactly \(r_C+r_T\) and \(r_T\le\sum_b r_b\).

The complete parent rank \(d\) is also computed through band gluing,
by taking one representative per band and all within-band address
differences. Its address contribution is between zero and \(k\).
The presence of any exclusive band forces all \(k\) address directions
into the parent span. If \(d_C\) denotes the rank of the parent core
spectrum and the parent tail is nonempty, the exact increment is
\(d-d_C=r_T+k-(d_C-r_C)\). If the tail is empty, \(d=d_C\).

| Statistic | Mac mini | MacBook |
| --- | ---: | ---: |
| Case-ID range | \([0,12000000)\) | \([12000000,20000000)\) |
| True parent bands | 169,588,244 | 113,063,848 |
| Uncertainty equality bands | 151,733,577 | 101,165,768 |
| Exclusive bands | 133,360,744 | 88,919,718 |
| Exclusive bands inside \(C\) | 37,033,264 | 24,675,413 |
| Exclusive bands outside \(C\) | 96,327,480 | 64,244,305 |
| Parent tail frequency points | 276,980,416 | 184,751,524 |
| Cases with \(r_T<\sum_b r_b\) | 5,045,578 | 3,365,622 |
| Matching CPU and Metal checksum | ebaa9aadab737935 | a6a078dd9926f9bb |

C uses direct parity-popcount evaluation of the child Walsh transform.
Metal recomputes the child spectra by butterflies. The full parent
transform is computed independently of the child spectra by both.
The Python checker reproduces the 8,192-case pilot, including all
115,708 true parent bands and six explicit signed truth tables. It also
checks both two-level binary recursion orders for the four-leaf controls.

In control 1, the child character labels \(4,4,1,5\) leave a
one-dimensional repeated core and two tail labels. The two tail labels
have the same quotient, so the actual tail rank is 1 while the child ranks
sum to 2. In control 3, the repeated labels 1 and 2 span \(C\).
The labels 0 and 3 are exclusive and lie in \(C\); the single exclusive
outside label 4 accounts for all four parent tail points. Control 4 has
no exclusive label and realizes one missing address dimension.

## Devices

| Measure | Mac mini | MacBook |
| --- | ---: | ---: |
| Project CPU workers | 12 | 15 |
| CPU workload wall time | 9.106 s | 1.878 s |
| Project CPU average logical cores | 4.72 | 13.94 |
| Project CPU peak logical cores | 9.12 | 15.00 |
| Whole-machine idle 0% samples | 48 of 48 | 13 of 14 |
| Metal device | Apple M6 | Apple M5 Pro |
| Completed Metal commands | 1,465 | 977 |
| Sum of Metal command durations | 4.065 s | 1.179 s |
| Whole-device GPU utilization median / peak | 63% / 89% | 61% / 100% |
| GPU utilization before launch | 7% | 0% |

The CPU idle denominator includes system tick samples with changes.
Whole-device GPU utilization can include other work. The project's
Metal executions are confirmed by the completed command receipts. Across
both host clocks, the actual four-resource overlap was 0.484 seconds
over 395 disjoint intervals. All four compute processes exited with
status 0 and a subsequent PID check found none.

## Reproduction

The hosts used macOS 27.2 (26B123) and 27.3 (26C5035r), Apple clang
21.0.0, and Apple Swift 6.4.2. Python 3.9.6 supervised both runs.
Python 3.15.0 ran the independent verifier, which requires 3.10 or above.

~~~sh
clang -O3 -std=c11 -pthread recursion_cpu.c -o recursion_cpu
swiftc -O RecursionMetal.swift -o recursion_metal
./recursion_cpu 0 8192 12 0 pilot_cpu.json pilot_cpu.csv
./recursion_metal pilot_gpu.json 0 8192 4096
python3 -B verify_recursion.py local macbook reproduced.json
cmp reproduced.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

To rerun the main job, compile the same sources on both hosts and use
`python3 -B run_host.py output 0 12000000` on the Mac mini and
`python3 -B run_host.py output 12000000 8000000` on the MacBook.
Once both supervisors write `ready.json` in fresh output directories,
write the same future Unix time to both `start_epoch.txt` gates.
The archived start time was 1794982266.510526.
