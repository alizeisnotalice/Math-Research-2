# R69 finite-valued means with dependent common labels

The physical coordinates are an address \(y\), kernel coordinates \(v\),
four common quotient coordinates \(u\), and five private coordinates \(z\).
Let \(g(y,v,u)=\mathbb E_z F(y,v,u,z)\). In frequency, the common
labels are those at zero private label. Contracting the kernel frequency
coordinates maps the common labels to at most five affine-independent
labels on \(u\). The support of the actual common spectrum is obtained by
taking the full \((v,u)\) transform for each address, then taking the
union across addresses. Kernel fibers can carry many different
frequencies of \(v\). Generated case IDs are deterministic and can
sometimes yield the same model.

For fixed \(y,v\), write \(t\) for the actual number of nonzero
coefficients in the \(u\)-spectrum of \(g\). This includes the zero label
if its coefficient is nonzero. Actual Boolean private slices determine
\(g\in\{0,\pm1/2,\pm1\}\). The R68 cost envelope \(c\) is \(9/4\), \(4/3\),
or \(0\) when \(|g|\) is \(0\), \(1/2\), or \(1\), respectively. The
experiment verifies
\[
t\le4,\qquad t+\mathbb E_u c(g)\le529/96 .
\]
It also computes the cost from the real private spectrum, including
zero-mean characters whose cost is 1. The implementation checks the cost
as an exact integer in units of \(1/192\). For \(t=4\), the absolute means
\(0,1/2,1\) occur exactly \(6,8,2\) times, giving \(192t+6\cdot27+
8\cdot16=1058\).

For affine-independent quotient labels, multiplying \(g\) by an additional
independent uniform sign makes its active characters jointly uniform signs.
The sign sums in \(4g\) are even integers of magnitude at most 4. Flipping
one sign forces the corresponding coefficient of \(4g\) to be an integer.
Choosing all signs to align with their coefficients bounds the sum of
coefficient magnitudes by 4. Thus \(t=4\) forces each coefficient magnitude
to be 1 in \(4g\) and gives the stated mean counts. The upper bound for
\(t\le3\) is \(3+9/4\).

Each tested private slice is an actual Boolean truth table. The cases use a
constant when \(|g|=1\), a signed quadratic phase on the first two private
bits when \(|g|=1/2\), and either a character or a signed quadratic phase
on the next two bits multiplied by the fifth character when \(g=0\).
The union of private labels at each address lies in the union of the
plane through zero and the translated plane in the respective independent
coordinate blocks. Both implementations check, on the actual support union,
the R68 condition that each four-point affine plane has a direction
represented by exactly two unordered pairs.

| Statistic | Mac mini | MacBook |
| --- | ---: | ---: |
| Case-ID range | \([0,6000000)\) | \([6000000,10000000)\) |
| Models with affine-dependent common spectrum | 5,999,917 | 3,999,952 |
| Fixed \((y,v)\) pairs | 168,000,056 | 111,999,904 |
| Cost-envelope equality pairs | 47,375,215 | 31,585,106 |
| Actual private-spectrum cost equality pairs | 46,032,196 | 30,690,649 |
| Maximum common-support affine excess | 71 | 71 |
| Maximum number of frequencies in a kernel fiber | 16 | 16 |
| Matching CPU and Metal checksum | 2503d8feacb7c08b | 47a9ef52067a01ff |

C computes each private Walsh support by parity-mask popcounts. Metal
independently transforms the real private truth table by butterflies.
The Python checker reproduces the 8,192-case pilot with 229,368 fixed
\((y,v)\) pairs and independently recomputes the affine ranks.
Model 0 has 80 common frequencies in five fibers, affine rank 8, and exact
actual cost \(529/96\) at every one of its 64 fixed \((y,v)\) pairs.

Two further exact Python controls use respectively eight and ten kernel
bits. Multiplying each of four quotient characters by a nondegenerate
quadratic sign on \(v\), possibly modulated by a kernel character, gives
four complete kernel fibers. Their actual common supports have 1,024
and 4,096 points with respective affine ranks 11 and 13. Both independently
give actual cost \(529/96\) at every kernel input.

## Device receipts

The CPU runs used 12 and 15 workers on the Mac mini and MacBook.
All workers accrued positive CPU time. Average project work used 7.52
and 12.77 logical core equivalents with respective peaks 8.27 and
13.21. Whole-machine CPU idle had median 0% on both hosts, with 43 of 43
and 13 of 14 valid full-worker samples recording idle 0%. The machines
completed 2,930 and 1,954 commands on Apple M6 and Apple M5 Pro. The
respective sums of Metal command durations were 29.752 and 49.139 seconds.
Two command queues on each device permit overlapping intervals.
Both GPUs had 100% median whole-device utilization, starting from
baselines of 9% and 100%. Whole-device telemetry includes other processes.
The command receipts establish this project's GPU activity. The common
overlap of both CPU execution windows and both sets of GPU command
intervals was 2.368 seconds in 161 intervals, as measured by the host
clocks. All four compute processes exited with status 0 and were absent
when their PIDs were subsequently checked.

## Reproduction

The two hosts used macOS 27.2 (26B123) and 27.3 (26C5035r), Apple clang
21.0.0 and Apple Swift 6.4.2. Python 3.9.6 supervised the main run;
Python 3.15.0 ran the independent verifier. The verifier requires Python
3.10 or newer.

~~~sh
clang -O3 -std=c11 -pthread fiber_cpu.c -o fiber_cpu
swiftc -O FiberMetal.swift -o fiber_metal
./fiber_cpu 0 8192 12 0 pilot_cpu.json pilot_cpu.csv
./fiber_metal pilot_gpu.json 0 8192 2048
python3 -B verify_fibers.py local macbook reproduced.json
cmp reproduced.json VERIFIED.json
shasum -a 256 -c SHA256SUMS
~~~

After compiling on both hosts, run
`python3 -B run_host.py output 0 6000000` on the Mac mini and
`python3 -B run_host.py output 6000000 4000000` on the MacBook.
Use fresh output directories, and after both supervisors write
`ready.json`, put the same future Unix time in their
`start_epoch.txt` gates. The archived start time is 1794979244.387525.
