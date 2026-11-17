# Fixed-hyperplane exact neighbour enumeration

Date: 2026-11-17, Asia/Shanghai. This is a finite specified neighbourhood,
not an asymptotic A28 conclusion.

## Scope and gate

The input is round005 `corrected-results.json`, row 1. Its sign code is
`C + span(g0,...,g7)`, with `dim C=7` and total dimension 15. The selected
hyperplane is `W=C+span(gj: j!=2)`, dimension 14. Complete enumeration of
its minimum words and exact cover confirms that W has 14 weight-8 words
and contains no eight-block partition. Among all 255 C-containing hyperplanes
of the original code, 192 pass this no-partition gate.

The first proposed deletion `j=0` fails the gate: it contains another such
partition and is already confined by the analytic closure theorem. An earlier
debug scan of that deletion also had an incorrectly transcribed H mask and
used the absolute phase in an H-membership statistic. It is rejected and not
used for any mathematical conclusion here. The integrated source uses an
input header generated from the JSON certificates, checks relative phases,
and deletes generator 2. A malformed JSON quote in the unexecuted intermediate
corrected source was fixed before the accepted remote run.

## Accepted remote run

`run_remote_neighbor.py` transferred the C++ source/header through one strict
SSH session to the authorized MacBook. It compiled with Apple clang, ran a
1,048,576-candidate pilot, then launched 14 `nice 10` processes on disjoint
intervals covering all 268,435,456 elements of RM(3,6)/W. Actual process IDs,
intervals, commands, compiler version, source hashes, exits and timings are
in `remote-neighbor-execution.json`.

The bent predicate is an exact integer 64-point Walsh transform. For every
candidate, all 128 slices modulo C must pass. Early rejection shortens a
candidate test; the candidate ranges remain exhaustive. Gray traversal visits
each quotient element once. All 14 processes and the SSH session exited 0.
The worker computation times were 0.632–0.660 seconds each; the full orchestration
including CPU sampling took 1.345 seconds, and the SSH/compile/pilot/full call
took 2.441 seconds. These are recorded elapsed times, not sustained throughput
estimates or a CPU benchmark.

The 15-core MacBook's first whole-machine `top` sample during the batch showed
0.0% idle. A28 process-CPU sampling peaked at 952.9% (about 9.53 cores); other
work and system activity used the remaining capacity. The next top sample
was already after the short worker phase and showed 61.38% idle. This batch
reached full machine utilization; it did not keep the machine busy after the
useful work finished.

## Independent local exhaustive scan

`independent_derivative_scan.cpp` uses all 63 nonzero derivatives being
balanced as the bent criterion, reverses the monomial order when constructing
the quotient complement, and reverses the slice order. Two local nice-10
processes independently covered all 268,435,456 candidates in 1.517 and
1.707 seconds respectively. Each executable first compared its derivative
criterion against integer Walsh on 4,096 deterministic cubic phases.
Both processes exited 0. Compiler/source metadata are in
`local-execution-record.json`.

The two different complements give different representatives but identical
classes modulo W: zero and the class of g2. `verify_neighbor_results.py`
checks the complete disjoint range ledger, source hashes, quotient equivalence,
and all 256 returned coset slices independently with Python integers.
Consequently the only 15-dimensional extension of this W inside RM(3,6)
is the original V. No global optimality is claimed.

## Reproduction

Compile locally using the existing compiler and data cache:

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -O3 -std=c++20 independent_derivative_scan.cpp -o /Volumes/data/CodexTemp/T2/a28/r006/independent_derivative_scan
nice -n 10 /Volumes/data/CodexTemp/T2/a28/r006/independent_derivative_scan 0 134217728 /Volumes/data/CodexTemp/T2/a28/r006/local-derivative-0
nice -n 10 /Volumes/data/CodexTemp/T2/a28/r006/independent_derivative_scan 134217728 134217728 /Volumes/data/CodexTemp/T2/a28/r006/local-derivative-1
```

The remote runner deliberately refuses to overwrite an existing completed or
running manifest. Before a replay, inspect that manifest and use a new task
subdirectory; do not delete a marker and restart a possibly live job. The
archived receipts are sufficient for the independent verifier and do not
require a remote replay.
