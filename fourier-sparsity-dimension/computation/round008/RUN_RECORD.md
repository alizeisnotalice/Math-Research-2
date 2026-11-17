# Round 8 execution and reproduction

Date: 2026-11-17 (Asia/Shanghai). All mathematical checks use integers,
binary linear algebra, or exact Fraction thresholds. No floating tolerance.

## Triangular low-weight hierarchy

check_triangular_formulas.py exhausts all masks in the standard triangular
families for k=1,2,3. It checks the minimum-word classification, its span,
intersection with the affine code C, and low-weight span dimensions at
integer and fractional thresholds. The integrated run writes
triangular-small-results.json.

For k=4, triangular_weight_scan.cpp exhausts all 2^31=2,147,483,648 masks.
Each is a 256-bit truth table, traversed in Gray order. It records a complete
weight histogram and the OR of active parameter coordinates at each weight.
For any weight threshold this OR supplies an upper bound on the dimension of
the low-weight span; the individually eligible coordinate basis vectors give
the matching lower bound. It does not pretend that a set of individually
legal directions is automatically jointly legal: the whole triangular family
is proved legal by its explicit triangular permutation construction.

The initial MacBook attempt was not dispatched: SSH timed out before
connection. Read-only rediscovery initially found no usable route. After the
user again requested full MacBook load, the existing alias became reachable.
The hostname and 15 logical CPUs were checked, and the new task directory was
confirmed absent before dispatch. No network, credentials, firewall or sleep
settings were changed. Both the initial failure and subsequent recovery are
distinguished in remote-connectivity-status.json.

While the connection was unavailable, run_local_triangular.py compiled the C++ source on the local
machine and launched eight nice-10 processes on disjoint ranges. All exited
0; the complete worker phase took 1.414 seconds. The actual compiler, source
and binary hashes, process IDs, ranges, commands, outputs and versions are
in local-triangular-execution.json. No claim of measured continuous full CPU
utilization is made for this short local batch. After connectivity recovered,
run_remote_triangular.py performed the same full 2^31 enumeration with 14
MacBook workers; it exited 0 and its complete histogram was independently
verified. That batch was too short for process-CPU sampling to measure a
useful peak; it is not the evidence for sustained full load.

verify_triangular_scan.py independently enumerates all 32768 triangular
permutations and counts the intercept choices combinatorially. If a slope
field is nonzero on a fibers, the histogram contribution at weight
q*a/2+q*b is 2^a binomial(q-a,b). The complete resulting histogram agrees
with the 2^31-mask scan. At q=16 the low spans have dimensions 24,28,30,31
for L=1,2,4,8. After adding C their residual dimensions are 4,1,0,0.

## Partial-coverage interface

verify_partial_coverage_interface.py uses the seven round007 dimension-15
extensions. Five have four E4 components and two have two E4 plus four
isolated components, giving 5*16+2*64=208 component subsets.

For every subset it verifies the slope-kernel construction, the exact
dimension identity, the supported residual code, and its atom bound.
The parent strengthened the worker version to check the nullspace equations,
membership in V of every distinct residual zero-extension, and the residual
Walsh amplitudes as well as their supports. The archived result comes from
this strengthened integrated run.

For all nonempty residuals it checks the seed support size, every coarsest
atom's uncertainty inequality, and containment of its Walsh support in D_R.
All K-slices are explicitly transformed for the 32 singleton-component
subsets only: 229376 slices. Other subsets use the exact linear-algebra and
atom checks, and are not reported as full K-slice brute-force checks.

Selecting three of four E4 components leaves a three-dimensional residual
code, not a one-dimensional code. This is a positive control for the extra
"all minimum components selected" hypothesis in the near-coverage theorem.

## Partial affine packing and explicit obstructions

affine_normal_search.py completely enumerated all codimension-two flats in
F2^4: 35 normal spaces, 140 flats, and disjoint-family counts
140,2730,6020,1505 for sizes 1..4. Maximum normal rank is 3. A direct proof
in proofs/008c-partial-codimension-two-zh.md establishes the same bound for
arbitrary ambient dimension, independently of this finite computation.

Bounded codimension-three searches used 1000 trials per seed at n=8,9
(seeds 2026111801/02) and 250 punctured-partition trials (seed 2026111803).
They found no rank-8 example. This is not a proof of the conjectured
general partial-packing bound. Total worker run time was 22.153 seconds.
The parent replayed the independent stored-witness checker. An initial
noncanonical-basis pilot duplicated normal spaces; the canonical RREF code
was corrected before the archived accepted run. No rejected pilot is used.

verify_interface_obstructions.py independently checks two explicit n=4
examples: a legal two-dimensional code whose union-support flip is illegal,
and a partial one-block cover whose slope rank is 2 rather than the invalid
full-cover bound 0. It saves complete input masks and the failing Walsh
vector in interface-obstructions.json.

## Sustained 14-process MacBook search

At the user's explicit request to use the MacBook fully, packing_search.cpp
ran a separate, substantive search for the unresolved partial codimension-3
packing cases. It starts from punctured depth-three decision-tree partitions
with seven independent query normals, then uses exact-validity add/remove/
replace moves, GL(3,2) row changes, and bounded simulated annealing. The
acceptance heuristic uses floating arithmetic; every reported flat and rank
certificate is checked with exact integer/F2 operations.

run_remote_packing.py first measured one million iterations at n=8 and n=9,
then launched 14 nice-10 jobs with distinct seeds. Seven jobs used 98,927,842
iterations at n=8; seven used 59,439,156 at n=9. Total accepted run length:
1,108,568,986 search iterations, not an exhaustive enumeration of families.
Each process ran 118.860–121.238 seconds. Orchestration took 121.736 seconds;
the SSH/compile/pilot/full call took 125.382 seconds.

A28 CPU sampling peaked at 1391.1%; many later samples were approximately
1350–1380%. Whole-machine top samples recorded 0.0% and 0.6% idle. A separate
live status check during the batch recorded 0.0% then 0.18% idle. This was a
sustained useful compute batch, not an idle stress test. All fourteen workers
and SSH exited 0 after their fixed budgets completed.

No rank-8 packing was found. verify_remote_packing.py independently rebuilds
all fourteen final certificates and both pilots from normal rows and right
sides, checks every point mask, cardinality, pairwise disjointness, rank,
source/runner hash, command parameters, unique seeds and completion ledger.
The exact reconstruction passed. Search failure does not prove the general
rank bound; only the separate low-codimension theorems are mathematical
upper bounds.

## Two-hole boundary and cross-machine comparison

compare_triangular_scans.py compares all 257 histogram and parameter-activity
bins between the complete local and remote executions; both agree exactly.
verify_rm_twohole_boundary.py uses the explicit ANF monomial masks 14,21,35.
Its weight-16 support contains none of all 11160 affine three-flats in F2^6.
The complement contains 940 eligible flats; a complete exact-cover search of
14441 states proves there is no six-flat cover. Therefore this is only a
counterexample to a degree/weight-only splitting argument, not to completion
of actual six-flat packings. All checks use integer arithmetic.

## Commands

Run from the repository root. Output directories below must be on the mounted
data volume; full-run directories deliberately refuse existing task markers.

```sh
python3 fourier-sparsity-dimension/computation/round008/check_triangular_formulas.py --out /Volumes/data/CodexTemp/T2/a28/r008/replay-small.json
python3 fourier-sparsity-dimension/computation/round008/verify_triangular_scan.py --receipt fourier-sparsity-dimension/computation/round008/local-triangular-execution.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-k4.json
python3 fourier-sparsity-dimension/computation/round008/verify_triangular_scan.py --receipt fourier-sparsity-dimension/computation/round008/remote-triangular-execution.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-remote-k4.json
python3 fourier-sparsity-dimension/computation/round008/verify_remote_packing.py fourier-sparsity-dimension/computation/round008/remote-packing-execution.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-packing.json
nice -n 10 python3 fourier-sparsity-dimension/computation/round008/verify_partial_coverage_interface.py --candidates fourier-sparsity-dimension/computation/round007/quotient-candidates.json --extensions fourier-sparsity-dimension/computation/round007/all-maximal-extensions.json --audit fourier-sparsity-dimension/computation/round006/common-mm-audit.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-interface.json
python3 fourier-sparsity-dimension/computation/round008/verify_interface_obstructions.py --out /Volumes/data/CodexTemp/T2/a28/r008/replay-obstructions.json
python3 fourier-sparsity-dimension/computation/round008/verify_affine_normal_result.py fourier-sparsity-dimension/computation/round008/r008-search.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-affine-witnesses.json
nice -n 10 python3 fourier-sparsity-dimension/computation/round008/affine_normal_search.py --trials 1000 --out /Volumes/data/CodexTemp/T2/a28/r008/replay-affine-search.json
python3 fourier-sparsity-dimension/computation/round008/compare_triangular_scans.py fourier-sparsity-dimension/computation/round008/local-triangular-execution.json fourier-sparsity-dimension/computation/round008/remote-triangular-execution.json --out /Volumes/data/CodexTemp/T2/a28/r008/replay-cross-machine.json
python3 fourier-sparsity-dimension/computation/round008/verify_rm_twohole_boundary.py --out /Volumes/data/CodexTemp/T2/a28/r008/replay-twohole.json
```

A full local 2^31 replay can use run_local_triangular.py with a fresh --work
directory and an explicit --out path. The remote helper must not be rerun
after an ambiguous observation failure without checking its unique task
directory and live processes. The initial failure occurred before SSH
connection; after recovery, both remote task directories were inspected and
the successful jobs were run once, to completion.

All launched local and remote calculation processes exited. This round did not start
MATLAB, Wolfram or GPU contexts. The workload is exact bitset/linear-algebra
work suited to CPU, and earlier GPU-produced inputs are still used by the
partial-coverage checks. External papers and compiled binaries are not
included in the public commit.

Archive size policy was checked against the GitHub official page on
2026-11-17: https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github
All archived files are individually below the project default of 10 MiB.
