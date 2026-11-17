# Round 7 exact computations

Date: 2026-11-17, Asia/Shanghai. All certificates use integers, F2 linear
algebra or Python Fraction, with zero floating-point tolerance. CPU timings
and utilization are measurements, not mathematical certificates.

## Scope and complete scans

Input: round006/common-mm-audit.json, row 1, derived from round005 GPU
corrected-results.json row 1. The header retains the actual eight original
relative generators. This round keeps indices 4,5,6,7, so
W=C+span(g4,g5,g6,g7) has dimension 11. RM(3,6) has dimension 42; the
quotient has dimension 31. The two leading comments of rm3_multineighbor.cpp
retain old wording about a hyperplane and seven generators; the executed
arrays, KEEP_FROM=4, dimension assertions and receipt specify four retained
generators. The source is archived byte-for-byte as actually run.

MacBook: run_remote_multineighbor.py sent the source and header using one
strict-host-key SSH session, compiled, checked a 2^20 pilot, and ran fourteen
nice-10 processes on disjoint intervals covering all 2^31=2,147,483,648
quotient elements. All processes and SSH exited 0. For each candidate the
16 W/C slices must all be bent; failures stop that candidate early without
skipping later candidates. The exact predicate is integer Walsh modulus 8.
The source/header hashes, actual compiler, Python version, commands, PIDs,
range ledger and outputs are in remote-multineighbor-execution.json.

Local: independent_multineighbor.cpp uses the equivalent criterion that all
63 nonzero derivatives are balanced, a reverse monomial complement and a
different slice order. Each executable cross-checks this predicate with Walsh
on 4096 deterministic cubic phases (seed 2026111706). Two processes cover
2^30 candidates each and exited 0, taking 16.3581 and 18.9007 seconds.
The source and binary metadata were collected at archival and are labelled
accordingly in local-execution-record.json; they are not a retrospective
claim that a pre-run binary hash was recorded.

Both complete scans returned exactly 464 classes modulo W. The independent
Python verifier checks all source hashes in the remote receipt, disjoint full
ranges and hit counts, local/remote class equality, seed/header consistency,
and all 7424 returned coset slices. The coarsest canonical representative set
K spans dimension 18 but has period space zero. K is not itself a vector
space, so its span is not a legal extension.

## Exhaustive finite subspace classification

search_quotient_subspaces.py uses the unique increasing least-outside-span
basis of each binary subspace. At each node the feasible set is the set of
v with v+U contained in K. A target dimension 5 space requires at least 32
feasible elements throughout. This complete DFS finds none (40 visited
nodes, 979 extension candidates).

classify_extensions.py independently uses deduplicated spans rather than
the canonical basis DFS. It retains feasible sets of size at least 16;
this cannot remove a subspace that can reach dimension 4. The counts at
dimensions 0..4 are 1,79,215,101,7. Each of the seven final dimension-4
subspaces has feasible set size exactly 16, independently excluding any
dimension-5 extension. Thus there are exactly seven dimension-15 extensions
of this W in RM(3,6), and none of dimension 16. The filename
all-maximal-extensions.json denotes maximum-dimension extensions in this
specific neighbourhood; it does not classify lower-dimensional spaces that
may be maximal by inclusion.

All seven have an eight-block minimum-switch partition. A separate verifier
checks membership and the round6 analytic closure certificate, along with
minimum-word components and genuine triangle completion from a partial code.
The parent strengthened the worker verifier to assert every full-slice result
and explicitly check all four half-flat atoms, their common direction, the
fixed fourth coset and all six E4 edges; the integrated rerun passed in 6.267
seconds. A graph-triangle XOR check alone was not treated as a certificate of
the fourth-coset identity.
See verify_extension_geometry.py and extension-geometry-verification.json.

## General local-support budget checks

check_weighted_local_budget.py runs on all 2299 stored n<=4 spectrum groups
from round004/remote-result.json. It verifies 62747 supported subspaces,
including 60448 proper supports and 59392 nonflat proper supports. All 95533
evaluation atoms satisfy their exact uncertainty inequalities and vanish
outside the common physical support; 41281 attain equality. These finite
checks support but do not replace the proof for arbitrary Boolean functions.

## CPU accounting and completion

The remote worker times were 4.90926–4.96131 seconds; orchestration including
sampling took 5.0122 seconds, and SSH/compile/pilot/full execution 6.2575
seconds. Process-CPU sampling peaked at 1399.6%, approximately fourteen cores.
The whole-machine top samples reported 0.0% and 5.7% idle. These establish a
short full-load batch, not continued full utilization after completion.
At the user's later status check, no A28 scan process remained; whole-machine
idle was 19.75%. Other work is not attributed to A28 or terminated.

This round's exact bitset predicates suit CPU execution. Round5 GPU-produced
inputs were used substantively. No idle GPU/MATLAB process was started to
inflate utilization. External paper PDFs and compiled binaries remain in
CodexTemp, outside the public commit.

## Reproduction

From the repository root (outputs go to the existing data volume):

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fourier-sparsity-dimension/computation/round007/verify_multineighbor.py --audit fourier-sparsity-dimension/computation/round006/common-mm-audit.json --results fourier-sparsity-dimension/computation/round007 --out /Volumes/data/CodexTemp/T2/a28/r007/replay-verification.json --candidates /Volumes/data/CodexTemp/T2/a28/r007/replay-candidates.json
PYTHONDONTWRITEBYTECODE=1 python3 fourier-sparsity-dimension/computation/round007/search_quotient_subspaces.py /Volumes/data/CodexTemp/T2/a28/r007/replay-candidates.json --target 5 --out /Volumes/data/CodexTemp/T2/a28/r007/replay-dim5.json
PYTHONDONTWRITEBYTECODE=1 python3 fourier-sparsity-dimension/computation/round007/classify_extensions.py /Volumes/data/CodexTemp/T2/a28/r007/replay-candidates.json --out /Volumes/data/CodexTemp/T2/a28/r007/replay-extensions.json
PYTHONDONTWRITEBYTECODE=1 python3 fourier-sparsity-dimension/computation/round007/verify_extension_geometry.py --candidates fourier-sparsity-dimension/computation/round007/quotient-candidates.json --extensions fourier-sparsity-dimension/computation/round007/all-maximal-extensions.json --audit fourier-sparsity-dimension/computation/round006/common-mm-audit.json --out /Volumes/data/CodexTemp/T2/a28/r007/replay-geometry.json
PYTHONDONTWRITEBYTECODE=1 python3 fourier-sparsity-dimension/computation/round007/check_weighted_local_budget.py fourier-sparsity-dimension/computation/round004/remote-result.json --out /Volumes/data/CodexTemp/T2/a28/r007/replay-budget.json
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -O3 -std=c++20 fourier-sparsity-dimension/computation/round007/independent_multineighbor.cpp -o /Volumes/data/CodexTemp/T2/a28/r007/replay-search
nice -n 10 /Volumes/data/CodexTemp/T2/a28/r007/replay-search 0 1073741824 /Volumes/data/CodexTemp/T2/a28/r007/replay-local-0
nice -n 10 /Volumes/data/CodexTemp/T2/a28/r007/replay-search 1073741824 1073741824 /Volumes/data/CodexTemp/T2/a28/r007/replay-local-1
```

The remote runner deliberately refuses to overwrite its completed/running
manifest. A fresh remote full replay requires a distinct task directory after
checking existing jobs; no marker should be deleted to force a rerun.
