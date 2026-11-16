# CPU forced-second-birth SIS

This is an independent CPU implementation of the sequential importance
sampler for the original free-boundary twin-peaks law.  It is a numerical
diagnostic only.  The code does not certify a conditional probability.

For every forced second-birth rank `tau`, the proposal chooses a root at rank
one, a current boundary vertex before `tau`, an unrevealed nonboundary vertex
at `tau`, and a boundary vertex afterwards.  A branch that has an empty
forced set contributes zero rather than being removed and renormalized.  The
product of forced choice-set sizes is the importance weight.  These are the
exact pointwise identities documented in the project proof, while this C++
program evaluates their finite random-sample means.

`original_law_sis.cpp` uses 12 threads and separately seeded `mt19937_64` instances (stream independence is not certified)
with a rejection-uniform index sampler.  Each batch records A1, its weight
ESS, and per-rank A2, Euclidean-square and original-Manhattan observables.
It saves complete descending orders for an external peak verifier.  The
external `verify_saved_labels.py` reconstructs labels and rechecks strict
local maxima and Manhattan distance without sharing the C++ peak predicate.

The direct arithmetic implementation is deliberately restricted to n<=8.
On this machine, `long double` has 15 decimal digits and max exponent 308;
while W and W^2 are finite through n=8, n=12 requires log-sum-exp or scaled
block reductions.  It exits before running n>8, so it cannot silently create
infinite W^2 diagnostics.

Executed commands used the Apple clang module cache on the data volume:

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang original_law_sis.cpp -o original_law_sis
./original_law_sis 8 1 200000 12 n8-strata
python3 -B summarize.py n8-strata
python3 -B verify_saved_labels.py n8-strata/saved-k2-labels.json
python3 -B audit_results.py
```

`n8-strata` consists of 63 ranks times 200,000 trajectories = 12,600,000
K2-proposal attempts.  A late, near-zero-mass stratum has ESS 8.79, but that
alone is not a valid global degeneration conclusion.  Recombining layers
with their actual A2 mass shares gives global A2-weight ESS 143.93 and a
maximum normalized contribution 0.0672.  The n=7 counterpart has global ESS
4598.08 and maximum contribution 0.00528.  `global_diagnostics.py` records
every layer's share and the exact recombination formula.  No stored
second-moment of W*EuclideanSquare exists, so these numbers do not diagnose
the Q numerator variance.

The two n=8 batches still differ in Q (0.32110 and 0.32917), so both remain
search diagnostics rather than reliable conditional-law estimates; this
conclusion rests on the ratio's unmeasured numerator variance and finite
randomness, not merely the minimum late-rank ESS.  All completed C++
processes exited after their runs; `.stderr` files preserve per-batch timings.

The larger independent `n8-strata-r002` run made 63 * 4 * 200,000 =
50,400,000 attempts in 13.40 seconds wall time (98.14 seconds aggregate user
time; 2.42 MB maximum RSS).  Its global A2-weight ESS is 1075.92 and maximum
normalized contribution 0.01826; its 90% A2 mass lies in ranks 2 through 8.
The individual minimum ESS 22.04 belongs to a low-mass tail.  Its aggregate Q
is 0.32917, whereas the earlier 12.6-million run gave 0.32110.  No claim about
the original asymptotic question follows from either numerical ratio.

归档说明：SOURCE-RUN-HASHES.txt保留原交接/运行哈希及历史路径，主Agent已核对原文件。归档对说明文字作适用范围修正后，当前文件以课题manifest-r010.json为准；不把历史绝对路径当可移植复现入口。
