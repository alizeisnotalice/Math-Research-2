# R19 exact computation record

All mathematical arithmetic is integer or GF(2), with zero tolerance. No floating-point approximation supports a dimension claim. The general conditional theorem is analytic; the RM(3,6) conclusion is a finite exhaustive computation with the reduction in `proofs/019b-six-variable-projection-audit-zh.md`.

## Local CPU

1. `macwilliams_rm36.py` (Python 3.15.0) constructs and checks the 22/42-dimensional dual monomial codes, enumerates 32,768 homogeneous quadratic parts representing all 4,194,304 RM(2,6) words, and obtains the full RM(3,6) distribution by exact integer MacWilliams sums. All coefficients are integral/nonnegative, their sum is 2^42, and the weight-16 coefficient is 232,081,500. Actual runtime 1.877525458 s, exit zero; source binding is in `macwilliams-rm36.json`. No random seed is used.
2. `verify_rm36_popcount.cpp` actually ran with 12 threads. It independently reconstructs the six representative fibers, uses population count rather than FWHT, and tests projection by direct point/coset collisions rather than a difference-set mask. The 12,582,912 quadratic/linear combinations give 4,220 weight-16 representative words; all projectable, with weighted total 232,081,500. Actual program runtime 0.00834533 s, exit zero. The timer includes basis/kernel construction; compiler, commands and hashes are in `LOCAL_CROSSCHECK_RECORD.md`. The orbit representatives/sizes are inherited from the separate exhaustive orbit computation, not independently classified by this checker.
3. `d15_graph_projection_audit.py` was rerun by the root after making input handling portable: its R6 audit already contains seed permutation and h fields, so legacy R5 absolute paths are provenance only. Both R6 input files are hash-bound. The final run checks 5,680 weight-16 occurrences, 5,076 distinct supports, and zero nongraph supports in the 20 archived dimension-15 codes. It checks each base bent phase and basis rank; common bentness of the full old codes remains the historical R6 certificate, not a new exhaustive recheck. Four R18 controls each have 16 admissible kernels. The deterministic negative-control search selects seed0, mask `0x1270021c6492004`, with no admissible kernel, and all 651 decisions agree with direct collisions. Exact runtime, Python/platform and source hash are in its final JSON; foreground exit zero.
4. `nonstar_quotient_freeze_local_audit.py` actually used four worker processes and exited zero. Its result records 0.897655709 s for the script; a separate invocation wall observation was about0.935 s. It builds one four-coordinate nonstar constant-rank-two cube and checks the four h cases. It verifies j=(2,1,1,0) and local kernel image dimensions (3,2,2,1). This is an **abstract affine-direction model**: fixed centers are set to zero and unused beta-kernel coordinates are omitted; it does not check actual Walsh centers or construct a global common-bent code. For h=3 the Y space is enlarged to dimension5, giving codimension6; for h=0,1,2, Y has dimension4 and codimensions2,3,4. This supports the normal/kernel bookkeeping only.

The preliminary D15 difference implementation incorrectly added repeated bit positions as integers, creating carries. It was corrected to a deduplicated set/OR mask and actually rerun. The old numerical values are invalid and are not archived as results. The root's portable final version was rerun again after copying. The complete RM scan and independent local checker do not rely on that defective implementation.

## MacBook CPU and GPU

`audit_rm36_weight16.cpp` exhaustively classifies all 2^20 cubic leading masks using three GL generators and checks six complete quadratic/linear/constant fibers. It uses 15 CPU threads for the fiber scan and stops kernel search once a valid kernel is found. Its 0.00991133 s timer starts **after** orbit construction, generator checks and control checks; wrapper wall time1.083731167 s includes subprocess and observation overhead. Neither is a sustained CPU utilization claim. Start/end host observations had about26% idle. Source/runner/binary hashes, compiler21, Python3.9.6, command and exit zero are in `full-run-receipt.json` and `RUN_RECORD_RM36.md`.

`R19.swift` uses native Metal int32 on Apple M5 Pro, with an independent Swift int64 CPU transform. All196,608 phases and12,582,912 signed coefficients agree; fiber counts are2604,784,384,448,0,0. The maximum explicit GPU output buffer is524,288 bytes. It reads the orbit input from the CPU classification, so GPU agreement verifies coefficients/counts, not orbit completeness or projection. The shell command succeeded and a postcheck found no owned process. Elapsed time was not instrumented; Swift6.4.2/macOS27.3 are the surrounding verified environment, not a new per-run version capture. Acceptance is the JSON `status == "passed"`, since this source does not return nonzero solely for a failed final comparison. See `RUN_RECORD_METAL.md` and `metal-result.json`.

## Resource lifecycle and scope

The user's requested local CPU saturation was observed in three consecutive host samples:0.0% idle on the12-logical-CPU host. Process inspection showed unrelated research workloads supplied most of this load. It is not attributed to A28. A28's12-thread check and four-worker local model actually ran, then exited; no ongoing A28 busy loop was left behind. MacBook CPU and Metal checks also completed and postchecks found no owned jobs. No MATLAB session was started. The computation is intentionally compressed; do not inflate a weighted coverage count into a claim of two hundred million individual kernel searches or persistent device saturation.

## Reproduction

Preserve archived outputs: copy the needed sources/inputs to a new directory under `/Volumes/data/CodexTemp/T2/a28/`. For the D15 checker, pass explicit paths to this repository's `computation/round006/weight8-partitions.json` and `common-mm-audit.json`. It no longer opens their legacy absolute `source_file` paths. All scripts require only the Python standard library or Apple C++/Metal toolchains.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 macwilliams_rm36.py
PYTHONDONTWRITEBYTECODE=1 python3 nonstar_quotient_freeze_local_audit.py --workers 4
PYTHONDONTWRITEBYTECODE=1 python3 d15_graph_projection_audit.py --partitions /absolute/checkout/fourier-sparsity-dimension/computation/round006/weight8-partitions.json --audit /absolute/checkout/fourier-sparsity-dimension/computation/round006/common-mm-audit.json
clang++ -O3 -std=c++17 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang audit_rm36_weight16.cpp -o audit_rm36
./audit_rm36 orbit 1 orbits-gpu.json
./audit_rm36 full 12 full-result.json
clang++ -O3 -std=c++17 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang verify_rm36_popcount.cpp -o verify_rm36
./verify_rm36 popcount-result.json 12
```

On the authorized MacBook, use its own existing remote cache path instead of the local data-volume path. Metal compilation/run and hash arguments are recorded in `RUN_RECORD_METAL.md`. The wrapper `run_rm36_full.py` resolves sources, binary and outputs relative to its own directory and uses15 threads; the direct executable commands above allow choosing a local thread count.

`verify_receipts.py` checks archived bindings, finite counts, independently reconstructed monomial inputs and the generator action on every basis monomial. It does not rerun the full scans or prove the general theorem. Root source inspection, two mathematical reviews, JSON/source/link checks and staged diff inspection complete this round's archival checks.
