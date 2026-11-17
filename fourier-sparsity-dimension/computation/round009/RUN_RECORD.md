# Round 9 execution and reproduction

Date: 2026-11-17 (Asia/Shanghai). Baseline 8d6847d3f0297493a3d6df82ce935a6af6750bb7.
All calculation uses exact Python integers/F2 or Fraction; zero tolerance.
The universal claims are proved in proofs/009-partial-affine-rank-zh.md and
proofs/009b-residual-budget-zh.md. Finite checks are implementation evidence.

## Recursive normal-rank certificates

The final source was strengthened after parent review to validate inputs,
rank-preserving deletion order, inclusion-minimal core, true intersections,
restricted right sides, extension disjointness, both Kraft budgets, private
dimension sums, recursion dimensions and full certificate consistency.
The n=4 check uses separately evaluated point masks to cross-check equation
consistency and extension disjointness/measure. This is recursive certificate
replay, not a second independent mathematical proof algorithm.

The local run used Python 3.15.0, nice 10, and took 8.367185125 seconds.
The full result is local-recursive-audit.json. Coverage:

- All 171 affine flats of codimension 0, 1 or 2 in F2^4 and all 11072 pairwise
  disjoint subfamilies including empty. Counts by size 0..4 are
  1,171,3165,6230,1505. Maximum normal rank 3.
- Every exhaustive case was recursively replayed, with a rolling digest of
  its full certificate; representative certificates are archived.
- 120 fixed-seed mixed-depth independent-query-tree leaf subfamilies in F2^31,
  seed 2026111709, depth at most five. Their complete inputs and certificates
  are archived; no 2^31 point enumeration was needed.
- Ten full/deleted-one-leaf tree controls, depths 1..5, attain ranks
  1,3,7,15,31. This includes proper partial packings attaining the bound.
- Sixteen round008 packing witnesses (14 jobs and two pilots), all rank seven.
  These old witnesses were actually rechecked this round; the old searches
  were not rerun or counted as new searches.

## MacBook execution and cross-machine comparison

The existing macbook-math alias was authenticated with strict host-key
checking. Hostname zzhMateBook-Pro-S.local and Python 3.9.6 were actually
checked; the new task directory was absent before dispatch:
/Users/zhengzhihao/CodexRemoteMath/a28-selector/r009/recursive-audit.
No configuration, credentials, network, firewall or sleep settings changed.

run_remote_audit.py sends the exact verifier and archived r008 input by
base64 inside SSH standard input. The remote task refuses an existing task
directory, writes its PID/command, runs nice 10 Python, saves completion,
then removes its running marker. The final verifier is Python 3.9 compatible.
Actual remote verifier time: 12.261471708 seconds. Total SSH call:
13.342897875 seconds. Both calculation and SSH exited zero. A subsequent
pgrep for the precise remote verifier path returned exit 1 with no output.
All launched calculation processes are finished.

compare_recursive_audits.py checks every certificate section across the two
machines, the input and source hashes, and the exhaustive rolling digest.
Only Python/time and the machine-specific input pathname are excluded.
The exact data agree; digest:
6e36b4cc10245763f7c55678946c8c86fb01367b31c08574774dde133ab21845.
The comparison result is cross-machine-verification.json.

This round used one remote replay process, not a new full-load search batch.
The round008 full-load evidence remains correctly labelled round008. Once the
partial-rank theorem ruled out the search target analytically, further
rank-eight searches would not advance the research. No CPU/GPU busy-work was
launched. No MATLAB, Wolfram or GPU context was started this round.

## Residual weight/quantization checks

check_residual_gap.py actually ran locally with Python 3.15.0, nice 10, in
0.266788291 seconds. It enumerated all 2048 RM(2,4) codewords (positive weights
4,6,8,10,12,16); all 256 Boolean three-variable truth tables (smallest possible
maximum absolute normalized Fourier coefficient 1/2); and all 65536 Boolean
four-variable truth tables, verifying algebraic degree two for the 896 bent
ones and even weights on 7168 coordinate three-faces. The result is
residual-gap-verification.json. The general RM second-weight formula is a
cited theorem, not inferred from this n=4 computation.

## Reproduction

Run from the repository root, with data mounted. Use a Python interpreter
already installed on the machine. Local environment was
/Users/zzh/.local/bin/python3; the scripts require only the standard library.
Temporary outputs below stay on the data volume.

```sh
nice -n 10 python3 fourier-sparsity-dimension/computation/round009/recursive_normal_rank_certificate.py --remote-packing fourier-sparsity-dimension/computation/round008/remote-packing-execution.json --out /Volumes/data/CodexTemp/T2/a28/r009/replay-recursive.json
nice -n 10 python3 fourier-sparsity-dimension/computation/round009/check_residual_gap.py --out /Volumes/data/CodexTemp/T2/a28/r009/replay-residual-gap.json
python3 fourier-sparsity-dimension/computation/round009/compare_recursive_audits.py --out /Volumes/data/CodexTemp/T2/a28/r009/replay-cross-machine.json
```

The remote launcher is a record of the exact dispatched job and deliberately
refuses to overwrite its existing directory. Before any new remote replay,
inspect that task's status and select a fresh task identifier; do not remove
old records or restart based solely on an observation timeout.

Only self-authored proof/code/results are archived. No external full text,
compiled binaries, caches, private SSH configuration or credentials are
included. Every archived file is below the project 10 MiB default. GitHub's
official size policy was checked in the immediately preceding round on the
same date; no file approaches its limits.
