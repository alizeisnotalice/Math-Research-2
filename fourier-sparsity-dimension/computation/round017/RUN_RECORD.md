# R17 actual computation and reproduction

All algebra and Walsh transforms use exact integers/GF(2); there is no numerical tolerance or random seed. These runs check finite instances. The arbitrary-dimensional result is proved in `proofs/017-double-block-rank-budget-zh.md`.

## Local exact enumeration

`verify_rank_two_families.py` uses Python 3.15.0 on macOS 27.2 arm64. The complete four-dimensional run enumerated all 2,825 linear subspaces of the six-coordinate alternating-form space and all 26,387 distinct affine cosets. Accepted family counts by dimension are 35, 595, 840, 105. It also independently reconstructs the rank-two set as all nonzero wedges.

The five- and six-dimensional runs enumerate every affine constant-rank-two family containing the fixed form `e0 wedge e1`, by exhaustively extending direction subspaces inside the set of allowed differences. Counts by direction dimension are respectively `(1,154,616,192)` and `(1,650,3000,1200,48)`. Only this anchored scope is claimed. All 48 dimension-four six-variable families have a common factor. Elapsed times from the source monotonic clock were 0.339889, 0.079200 and 4.395483 seconds. All three foreground commands exited zero.

`verify_nonparallel_blocks.py` constructs phases on bits `(x,y,z,t,w,u)` with base `xy+zw+tu`. Its first partition is `w=0,u=c` or `w=1,t=c`. The code includes global affines, all four block indicators, and `xz,yz`; its dimension is 11. The second case omits `xz,yz` and splits the first block further by `x`, obtaining three size-16 and two size-8 blocks and dimension 10. The script verifies affine block structure, distinct directions, exact full bentness, disjoint covering of branch frequency supports, one/four intrinsic spectral points, and intrinsic quadratic rank zero/two. All 3,072 phases and 13,312 branch instances passed in 0.950781 seconds, foreground exit zero.

## Independent MacBook CPU enumeration

`enumerate_affine_rank2_4x4.py` was separately implemented and run with `/usr/bin/python3` 3.9.6 on the MacBook. All 26,387 affine spaces were covered, agreeing with the local counts. A dimension-six star family in ambient dimension eight and an embedded rank-four-difference line were also checked. Final elapsed time 0.293137 seconds; final exit zero; no owned process remained. The initial implementation's overly strong test on the rank-four-difference branch was corrected before this final result. See `RUN_RECORD_AFFINE_RANK2.md` for the failure boundary and final hashes; `remote-enumeration-result.json` binds the actual final source.

## Native Metal and independent signed CPU comparison

`NonparallelMetal.swift` reads the exact full truth-table masks in `nonparallel-cases.json`; there is no hard-coded reference phase. Apple M5 Pro native Metal evaluated all 3,072 phases. Swift CPU int64 independently computed all 64 Walsh coefficients for every phase, giving 196,608 exact signed comparisons, all equal; every amplitude is 8. Basis ranks and the supplied point-set partitions were also checked. The Swift block check does not check local Walsh supports; the Python check above does.

Swift 6.4.2, macOS 27.3, Metal execution width 32; maximum output buffer 524,288 bytes. All command buffers completed without error. The measured 4.390501 ms is the sum of host submit/wait/readback intervals, not pure kernel time or GPU utilization. Compiler startup and CPU comparison are not included in that number. The foreground runner completed and the remote postcheck found no verifier or Swift compiler process. See `RUN_RECORD_NONPARALLEL_METAL.md`, environment, manifest, and report for bindings. The compiled binary is regenerable and is not archived.

## Independent joint-cylinder audit

`check_cylinders.py` imports no project helpers. It reconstructs local ANFs, the non-common-factor branches, the joint freezing kernel, all eight cylinders in each example, and their normal spaces. It checks all 512×64 and 1024×64 parameter/frequency pairs against direct integer partial Walsh sums. The first case has b=4, r=2, dim V0=9, normal and projected ranks 4, and kernel dimension 5. The second has b=r=0, dim V0=10, normal and projected ranks 5, and kernel dimension 5. Local and joint restriction bounds also pass.

Source review found an ordinary integer sum where GF(2) XOR was required in the provisional implementation. The corrected source, with additional explicit assertions and self/input hashes, was rerun; only that v2 output is archived. Final monotonic elapsed time 12.059141 seconds, foreground exit zero, and the observed process was absent after completion. See `RUN_RECORD_CYLINDERS.md`. The provisional output is not used as evidence.

## Resources and scope

At 18:22:40 and 18:22:42 Asia/Shanghai, the local 12-core host was at 0% idle, principally occupied by other projects; the load average was 57.88. No other project's process was changed. R17's local computations are short single-process exact checks; they are not claimed to saturate the host. The remote CPU enumeration was also short and single-process. All tasks launched for this round finished, and no MATLAB session was opened. MacBook remains available for the main large computations when the research requires them.

## Reproduction

Preserve archived outputs by copying these sources to a fresh directory on the data volume, for example `/Volumes/data/CodexTemp/T2/a28/reproduce-r017`. Use `PYTHONDONTWRITEBYTECODE=1` for the Python commands. With this directory as the current directory:

```sh
python3 verify_rank_two_families.py --n 4 --all-cosets --out local-n4.json
python3 verify_rank_two_families.py --n 5 --out local-n5-anchored.json
python3 verify_rank_two_families.py --n 6 --out local-n6-anchored.json
python3 verify_nonparallel_blocks.py
python3 enumerate_affine_rank2_4x4.py
python3 check_cylinders.py --input nonparallel-cases.json --output cylinder-audit.json
```

The independent enumerator writes `enumeration-result.json`; the archived original is named `remote-enumeration-result.json`. For native Metal, place the Swift source, shell runner and exact input JSON in the remote task directory named in `run_nonparallel_metal.sh`, then run that shell script. It uses the remote shared compiler cache. Reproducing local-only checks does not require Metal.

`python3 verify_receipts.py` verifies the archived source/input/result bindings and coverage counts, without rerunning the finite scans or certifying the general proof. JSON formatting, Python syntax, shell syntax and Git diff checks are also checked before archival.
