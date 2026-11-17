# R18 exact computations and reproduction

All mathematical operations below use integer Walsh sums or GF(2) arithmetic. There is no floating-point tolerance or random seed. Timings are measurements, not mathematical certificates. The infinite family and dimension bounds are proved in `proofs/018-quadratic-graph-blocks-zh.md`.

## Local Python

Python 3.15.0, macOS arm64. `build_graph_cases.py` independently constructs the two codes from the explicit formula, creates exact hexadecimal truth-table masks, and checks all base partial Walsh supports. The r=0 code is completely checked (512 phases); the r=1 code is constructed with rank14 and passed to the independent remote checks. Recorded elapsed time 0.010862 s, exit zero.

`n4_double_switch_audit.py` was rerun by the main agent after adding its source hash: it checks all 65,536 truth tables, all 140 affine two-flats, and every legal weight-eight switch. There are 896 bent phases and 390 weight-eight switches; all390 admit two legal minimal pieces, while only30 supports are affine. Final elapsed time 0.482932 s, exit zero. These statements are for the specified standard four-variable reference phase only.

`nonaffine_four_block_audit.py` is an independent implementation. It constructs the six-variable witness, checks all512 common phases and the four explicit partial-Walsh support formulas, and generates all1,395 linear three-spaces and all11,160 affine three-flats. Every witness block contains zero of all these flats (not merely zero of the568 that are legal for this base phase). Its final run cross-checks the exact input JSON and records source/input hashes. Foreground exit zero; elapsed time is in its JSON.

`verify_rank2_graph_cylinders.py` independently evaluates the epsilon-free cylinder formulas for a rank-two graph example: all64 phases and all1024 parameter/frequency pairs pass. Four cylinders have codimension2, total and projected normal rank3, kernel dimension3, and local restriction ranks2,2. It also checks all651 rank-two alternating forms in dimension6; compatible counts against fixed forms of rank2,4,6 are90,20,0. Elapsed time0.098341 s, exit zero. It uses the archived exact basis/FWHT helper, whose hash is bound.

## MacBook CPU and Metal

The independent C++ source constructs the family from its formula, rather than reading the Python masks. With15 threads it checked all512 and16,384 phases for r=0,1; zero Walsh errors. Each run also checked2048 individual partial-block transforms, correct weights, support sizes2q, nonzero amplitudeq, and disjoint coverage. The r=1 direct partial checks range over all256 global linear phases; block indicators and the constant only negate the respective restricted signals, covering the other support cases analytically. Final program timings are0.000644917 and0.00917667 seconds. These exclude startup and do not establish sustained utilization. Compilation, versions, exact commands, hashes and exit status are in `cpu-receipt.json` and `RUN_RECORD_NONAFFINE_CPU.md`.

Native Metal on Apple M5 Pro reads exact64/256-bit phase masks from `graph-cases.json`. All16,896 phases were checked; independent Swift CPU int64 transforms compared all4,227,072 signed coefficients. The reported38.518124 ms is host submit/wait/readback only, excluding CPU generation and comparison; it is not GPU utilization or a performance comparison. Maximum output buffer is2MiB; no peak-process-memory claim is made. See `RUN_RECORD_GRAPH_METAL.md`, environment and manifest. Both remote tasks completed and their postchecks found no owned process. No MATLAB session was opened.

The CPU preliminary r=1 partial-support failure came from OR instead of XOR when adding a linear phase. The provisional Swift signal133 came from using variable count6/8 in place of truth-table length64/256. The main agent identified both source errors; the corrected implementations were actually rerun. Only final source-bound results are archived here; the diagnoses are retained in the run records.

## Reproduction

Copy sources and inputs to a fresh data-volume directory, such as `/Volumes/data/CodexTemp/T2/a28/reproduce-r018`, to preserve archived results. Set `PYTHONDONTWRITEBYTECODE=1` for Python. From that directory:

```sh
python3 build_graph_cases.py
python3 n4_double_switch_audit.py
python3 nonaffine_four_block_audit.py --input graph-cases.json --output nonaffine-four-block-audit.json
python3 verify_rank2_graph_cylinders.py
```

The C++ source needs C++17 and pthreads; its exact remote command and flags are recorded in `cpu-receipt.json`. Run the executable with `r threads output.json`, where r is0 or1 and threads is1–15. Metal reproduction requires a compatible Mac: place the Swift source and input in the dedicated remote directory named in `run_graph_metal.sh` and run that script, which uses the shared remote compiler cache. Its generated `report.json` is archived here as `graph-metal-report.json`.

`python3 verify_receipts.py` checks archived source/input/output bindings and finite coverage. It does not rerun the phase scans or certify the general mathematical theorem. Source syntax, JSON, shell syntax, file sizes, links and staged Git diff are checked before archival.
