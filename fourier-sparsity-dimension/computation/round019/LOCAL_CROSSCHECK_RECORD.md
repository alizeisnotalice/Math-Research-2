# R19 local independent popcount/coset crosscheck

This locally compiled checker intentionally does not call the remote CPU
program's FWHT or difference-set implementation. It hardcodes the six cubic
representatives and their orbit sizes, independently rebuilds every degree-3,
degree-2 and linear ANF truth table, and directly computes the weight of every
quadratic-plus-linear phase by 64-bit population count. A weight of 16 or 48
selects the appropriate constant, yielding a weight-16 support.

For every candidate support it searches the list of all 651 two-dimensional
kernels, stopping at the first valid one, using only direct collision tests: for each selected point x and each of the three
nonzero kernel elements h it checks whether `x xor h` is selected. No
difference-set mask is constructed or consulted.

## Actual local run

- Source: `verify_rm36_popcount.cpp`, SHA-256
  `e47e6493cf0863cf0b95be354ac0cf25a0a3ddf6a9a8145b052337d7e36a0e1e`.
- Compiler: Apple clang 21.0.0 (clang-2100.4.1.1), arm64-apple-darwin27.2.0,
  POSIX threads.
- Compile command:
  `/usr/bin/clang++ -O3 -std=c++17 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang verify_rm36_popcount.cpp -o /Volumes/data/CodexTemp/T2/a28/r019/verify_rm36_popcount`.
- Binary SHA-256:
  `acfdca3407383ef646794f56b5568162c5d8d9177bb2603a289cf64acd4d9e16`.
- Run command:
  `/Volumes/data/CodexTemp/T2/a28/r019/verify_rm36_popcount /Volumes/data/CodexTemp/T2/a28/r019/popcount-result.json 12`.
  It returned exit 0. The program starts its timer before its ANF/basis and
  kernel constructions and reports 0.00834533 seconds.

The program has no per-codeword storage: its largest explicit table is 32,768
64-bit quadratic truth tables (256 KiB), plus 651 small kernels and 12 small
worker-local scalars. Thus its explicit dynamic data are far below 64 MiB.
The host was already saturated by unrelated work (start/end top: 0.0% idle;
only 183/147 MiB system memory unused), so no separate peak CPU or RSS claim
for this sub-second process is made.

## Exact agreement

The weighted weight-16 count and weighted projectable count are both
232,081,500; the weighted nongraph count is zero. Fiber counts in orbit order
are `(2604,784,384,448,0,0)`, exactly matching the independent remote scan.
The output SHA-256 is
`3f6752b2a404c148bbcbc94adcc0389ef8aa328f8d5cc626407660c95c1e3904`.

This is an independent finite implementation check, not a proof beyond
RM(3,6).
