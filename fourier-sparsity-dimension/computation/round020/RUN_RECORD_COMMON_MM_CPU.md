# R20 common-coordinate MM pairs: final MacBook CPU receipt

The final finite sample uses the fixed seed `20261120`. For every
`k=2,3,4,5,6` and each of `c4_t0`, `c2_t1`, `c0_t2`, it constructs eight
independent cases, for 120 total. The corrected source randomizes the relative
offset `r=h+j` on every changed fiber (four values in c4, two in c2); the
unchanged full-fiber pattern supplies t=0,1,2 respectively.

For each case the program computes actual truth tables for f and g, runs an
integer FWHT on each to test bentness, checks `|supp(f+g)|=2q`, constructs the
case-specific `(k-1)`-dimensional kernel from the proof recipe, and checks
directly that no support point collides with a translate by a nonzero kernel
vector. This is a finite crosscheck, not a proof for all permutations/offsets.

## Actual final run

- Host/work directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r020/common-mm-cpu`.
- Compile command:
  `CLANG_MODULE_CACHE_PATH=/Users/zhengzhihao/CodexRemoteMath/a28-selector/cache/clang /usr/bin/clang++ -O3 -std=c++17 -pthread verify_common_mm_pairs.cpp -o verify_common_mm_pairs_final`.
- Source SHA-256:
  `6e416b99f3f816985105d3fdfe69e827971f4aabb375a4f5d33d71f5ee74a150`.
  Binary SHA-256:
  `5fc4b918a5b708315c59487195f9c05bb1b29818331f4cab15577fb3cfa8b2d2`.
  Compiler: Apple clang 21.0.0 (clang-2100.4.1.1), arm64-apple-darwin27.3.0.
- Smoke command `./verify_common_mm_pairs_final smoke smoke-final.json`:
  15 cases, exit 0, output SHA-256
  `59fc50eb5785d8d7630e125f590d4fd34a65426929ac89f7c1441907da42d523`.
- Final command `./verify_common_mm_pairs_final full cases-gpu-final.json`:
  120 cases, exit 0, program elapsed 0.00695383 seconds, all checks passed.

The final GPU input is [cases-gpu-final.json](cases-gpu-final.json), SHA-256
`5e7303916cc3029c7088955e5bed677071878cf717ed4578f078c46d02641789`.
It explicitly records `n=2k`, `N=2^(2k)`, and the little-endian word indexing
`x+(y<<k)`. It was compared against root's separately produced
`local-final-cases.json`: after deleting only `elapsed_seconds`, every field
of every one of the 120 cases is equal. This file is distinct from and did not
overwrite root's local output.

The standalone portable metadata is [cpu-final-receipt.json](cpu-final-receipt.json).
The post-run remote process query found no owned R20 process. The task was too
short to claim sustained CPU utilization.
