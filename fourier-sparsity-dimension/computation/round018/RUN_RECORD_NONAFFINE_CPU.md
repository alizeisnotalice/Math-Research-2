# R18 non-affine four-block family: final exact CPU check

## Implemented finite object

For `r=0,1`, the verifier implements

`Q=x1*x2+x3*x4`,

`F6=x1*x2+x0*(x1+x3)+(x5+Q)*(x2+x4)`, and

`F=F6+sum(u_k*v_k)`.

Blocks fix `x0=a`, `x5=Q+b`, and `v=c`, with `z,u` free. The code constructs
V from all block indicators, the constant, and all coordinate functions; does
GF(2) elimination; enumerates every element of `F+V`; and applies an integer
FWHT. It verifies both pointwise block partitioning and every block's weight.

For dual partial supports it transforms `(-1)^f 1_A` for each block A. At r=0
it checks all V phases. At r=1 it directly checks all global linear phases;
the omitted global constant preserves support, while adding a block indicator
only negates the restricted signal on that block, so also preserves support.
For each direct check it asserts disjointness/coverage, support size `2q`, and
nonzero partial-Walsh magnitude `q`.

## Final result

| r | rank V | all phases | expected abs Walsh | blocks | partial block checks |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 9 | 512 | 8 | 4 × weight 16 | 512 × 4 = 2048 |
| 1 | 14 | 16,384 | 16 | 8 × weight 32 | 256 × 8 = 2048 |

Both final runs have zero bad phases, bad Walsh terms, maximum Walsh deviation,
block partition failures, block weight failures, partial support failures,
partial support-size failures, and partial nonzero-magnitude failures. They
returned exit code zero and a post-run query found no owned process.

## Audit correction

The first source revision incorrectly used `setbit` (OR) rather than a bitwise
XOR when adding a nonzero r=1 global linear phase to F. Its report of 255/256
partial-support failures was therefore an implementation artifact, not a
counterexample. The old source SHA-256 was
`7fd0e5276f0930c46937eea77e8c91c5cb4958587f0fbee818ed559636323071`, and
the invalid preliminary r=1 output remains at
`/Volumes/data/CodexTemp/T2/a28/r018/r1-result.json` (SHA-256
`02810a32ac7d37aaee0bc22228684f4fcaae5fcfc91a1cdcc34c5657504e9758`).

The corrected final source SHA-256 is
`f134e419555bf48f5edb539973a630ae7903106cc72f4503fac9c085845ad0b1`; the
Apple clang 21 binary SHA-256 is
`efd5e9f75d36b9f98ca45836e171cb27b511f4f0c1c3b371463234afbcb01c2a`.
Remote final outputs are r=0 SHA-256
`6da257a0c3571cfc82f88b141f7633989f11a9d141c9f5e9b8d1b26918ee4639` and
r=1 SHA-256
`ea65f06d48a997f812a383fda04bfdcbdbd8c8ce13e9f052183f529ad1d9a477`;
they were retrieved to `/Volumes/data/CodexTemp/T2/a28/r018/`.

## Reproduction receipt

Actual final compilation command:

`CLANG_MODULE_CACHE_PATH=/Users/zhengzhihao/CodexRemoteMath/a28-selector/cache/clang /usr/bin/clang++ -O3 -std=c++17 -pthread verify_nonaffine_blocks.cpp -o verify_nonaffine_blocks_final`

The compiler was Apple clang 21.0.0 (clang-2100.4.1.1), targeting
`arm64-apple-darwin27.3.0`, with POSIX threads.  No Python driver ran the
scan: both were foreground remote-shell commands with 15 C++ threads:

- `./verify_nonaffine_blocks_final 0 15 r0-final.json` — exit 0, program
  elapsed 0.000644917 s.
- `./verify_nonaffine_blocks_final 1 15 r1-final.json` — exit 0, program
  elapsed 0.00917667 s.

The read-only post-run host query found Python 3.9.6, but it was not used as a
driver.  [cpu-receipt.json](cpu-receipt.json) binds the full commands, flags,
versions, thread count, final source/binary hashes, output hashes, elapsed
times and actual exit codes.  Its fields were transcribed from the completed
foreground run receipts and checked against the current read-only remote hash
query; this supplement did not re-execute the phase scans.

Remote work directory: `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r018/nonaffine-cpu`.
The short workloads had only adjacent aggregate host readings, about 28% idle;
these do not show sustained CPU utilization by this verifier.

This is exact finite verification for r=0 and r=1 only, not a proof for the
proposed infinite family.
