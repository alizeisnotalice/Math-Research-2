# R19 RM(3,6) weight-16 support audit

## Scope and exact finite reduction

The computation covers every Boolean polynomial of degree at most three on
six variables only up to an invertible linear change of variables.  It first
enumerates all `2^20` homogeneous cubic coefficient masks under the action on
the degree-three quotient.  The coordinate cycle, the swap `x0,x1`, and the
transvection `x0 -> x0+x1` generate this action: the first two generate all
coordinate permutations and their conjugates give all elementary
transvections; elementary row reduction generates `GL(6,F2)`.

The exhaustive lower-degree fiber is then retained for each orbit
representative: all 32,768 quadratic masks, every one of the 64 linear
shifts, and the uniquely appropriate constant when the Walsh coefficient has
absolute value 32.  Lower-degree terms produced by a linear substitution of a
cubic are therefore absorbed by this complete fiber, rather than discarded as
an assumption.

For every resulting weight-16 support A, the code forms the difference set
with a 64-bit **bitwise OR** mask (including zero).  It tests all 651 two-
dimensional F2 kernels K using `K ∩ (A+A) = {0}`.  A deterministic positive
and deterministic negative 16-point control each cross-check every kernel's
difference-mask decision against direct coset-collision enumeration.  No
nongraph witness occurred, so no 651-kernel witness list was emitted.

## Bounded orbit preflight

The orbit-only run returned exit 0 after validating `cycle^6=swap^2=
transvection^2=id` on every one of the 1,048,576 cubic masks and checking the
651-kernel/control construction.  It found six orbits, with sizes
`1, 1395, 54684, 468720, 357120, 166656`, summing to `2^20`.

The portable GPU handoff includes the orbit representatives' 20-bit cubic
coefficient masks, 64-bit truth masks, sizes, and the explicit 15 quadratic
truth masks/basis order:
`orbits-gpu.json`, SHA-256
`feba42b40fa4185cba8062ddbf662221ab889418028fc8775bc1681fb0a3eab7`.

## Full result

The 15-thread C++ run covered all six representative fibers, each with 32,768
quadratic masks.  Its per-representative counts were:

| cubic representative | orbit size | weight-16 fiber count | nongraph fiber count |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 2,604 | 0 |
| 1 | 1,395 | 784 | 0 |
| 20 | 54,684 | 384 | 0 |
| 1,064 | 468,720 | 448 | 0 |
| 1,536 | 357,120 | 0 | 0 |
| 12,640 | 166,656 | 0 | 0 |

The weighted total is exactly `232,081,500`, matching the independently
provided MacWilliams computation, which is a hard machine assertion in the
final C++ result.  Every counted support was projectable in this finite scan;
the weighted nongraph count is zero.

## Run receipt

- Remote host/work directory: `zzhMateBook-Pro-S.local`,
  `/Users/zhengzhihao/CodexRemoteMath/a28-selector/r019/rm36`.
- Final source SHA-256:
  `fb6ac1006f19240afbb8e938ae53516574813f4d2740b911fcfab0451f76c3c4`.
  Apple clang 21 binary SHA-256:
  `d2d508247cbd39c792ff31439c0512805ffde72a2dd3c9b84311f6ba0f722a50`.
  Python wrapper SHA-256:
  `5ac122eb6170dfc99ee54d28f233d4405e3cbd273d892502c5680321fa735782`.
- Actual command: `audit_rm36_weight16 full 15 full-result.json`; the
  foreground wrapper waited for its owned child, which exited 0.
- C++ elapsed time was 0.00991133 seconds; wrapper wall time was 1.083731167
  seconds.  It was too short for an interval CPU sample: start/end aggregate
  host samples were about 26% idle, not evidence of sustained A28 use.
- Compact result: `/Volumes/data/CodexTemp/T2/a28/r019/full-result.json`,
  SHA-256 `9b590e790327d04187fe9bae289cd3b05a47299819f5b23dd13cf4690ab3473c`.
  Run receipt: `full-run-receipt.json`, SHA-256
  `f8832e0f360efc2b427da1366b1686eb6842f726a35c758ebb0c285e0296772c`.
  A post-run remote process query found no owned runner or scan process.

Two early bounded CLI checks terminated before any usable orbit JSON because
of argument-count/path indexing errors.  They were corrected before the valid
orbit preflight and full run; they are not computation results.  This finite
enumeration is evidence for RM(3,6), not a proof of a general projection
property in higher dimensions.

### Timing and kernel-loop scope correction

The `0.00991133` C++ elapsed field in `full-result.json` starts **after** the
remote executable has built the GL generator maps, validated their identities,
enumerated all cubic orbits, constructed the 651 kernels, and executed the
fixed control collision checks. It measures only the final
representative-fiber quadratic/linear/kernel loop. The orbit-only JSON's own
elapsed field likewise starts after its orbit construction, so neither field
is a full end-to-end BFS time. This clarification changes neither the
preserved remote source hash nor the completed output hashes.

For each actual weight-16 support, the remote kernel loop scans the fixed list
of 651 kernels only until it finds one collision-free kernel, then stops. The
reported `kernel_tests_examined` values are therefore early-stop totals, not
651 multiplied by the number of supports. The all-651 direct collision check
was used for fixed positive/negative controls, and would be serialized as a
complete witness list only if a nongraph support occurred.
