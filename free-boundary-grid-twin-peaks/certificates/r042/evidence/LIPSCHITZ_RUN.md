# Remote exact Lipschitz curved-cut run receipt (r075)

Remote task: `/Users/zhengzhihao/CodexRemoteMath/a022-curved-cut-overlap-pruning-r075`.

## Frozen inputs and exact gate

- Pruned source SHA-256: `eebe4e3d9844851f0d2a69c36d247a4693c3a7d6f65203b6a979bd92d1b3ea94`
- Apple-clang binary SHA-256: `f917681b5f871a25ef7e77012d3b4ad3a5d12d8b953266053f2e5d4b57c91508`
- Lipschitz orbit manifest SHA-256: `57169550ea416bb740067fdda91ceaddbfab58c44833f2803eb7a8d33e0c552b`

The required no-bridge gate compared the pruned program with the historical
unpruned r073-v2 program for curves `22322` and `44344`. Both returned zero.
The pruned program reported `core_bytes=0` and a geometric-zero certificate;
the unpruned comparison allocated its 536,870,912-byte DP table.

## Production and independent audit

All 4,134 symmetry-orbit jobs were newly run; r073 output was not reused. The
remote final receipt reports `completed=4134`, no failures or unlaunched jobs,
and `all_children_waited_reaped=true`. The local audit checked every raw JSON,
per-job receipt and resource-file hash, curve identity, count range, frozen
hash triplet, and manifest multiplicity. It found 178 diagonal and 15,753
off-diagonal original pairs; 1,051 orbit jobs had a strict geometric-zero
certificate and 3,083 had nonzero exact counts.

\[
S_1=16756357278752626754500,
\qquad
S_2=221954860157521456430308.
\]

Using the supplied exact \(A_2(5)=7063261481362958280000\), the finite
second-moment lower bound is

\[
\frac{S_1^2}{A_2(5)S_2}
=\frac{1220763083709591890762320154493866101175}
{6816196584138796738763603366326902311088}
\approx0.17909739964811052400.
\]

Maximum recorded per-job RSS was 538,542,080 bytes. At 2026-11-17T14:23:20Z,
the remote task had no matching process; its driver (99941) and caffeinate
process (99942) were both absent.

This is an exact finite \(n=5\) count and a finite conditional lower bound for
the stated Lipschitz curved-cut event family. It neither supplies a uniform
normalization estimate nor proves a positive-probability infinite subsequence.
