# Round 006: common MM splitting and 8-block partitions

## Inputs and exact question

This audit reads the 16 `rows` of `corrected-results.json` and the four rows
of `beam-continuation-4seed.json` from the A28 round005 GPU worktree.  Every
row has `dimV=15`, and is represented by its MM seed phase
\(q(u,v)=u\cdot\pi(v)+h(v)\) and eight relative-sign generators.  The input
SHA-256 values are respectively
`f98be2615e82a7f2df42c8ae51068738d23bedd62636733f683285687afa7b86` and
`36e8778d5378e570bfe9ae7aede5d7d5227d6df7a96444d4f1883950aafa4d05`.

For each of all 1,395 three-dimensional subspaces \(U\leq\mathbb F_2^6\),
the primary audit asks whether \(q\) and every one of the eight generators is
affine on every (U)-coset.  It tests the exact Boolean identity
\(D_aD_bf(x)=0\) for all 64 (x) and the three unordered pairs from a basis
of (U).

## Actual result

- 15 of the 20 recorded certificates have exactly one compatible common
  three-space (U).
- Five do not: corrected-results rows 1, 3, 5, 12 and beam-continuation row
  2 (zero compatible (U) among all 1,395).

Thus these twenty examples do **not** all lie in a common fixed MM split.
This is a finite, certificate-specific finding, not a classification of all
15-dimensional Boolean sign-code spaces.

For every one of the five no-common-(U) rows, the complete 32,768-word
enumeration of (C) plus the eight generators contained exactly 24 masks of
weight 8.  An exact-cover search found eight pairwise-disjoint masks covering
all 64 coordinates.  For the other 15 rows the eight cosets of their unique
compatible (U) supply the partition via equality in the 15-dimensional
affine-on-cosets parameter count.  The stored partitions were independently
checked for weight, disjointness, coverage, and (when found by exact cover)
membership in the 15-dimensional code.

## Verification

The local independent verifier enumerates all 1,395 subspaces by independent
triples rather than the primary incremental span construction, then tests
literal 8-point affine truth tables on every coset.  It agrees exactly with
the primary audit, including all zero and singleton counts.

MacBook performed the same independent verifier through a single SSH stdin
transfer, at `nice -n 10`, using Python 3.9-compatible popcount code.  The
remote receipt hash is
`c56017fa2dd09639de43455a83937762cd9047ab3427371e4b93a366bfbda036` and
the verifier source hash reported by the remote command is
`ac110ec1e768843404541cae1bc895cbcb7204335e84d4801dcf5969f19ca5ff`.

## Remote execution receipt for switch-closure verifier

The subsequently added `verify_switch_closure.py` was also executed on the
MacBook, in a separate one-SSH stdin run, rather than merely having its JSON
compared locally.  The host-stamped
[`remote-switch-closure-receipt.json`](remote-switch-closure-receipt.json)
records hostname `zzhMateBook-Pro-S.local`, exit code 0, Python 3.9.6,
macOS arm64 platform, exact `nice -n 10` command, 0.08701491355895996 seconds
wall time, input hashes, and the complete result.  It reports that all 20
partitions passed membership reconstruction, block affineness, dual Walsh
partition, and joint normal/slope rank 7.

The original main verifier source hash is
`ce2d69ef84c713a78b4c587b59698a7ca7a0b5c6cd9153f72f892bcc9d365627`.
Because the remote default Python lacks `int.bit_count`, the executed adapter
[`run_switch_closure_py39.py`](run_switch_closure_py39.py) makes only two
asserted semantic-preserving substitutions to `bin(...).count('1')`; its hash
is `b8e30e0aeec242665ca0ece4cd44bd86664c964970696f1460b382f3db19dc06`.
The receipt-file hash is separately
`7b3841023e5e3d1c5a4a0daecd89b0257ad5929427f732a11222e212ef977fb8`,
and must not be mistaken for either verifier source hash.

## Artifacts and reproduction

- `audit_common_mm_split.py` — primary all-subspace finite audit.
- `verify_mm_cosets.py` — independent triple/coset-table verifier.
- `find_weight8_partition.py` and `verify_weight8_partitions.py` — complete
  codeword enumeration, exact cover, and no-search partition verification.
- `common-mm-audit.json`, `common-mm-audit.local-verify.json`,
  `common-mm-audit.remote-verify.json`, `weight8-partitions.json`, and
  `weight8-partitions.verify.json` — results and receipts.

```sh
TASK=/Volumes/data/CodexProjects/T2-A28-compute/fourier-sparsity-dimension/computation/round006
TMP=/Volumes/data/CodexTemp/T2/a28/r006
nice -n 10 python3 "$TASK/audit_common_mm_split.py" \
  --input /Volumes/data/CodexProjects/T2-A28/fourier-sparsity-dimension/computation/round005/gpu/corrected-results.json \
  --input /Volumes/data/CodexProjects/T2-A28/fourier-sparsity-dimension/computation/round005/gpu/beam-continuation-4seed.json \
  --out "$TMP/common-mm-audit.json"
nice -n 10 python3 "$TASK/verify_mm_cosets.py" "$TMP/common-mm-audit.json" --out "$TMP/common-mm-audit.verify.json"
nice -n 10 python3 "$TASK/find_weight8_partition.py" "$TMP/common-mm-audit.json" --out "$TMP/weight8-partitions.json"
python3 "$TASK/verify_weight8_partitions.py" "$TMP/weight8-partitions.json" "$TMP/common-mm-audit.json" --out "$TMP/weight8-partitions.verify.json"
```
