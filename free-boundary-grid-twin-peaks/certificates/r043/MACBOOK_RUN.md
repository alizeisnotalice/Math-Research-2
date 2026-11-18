# Exact n=7 curved single-region root-band run (r076)

The frozen source accepts only the two reviewed domains: n=5 with widths
`{2,3}` and n=7 with widths `{3,4}`. It computes the two complementary
one-region root-band counts separately and uses a 256-bit product only for the
full event count.

## Gates and provenance

- Source SHA-256: `8b3426ae8069bd19d0cf04df2daf871f780e50c7b66374a7f2463ba6ca0ba87b`
- Remote Apple-clang binary SHA-256:
  `16385f3223f3f28d42ff0b4c87ea24313b51f7b54e868346d13f5069319c6124`
- Both local and remote fresh n=5 gates matched all 32 reference profiles.
  The remote gate receipt is `remote/gate/n5-all/audit.json`.
- The exact full-grid normalization is read from
  `computation/results/n7-exact-r005.json` in the public project, SHA-256
  `ff7974ff8487a339516fee95324674efd873abd9a36b8fc50b0f8121c632f8d6`,
  with
  \(A_2(7)=137077928808481869759712347602605837601908319657172480\).

## Production and audit

All 128 profiles in `{3,4}^7` completed on the remote MacBook, with no failed
or unlaunched profile and with every child waited and reaped. The local audit
checked every raw output, receipt and resource-file hash; it also checked the
complement-pair equality \(C_h=C_{7-h}\), the reviewed region/independent-set
bounds (maximum 28 vertices and alpha 14), and the exact normalization-file
hash.

The exact sum of event counts is

\[
S_1=162492250696126049348416111751504718204933980525701248.
\]

For \(X=\sum_h1_{E_h}\), one has \(0\le X\le128\), so the only union bound
asserted here without overlap information is

\[
\mathbb P_2(X>0)\ge\frac{S_1}{128A_2(7)}
=\frac{9003338358606274897407807610344897950184728530901}
{972183892258736664962498919167417287956796593313280}
\approx0.00926094171102572611.
\]

The quotient \(S_1/A_2(7)\) is the finite expectation \(\mathbb E_2X\), not
a union probability and is not reported as one. This is still a finite n=7
certificate; it does not establish an infinite-subsequence result.

At 2026-11-17T14:40:31Z there were no matching task processes; the controller
(9188) and its caffeinate process (9189) were absent. The maximum per-profile
RSS recorded in the receipts was 20,561,920 bytes.
