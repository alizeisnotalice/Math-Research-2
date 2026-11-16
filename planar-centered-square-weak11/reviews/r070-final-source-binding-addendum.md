# R70 final version and recovery binding

Finalized 2026-11-16 in pcw-sol-cover-structure. This addendum binds the
mathematical reading and existing exact results described in
[the total review](r070-independent-mass-radius-prefix-and-area-review.md).
It does not turn matching hashes into a new mathematical verification.

The newly executed byte/provenance checker checked six frozen component
manifests, 51 payload/dependency/source rows, and the round report. The
51 rows include repeated inherited-source bindings and are not a count
of distinct new files. Exact per-path pins are saved in
`computation/results/round070-independent/source-bindings.json`.

| Component | Frozen manifest SHA-256 |
| --- | --- |
| Main author proofs/domain check | 98cfb814f88d818197da28dd50401e3628585a92de56218243f47175b3752132 |
| Three production controls | 4c24d10ac335c441d80ac0ff4f41842a445bbfd5afa80f79054a16f2b9bc90c3 |
| Balanced four-corner positive class | 14f100e984dee7baa2295512a137ede1569d25247440fda27a91493e138fc7a3 |
| Separate all-H/prefix independent check | 8ac6967455e4fe73ecc15338f69cbd0c3d232770928fb65db037d956db3595d3 |
| Root finite receiver-point grid | 05e780669ce7e5cab88e5d61e9d4d18e7efbd155d0fe5f818eff67f71b0c7d6b |
| Root Wolfram scalar checks | 84274868a3a567fcde340e38c698fa3781c70b2bf31c3308f001c2721dad4f6e |

The captured-mass author proof is frozen at
`43054eb9d7ec35b8bf8cd27067fb2360ca4ba13c2de012b7bdb8cc890030989d`;
the four-corner obstruction proof at
`94226fade8d92d444318d3e76e4c24d6e46b5465279af10546b58275bc6f50c8`;
the final report `rounds/070.md` at
`c37e6fe36ee153c5b3d7a0f754135d654f3b416abc4fc849c32a1a3c22d328fd`.
The final audit integration manifest separately pins this review, the
independent grid proof, scripts, inputs and actual results.

The saved original 2026-11-15 independent runs are archived unchanged as
`grid-clamp.json`, `saved-cells.json` and `prefix-n17.json`. Their embedded
source hashes match the final checkers. The leaf replay also matches every
production payload/dependency pin. It independently checked 36414 leaves,
65805 literal nonempty subsets, 88 proper closures, 166551 leaf-label
tests and 2745792 compressed area cells. Balanced actual Delta was
independently matched; balanced true P, auxiliary M-minus-I/old-credit
fields, and historical actual P/Delta have only the source/code or reused
result binding specified in the total review. The n17 replay checks all
289 sources at 256 points, 5186 grouped events and 179 passing centers;
larger n were not independently replayed.

`recovery-provenance.json` is explicitly not a process receipt. The three
old output files record elapsed times and Python 3.13.18; contemporaneous
execution reports recorded successful completion. Separate old PID/exit
receipts were not found, and recovery neither reran the computations nor
re-observed their old process states. Root's persisted receipts are read
as historical assertions: prefix receipt has exit0 and nine PIDs recorded
absent (parent plus eight result-returning workers, twelve requested),
while Wolfram receipt has exit0/PID absent and all 14 saved tests pass.
Root's corrected README describes the closed process pool. The old
root-prefix manifest ending `46705` is superseded by the final pin above.

The new 2026-11-16 Fraction scalar run for proof section 6 has its own
`finite-prefix-area-receipt.json`: Python 3.15.0, actual exit0, owned PID
3903 absent after wait. It verifies the positive common-m prefix margin
`78401/8000020000` and the area certificate `2401/625`, with the input
specified by n=400001 and a uniform-grid formula. Root also independently
read the analytic proof of every bound and the Omega/C_n containment.
No source, prefix or subset enumeration of the 160000800001-source input
was performed. There is no actual halo or deletion-area certificate.

From the integrated topic directory, recheck these byte bindings using
an explicit new output (all source packets must already be integrated):

    python3 -B computation/review_r070_source_bindings.py --spec computation/inputs/r070-source-bindings.json --packet computation/results/round070-compute-integration-manifest.json --packet computation/results/round070-core-integration-manifest.json --packet computation/results/round070-balanced-four-corner-integration-manifest.json --packet computation/results/round070-allH-prefix-independent-manifest.json --packet computation/results/round070-root-prefix-grid/root-prefix-grid-manifest.json --packet computation/results/round070-root-symbolic/root-symbolic-manifest.json --report rounds/070.md --review-root . --archive-root .. --output /Volumes/data/CodexTemp/T2/r070-replay/source-bindings.json

The general C2 problem, actual outside-source compatibility and a general
W/Wprime/deletion payment remain unresolved. The goal stays active; the
earlier general bounds remain unchanged. This independent workspace
creates no commit or push. Root owns publication and remote verification.
