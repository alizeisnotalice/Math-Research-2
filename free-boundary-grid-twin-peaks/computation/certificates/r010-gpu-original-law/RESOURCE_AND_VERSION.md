# Direct-original-law SIS handoff resource record

- Host: macOS Metal default device, as emitted in each JSONL metadata row.
- Compiler: the host `swiftc -O`; Swift/Clang module cache was explicitly
  `/Volumes/data/CodexTemp/cache/swift-modules`.
- Arithmetic: Metal `float` log weights, host `Double` log-sum-exp;
  deterministic `mix32` then `xorshift32` streams.  This is not ideal random
  arithmetic, so the identity proofs describe the target algorithm rather than
  certify floating-point samples as unbiased.
- Direct-basic runs: r002, r003.  Old `run-r001` is deliberately excluded:
  it contains pre-fix invalid JSON (`nan`) and is not a source of results.
- Stratified pilot: `tau-r001`, seed domain `0xA022D***`.
- Stratified independent productions: `adaptive-n7-r001` and
  `adaptive-n8-r001`, seed domain `0xA022E***`.  The allocation files mark
  `pilot_not_combined: true`; production totals exclude all pilot paths.
- GPU processes were checked absent after both adaptive results were complete.
- Exact anchors: tree audits give \(A_1(2)=16,A_2(2)=8\),
  \(A_1(3)=29568,A_2(3)=143808\), and
  \(Q_3=1100/2247\) (Euclidean-square convention).  The production n=7
  reference is encoded by its executable and report, not used for n=8.

`sha256sums.txt` is generated from the package after explicit file selection.
All JSONL raw records and summaries are numerical evidence only; no probability
bound or asymptotic statement follows from them.
