# r014 larger-grid packed-mask GPU audit

This isolated package checks the two-peak SMC GPU transition on the first
sizes that use three, four, and five 32-bit occupancy words.  It is a kernel
correctness test, not statistical evidence about the estimator.

## Actual result

`gpu_large_audit` ran 30 mode/case combinations on Apple M6 and returned zero
failures.  For every case, each rank was checked against an independent Swift
CPU implementation, then against the same GPU observer path submitted either
in one command buffer or one rank per command buffer.  The exact rankwise
comparisons are: all allowed counts, group sums, CDFs, totals, selected
parents, dead status and every `State(m0,...,m4,root1,root2)`.  The terminal
total, `Z`, and the implementation's terminal `U` were checked too.

The cases are `n=9,10,12`, each with `tau` equal to `2`, `floor(N/2)`,
`N-1`, `N`, and `N+1`, with 17 particles, in both multinomial and
integer-stratified modes.  `tau=N+1` is the one-peak choice rule: the audit
checks all transitions and `Z`; the retained two-peak terminal `U` has no
one-peak interpretation.  Successful `n=12` records have `final_m4=65535`,
which explicitly covers vertices 128 through 143 in the fifth mask word.

The production-style host function was also invoked once for every mode/case
and matched the instrumented single-command result in `dead`, `Z`, and `U`.
In addition, the separately compiled, uninstrumented `production_r013.swift`
CLI actually ran `n=9,10,12`, each with `tau=2`, `P=17`, one
integer-stratified replicate and seed-domain `1`.  These three terminals
exactly match the corresponding audit stratified records (the audit assigns
the stratified branch its `+1` seed suffix); see
`original-terminal-check.json`.  Its full-GPU `respropFull` kernel contains
neither observers nor an audit `parentOut` argument.  The instrumented audit
Metal source does contain copy observers and an audit-only `parentOut`
argument, so it is deliberately different from production.

## Reproduce

```sh
mkdir -p /Volumes/data/CodexTemp/T2/a022-large-grid-kernel-audit-r014
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules \
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  audit.swift -o gpu_large_audit
./gpu_large_audit > /Volumes/data/CodexTemp/T2/a022-large-grid-kernel-audit-r014/replay.stdout 2> /Volumes/data/CodexTemp/T2/a022-large-grid-kernel-audit-r014/replay.stderr
jq '{device,cases,cpu_single_stepwise_comparisons,production_terminal_comparisons,failures}' /Volumes/data/CodexTemp/T2/a022-large-grid-kernel-audit-r014/audit-results.json
```

Observed summary: 30 cases, 90 full CPU/single/per-rank comparison sets, 30
production-host terminal comparisons, and zero failures.  `audit.stderr` is
empty.  The audit executable exited after the run.

## SHA-256

Production host source: `cdcd733e41bae122666b03c675123963b7a5ae9a26643dfd458493d156063523`.
Its extracted Metal string is `7474e22aee24049eb8692db789f32dff7b46e5c1c9d31f569622a7278c47da39`.
The audit source is `f3629acd8a73513e3ac90b4fceee298dfce6e004ea5c7b17cbe28eb4faccb061`,
and its different observer-bearing Metal string is
`dc328554d52c9b134bc773f1410b5a5d7aec3ae652708a2520930d6559f5809e`.
The result JSON is `1df5375c541bda43163e277058fe7c8319d11d1d5857e11cb32a208741b91a8c`.
These two Metal hashes use exactly the byte string returned by Swift-style
`split('let metal = """', 1)[1].split('"""', 1)[0]`, including the leading
and trailing newline.  The earlier `9dd4...` value omitted the leading newline
and is not a canonical source hash.

The frozen audit source writes its JSON to the absolute temporary path shown above. The uninstrumented production source is archived one directory up as `../production_r013.swift`. Its three CLI invocations use `n 17 1 2 2 OUTPUT full strat 1`.
