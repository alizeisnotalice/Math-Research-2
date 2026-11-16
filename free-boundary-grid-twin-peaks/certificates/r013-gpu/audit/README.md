# r013 full-GPU correctness audit (actual run)

This directory is an isolated, non-repository audit of the r013 original-law
two-peak SMC Metal implementation.  It was run on **Apple M6** with Apple
Swift 6.4.2.  It is numerical correctness evidence, not a proof about the
SMC estimator or its statistical properties.

## Result

`audit-results.json` records 32 actual runs and reports **zero failures**.
For every run, the audit compared, rank by rank:

1. an independent CPU implementation of the specified SplitMix64/rejection
   and keyed-RNG law;
2. a single-command-buffer GPU execution of the newly instrumented
   `respropFull`; and
3. the same GPU kernels submitted one rank per command buffer.

The exact checks at every rank were the allowed-count vector, 256-wide group
sums, group CDF, total, selected parent vector, dead flag, and every particle
state (`m0...m4, root1, root2`).  At termination the audit also compared total,
`Z`, and `U`.

Coverage is both integer-multinomial and integer-stratified modes for: `n=2`,
every `tau=2...4`, `P=513`; `n=3`, every `tau=2...9`, `P=513`; and `n=7`,
`tau in {2,8,47,48}`, `P=17`.  The `P=513` runs cross the 256-particle group
boundary.  Natural-death cases occur in the requested ranges; additionally an
explicit dead-at-rank-3 case (`n=3,tau=5,P=513`) covered the `dead[0]` branch.
The forced-dead hook is audit-only and is never present in the production
path.

The audit makes the only functional `respropFull` instrumentation requested:
the extra `parentOut` write (and its buffer binding).  Copy-only observer
kernels snapshot GPU-produced values after a rank; they do not feed selection
or mutation.  The audit also compiled and ran the unmodified copied r013
source.  Its 30 normal requested mode/case combinations exactly matched the
instrumented single-buffer run's `dead`, `Z`, and `U`.  The unmodified program
does not expose a terminal total, so that total is compared between the
instrumented single-buffer and per-rank actual GPU paths.

## Files

- `tau_smc_r013.swift`: verbatim copied original source.
- `r013_gpu_rank_audit_prior.swift`: active final audit program, based on the
  supplied parent-instrumentation extension and extended with CPU reference,
  snapshots, dead-branch coverage, and checks.
- `gpu_full_audit`: compiled audit executable.
- `audit-results.json`: machine-readable result records for all 32 runs.
- `tau_smc_r013_original` and `original-*.jsonl`: independent unmodified
  source runs used for the 30 terminal-output checks.

## Reproduction command

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  /Volumes/data/CodexTemp/T2/a022-gpu-full-audit-final/r013_gpu_rank_audit_prior.swift \
  -o /Volumes/data/CodexTemp/T2/a022-gpu-full-audit-final/gpu_full_audit
/Volumes/data/CodexTemp/T2/a022-gpu-full-audit-final/gpu_full_audit
jq '{device,cases,cpu_single_stepwise_comparisons,production_terminal_comparisons,failures}' \
  /Volumes/data/CodexTemp/T2/a022-gpu-full-audit-final/audit-results.json
```

Observed summary: 32 mode/case runs; 96 full CPU/single/per-rank comparison
sets; 30 production-style terminal comparisons; 0 failures.  The separate
unmodified-r013 terminal comparison also passed all 30 exact `dead/Z/U`
records.

## SHA-256

```text
c07f8f091713fa9e939cd53ba26d219ec4a53207499313a63de4015010295590  tau_smc_r013.swift
593e7b104cd355b043a1e9bd4076a2763006d69a91b0061ab50653932997d768  supplied extension-r014/r013_gpu_rank_audit.swift
1b76ff68680abe299eadc4f252c7057ece1491720d9368ab0228b3cba8eeeed5  r013_gpu_rank_audit_prior.swift
ea2b7e463ff066a00bc2afa3c0a2839cfbc402f7c40afc7a50979aaca79697c1  gpu_full_audit
d95b35b28d01da5b753d054435643bff0530a4ab95393390b7b544f833f44d98  audit-results.json
73005e7b0df19ba77a484b328411eb4104c03f101cbcabed9edbdddaf589d518  tau_smc_r013_original
```
