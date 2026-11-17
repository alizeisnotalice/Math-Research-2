# Calibration receipt

Frozen r026 trace source: SHA256
934275c0da831ea81c22ac78fcd4175517723afef9343553796f77e71d3bce63.
The observation-only calibration derivative has SHA256
3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6;
its binary has SHA256
5de42d7ec2eb946cfb805b59df658cf2f25d87e6018d4561d3e7a7b3a2a9d598.
The derivative adds host serialization of GPU terminal payload m statistics;
it does not change prefix sampling, RNG, resampling, MH, or terminal payload
construction.

Completed bounded runs: n=3 P=4096 R=64 L=16 domain 9e260b0000000001,
and n=7 P=16384 R=64 L=16 domain 9e260b0000000002. Both exited 0. A
current derivative n=3 P=17 trace was independently audited by
audit_derived_trace.py, including all four new m fields from its GPU payload.

aggregate_calibration.py hard-fails absent/mismatched n=3 or n=7 r023
top2-D4 baseline, and records the actual source hash, normalizer estimates,
relative errors, A2/B comparison, ESS and max share. Audit success concerns
data/source/baseline consistency only, not statistical accuracy. All
calibration processes were checked absent after the runs.
