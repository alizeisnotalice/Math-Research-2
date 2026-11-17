# r025 rank-biased fixed-prefix audit receipt

Date: 2026-11-17, local machine.

Executed:

```sh
python3 -B verify_rankbiased_prefix.py
```

The run used exact `fractions.Fraction` arithmetic and completed with exit
status 0. It enumerated all valid ordered prefixes at `n=2, r=1,...,4` and
`n=3, r=1,...,5`. For each state space it checked row sums, target support,
and every nonzero detailed-balance equality, first with the r025 rank weights
and then with an independently chosen fixed nonuniform rank-weight vector.
The largest enumerated case was `(n,r)=(3,5)`, with 3808 states and 62936
nonzero transition entries. The full result is `exact-prefix-balance.json`
with SHA-256
`51eaf1b22c7470cf8c9b43991309266d762c4b44a8f6ac53581591ebf291b321`.

The Swift host program was compiled but not executed; no Metal kernel or GPU
task was started. Compiler command:

```sh
TMPDIR=/Volumes/data/CodexTemp/T2 XDG_CACHE_HOME=/Volumes/data/CodexTemp/cache \
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-r025 \
swiftc -O -o tau2_rankbiased_prefix_r025 tau2_rankbiased_prefix_r025.swift
```

It exited with status 0. Source SHA-256 is
`60746225d734d03a5567c0b2e45581e89cd8a39599605f6e1a4397b612e23334`; binary
SHA-256 is `94cd9a14484b29d08662bbf7a22193b789c20503e9411cdd0cea7c4b09bb779d`.
The binary was not retained as a research artifact for a future commit.

An independent source-only Metal compile was attempted with `xcrun metal`,
without running a GPU kernel. This machine's active developer tools do not
provide that compiler (`xcrun: unable to find utility "metal"`, exit 72).
The successful `swiftc` build checks the Swift host program, but Metal source
inside the runtime string is compiled only when a future Metal library is
created; it is therefore not a substitute for a Metal compile or a GPU trace.

This is a kernel-invariance and build audit only. It supplies neither a
mixing result nor an SMC accuracy or asymptotic-probability result.

## Runtime Metal trace

After the CPU audit, an isolated small GPU trace was run while the separate
production task remained in its own process:

```sh
/usr/bin/time -l ./tau2_rankbiased_prefix_r025 \
  3 17 1 2 2 n3-p17-r1-l4-trace.jsonl 9e25090000000001 strat 4 trace
```

This invokes `MTLDevice.makeLibrary`, so it actually compiled the embedded
Metal source at runtime. It exited 0 in 3.77 seconds; the resource account
records 40,239,824 bytes peak footprint, zero swaps and zero signals. Raw,
trace and resource SHA-256 values are respectively
`b9034d0795cecd32f4799f26404168f5b1a616a536354bb736e26ddeefdeb11a`,
`dc6b27ad38ace87bd3660532896e102d0a2e8edc8d7cdd7f2fc532a5a07d21ff`, and
`199a3f016d80b257bc2b05ccb9ec7f84144f0104064a63546619407d674f4dd4`.

`verify_rankbiased_trace.py` independently replays every resampling, append,
rank-biased proposal, acceptance/rejection and occupancy/boundary state with
integer SplitMix64 arithmetic. It exactly matches all 9 ranks and all 17
particles. The trace contains 146 root--nonroot active proposals (48
accepted and 98 rejected), 217 root--outside proposals (134 accepted and 83
rejected), and all other active/rejection categories too. The verifier output
SHA-256 is `378a058994c2719c46028de755928ff0ef82dbdf4a59129b3305f3770c3bc7ca`;
the verifier source SHA-256 is
`a50956046a5b20beb5a4c76f2b40dc1d5778f13ca76bb4448fa0eaaaecfda238`.

The trace has nine rank records because \(N=9\). Moves occur only after
resampling and appending at ranks 1 through 8: this is exactly
`17 * 8 * 4 = 544` proposal attempts. Rank 9 is terminal-only and stores an
\(N-1\)-vertex prefix. The final verifier appends its unique omitted vertex
only for audit purposes, verifies births `(1,2)` for every full word, and
reconstructs all `L`, `L*d^2`, `L*F0`, `L*G0`, L-distance histogram, `L^2`,
and root fields against both trace and raw JSON. It passed; the updated
verifier/output SHA-256 values are
`a35300d615c139eb0c4b07f85d8dde8fb8ff3bba2a79941679c98bb910bb010c` and
`6d4b0abbdd92acaf65da2ce50db07b6556ca883193d79f5fe1649a29cf671f1c`.
