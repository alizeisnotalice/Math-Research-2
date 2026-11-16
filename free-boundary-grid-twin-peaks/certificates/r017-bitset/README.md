# 第17轮：位集 GPU 核独立审计

This is the curated archive of independent temporary audits. It does not assert any probability estimate. Full reproducible rank traces remain in CodexTemp; their inputs, sizes and original hashes are in trace-inputs.json.

## Frozen source identities

- Formal host/Metal source copied to `bitset_smc_frozen_final.swift`:
  `a63ce678ad85153fb29da3e3e74988615f4c54338f1bc0a4ad387925767a6172`.
- Its embedded Metal source (the bytes strictly between `let metal = """` and
  the following `"""`, excluding those delimiters) has SHA-256
  `6c59a560e34e746633fa5e6c5ba886c759786d9da2df4d313af9a32258a32181`.
- Rank-trace source used for the trace runs:
  `bf8da66017096d144899a6d92e12cd36a326f20800215cd4b6faf17513d17f95`.

The final host source differs from the immediately preceding copy only in the
algebraic spelling of `u`; the embedded Metal bytes are unchanged.  Thus the
trace audit is an audit of the same GPU kernels, while the terminal CLI tests
below use the final formal host source itself.

## Actual checks

`verify_rank_traces.py` is an independent Python implementation of the grid
neighbour boundary, candidate set, SplitMix/rejection RNG, stratified parent
selection, mutation, root updates, and death copying.  It checks every emitted
word of all 32 occupied and all 32 boundary words, roots, allowed vector,
group sums, CDF, parent and child state.  Its JSON certificate is
`independent-trace-audit.json` (SHA-256
`804ebcb6ac393b538421dc1037b45572f24e70e67d96962825169bd67ef1b1e6`):

- n=12,16,24,32 with tau=2;
- n=16,24,32 with tau=floor(N/2), N-1 and N;
- n=32 with tau=N+1.

All 14 traces passed.  In particular, the N-1 and N cases reach a real
zero-weight/dead state and the audit verifies the sentinel parent and exact
state copying through the remaining rank; this is distinct from floating-point
underflow.

The formal final CLI was independently compiled as
`bitset_smc_frozen_final` and run on n=32 at tau=512,1023,1024,1025.  Its
`terminal-log-moment-audit.json` (SHA-256
`5234297cb0df5da37f25ff06a219b85572b71c4f4cd9d742e8631a55ee293d21`)
records:

- tau=512: `dead=false`, `Z=0.0`, finite `logZ=-1394.267...`;
- tau=N-1,N: `dead=true`, `Z=0.0`, `logZ=null`;
- tau=N+1: `dead=false`, `Z=0.0`, finite logZ and
  `U`, `U_over_Z`, `scaled_U` all `null`.

The final n=12, tau=2, P=17, domain=1 CLI output matches the frozen old-kernel
terminal record exactly in `tau, replicate, seed_hex, Z, U, logZ, dead`.
For a live positive-range run it also directly checks `U=Z*(U_over_Z)`.

The terminal accumulator bound at n=32 is documented in the JSON certificate:
per particle at most `N*2(n-1)^2=1,968,128`, per 256-particle group at most
`503,840,768`, and under the host guard `P^2*N <= 2^64-1`, total raw bound is
`264,157,666,605,056`, all safely below `2^64`.  Cross-tau normalization must
therefore use log-sum-exp over finite `logZ` values, never a direct sum of
underflowed emitted `Z` values.

## Reproduction

From this directory run `python3 reproduce_rank_traces.py`. This compiles into
`/Volumes/data/CodexTemp/T2/a022-r017-rank-reproduction`, regenerates the 14
input cases, and invokes the independent verifier there. It writes no binary,
module cache or large trace into the repository. Use `--out` to select another
location under the data volume. The original trace audit was actually rerun
by the main Agent; this new convenience driver is an archival wrapper.

`python3 audit_terminal_log_moments.py` checks the committed short CLI records.
The n12-final-alltaus.jsonl and old12.jsonl files additionally compare
all 143 birth strata at P17 under the same seed, not just tau2.

The production host is `bitset_smc_checked.swift` (SHA256 65aefc323b20995564df7a65656b4bfffbc1570487f8be132fafa916225c2f67).
It adds `cb.status == .completed && cb.error == nil` after waiting; its Metal
bytes equal those of the frozen audited source. A live `Z=0` from underflow
retains finite logZ; an actual dead run has null logZ. Aggregate across birth
strata by summing stratum means, not by averaging over the number of strata.

Actual short audit processes exited. The million-particle production follow-up
is separate; no general mixing or asymptotic claim follows from kernel audits.
