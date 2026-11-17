# Portable exact replay

Use an existing Python3 interpreter. No third-party packages or randomness
are needed. From the repository root, choose a fresh external output directory:

```sh
python3 -B free-boundary-grid-twin-peaks/certificates/r033/replay.py \
  --out /Volumes/data/CodexTemp/T2/a022-r033-new-replay
```

The runner copies scripts that write beside themselves into that output
directory. It runs six audits and checks exact JSON. For the independent path
audit only, wall time/version/platform fields in `run` are excluded from
comparison; the mathematical payload and source hash must still match.

It replays full n3 A2 encoding, an independent B-only encoding, the Fraction
coupling check, n3 permutation root-pair counts, n4 subset DP, and n7 CRT from
saved modular outputs. It does not rerun the large n7 production calculation.
The original C++ production source and its smaller smoke are documented in
n7-global-h45/REPLAY.md. All child exits are waited for and recorded.

Each result has limited scope: the general path injection and coupling are
proved analytically; small computations verify their implementation and edge
cases. The n7 certificate is an exact finite integer reconstruction; it is
neither an asymptotic estimate nor an all-layer tail proof.
