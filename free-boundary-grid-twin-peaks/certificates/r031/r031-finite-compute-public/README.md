# r031 finite-compute public reproduction subpackages

This directory is a copy-ready, deliberately narrow public payload for two
finite computations. It contains no `n=7` production result, no production
binary, no trace, and no telemetry flood.

| Package | Replay result | Mathematical scope |
| --- | --- | --- |
| `short-gap-h45/` | 124 four-prime residues / 31 points match complete r029/r030 L references | finite `n=3,4,5` implementation check only |
| `merge-rb-n32/` | raw hash, frozen f1593 audit, and recorded-field comparison pass | finite n=32 input/invariant audit only |

Start with each package's `REPLAY.md`. `SHA256SUMS` lists every shipped file
and `COPY_WHITELIST.txt` is the exact r031 certificate-copy list. The replay
commands require only Python 3; A additionally documents a C++20 build.

No numerical item here establishes an asymptotic result, a probability claim,
a cone condition, SMC accuracy, mixing, or a solution of the underlying
problem.
