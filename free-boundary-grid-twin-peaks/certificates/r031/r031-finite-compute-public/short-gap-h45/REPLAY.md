# Public finite replay: h1--h5 modular rolling implementation

This package contains the frozen final C++ source and the 124 saved small
checkpoints for 31 puncture points: every vertex for `n=3,4`, and the six D4
representatives `0,1,2,6,7,12` for `n=5`.  Each point has all four prime
residues plus its saved stdout, resource receipt, and exit status.  It does
not contain, start, or authorize an `n=7` run.

## Requirements and portable command

Use Python 3 and a C++20 compiler.  Supply checked-out / archived r029 and
r030 certificate directories explicitly; no source file in this package names
a private absolute path.

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread \
  source/short_gap_h45_mod.cpp -o replay/short_gap_h45_mod
python3 source/verify_h1_h5_portable.py \
  --source source/short_gap_h45_mod.cpp \
  --checkpoints raw/checkpoints \
  --r029-tail /path/to/r029/r029-tail-full-n5 \
  --r030-geometry /path/to/r030/r030-gap-geometry \
  --output replay/finite-crosscheck.json
```

The verifier CRT-reconstructs every `(h1,...,h5)` from the four saved residues,
compares it to the complete r029 L histogram, and also reads r030's complete
n=3/n=5 geometry tables and requires their matching L/A2/B fields.  It does
not recompute any checkpoint; it verifies saved finite data.

## Scope

Successful replay establishes only that these finite saved checkpoint fields
agree with the specified independent finite references.  It is not a proof
for `n=7`, a limiting result, or a general theorem.

## Public whitelist

* `source/short_gap_h45_mod.cpp` — final source, SHA-256
  `08c31508f6d3cb0ae50f74973aa372a32e3bac8439d9d6bfd51d4bf6c0312fc8`.
* `source/verify_h1_h5_portable.py` — path-parameterized public verifier.
* `raw/checkpoints/` — 124 JSON residues and their 372 stdout/resource/exit
  receipts.
* `reports/` — original finite report copies and newly produced replay result.
* `PROVENANCE.json`, `SHA256SUMS`, and this file.

The historical binary is intentionally omitted.  Its observed SHA-256 is
recorded in `PROVENANCE.json`, not presented as a reproducible public build.
