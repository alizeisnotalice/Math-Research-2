# Portable replay

This package contains the frozen final modular source, four residue
checkpoints for each n from 3 through 7, and the small certificates needed
to check the CRT width and the four moduli. It intentionally contains no
compiled executable, cache, old implementation, standard-output trace, or
n=8 result. The n=8 plan is a memory-guard calculation only.

Run these commands from this directory.

```sh
python3 -B pratt_verify.py > verify-pratt.json
python3 -B prime_audit.py > verify-primes.json
python3 -B crt_verify.py modular-runs > verify-crt.json
cmp pratt-verify.json verify-pratt.json
cmp prime-audit.json verify-primes.json
cmp crt-verification.json verify-crt.json
```

For the remote n=8 finite completion result, follow
`n8-remote/REPLAY.md`. Its model payload (5,803,578,448 bytes) is distinct
from the observed whole-process RSS and peak-memory-footprint values recorded
in the receipt.

For an independent small finite execution, compile into a temporary
location and write a fresh checkpoint there. The program refuses to replace
an existing checkpoint.

```sh
clang++ -std=c++20 -O3 -pthread completion_layer_mod.cpp -o /path/to/completion_layer_mod
/path/to/completion_layer_mod 7 0 6 /path/to/n7-p0.json
```

The final source hash before archival copying is recorded in
`UPSTREAM_SOURCE_SHA256.txt`. The 20 recorded passes and CRT check establish
only the stated finite n=3..7 equalities. They do not execute n=8 and make
no asymptotic claim.
