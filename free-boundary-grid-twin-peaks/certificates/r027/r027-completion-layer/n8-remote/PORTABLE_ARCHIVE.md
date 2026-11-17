# Remote n=8 modular layer passes

This receipt contains four independently checkpointed modular residue passes
for the selected top-row completion count at `n=8`, after a four-prime n=7
bridge and boundary/plan control tests.  It is **not** an integer count
certificate by itself: the four residues require independent CRT and bound
verification.

Frozen source SHA-256:
`221ae28c1e074125b2a2fbc9794668c559d3b95f68e6f66de5740304746a6684`.
The output JSON files are `remote-results/n8-p0.json` through `n8-p3.json`.
Each has its own start/exit and `/usr/bin/time -l` resource receipt.

Run `python3 audit_remote.py` to check source/output/resource hashes,
per-prime fields and field ranges, n=7 bridge residues, control tests,
checkpoint exits and zero swaps.  This audit intentionally does not perform
CRT reconstruction; use the independent verifier supplied by the main line.

The model payload peak is 5,803,578,448 bytes (about 5.405 GiB), bounded by
the program's 7 GiB allocation-model guard.  The resource receipt's RSS and
peak-footprint values are whole-process observations and must be reported
separately; they are not covered by that modeled-payload guard.
