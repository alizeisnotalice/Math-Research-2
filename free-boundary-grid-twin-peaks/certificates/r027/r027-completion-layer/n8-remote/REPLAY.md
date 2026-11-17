# Portable n=8 receipt and CRT replay

Run from `r027-completion-layer`.

```sh
python3 -B n8-remote/audit_receipts_portable.py > n8-remote/portable-audit.json
python3 -B n8-independent/verify_n8.py \
  pratt-certificates.json \
  n8-remote/remote-results/n8-p0.json \
  n8-remote/remote-results/n8-p1.json \
  n8-remote/remote-results/n8-p2.json \
  n8-remote/remote-results/n8-p3.json > n8-independent/verify-replay.json
cmp n8-independent/n8-crt-verified.json n8-independent/verify-replay.json
```

The receipt audit verifies the archived source hash, bridge and n=8 fields,
command/exit/resource hashes, zero swaps, and the completed parent process.
It intentionally omits the nonportable executable, so it does not verify the
remote binary hash. The independent verifier checks the Pratt tree, a fresh
6x8 row-transfer histogram, two CRT formulae, and the proved `H <= 56!`
uniqueness bound. This is an exact finite top-row completion count only, not
a full `A2(8)` count or an asymptotic conclusion.
