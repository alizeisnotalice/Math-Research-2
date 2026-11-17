# Independent replay

The audit is pure Python integer arithmetic over archived files. It does **not**
compile or run an enumerator, and it does not claim to verify the historical
binary, which is intentionally absent from this public package.

After placing this directory in the repository as
`free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public`, run from the
repository root:

```sh
python3 free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/audit_orbits.py \
  --raw free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/raw \
  --r29-verified free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/r29-verified \
  --center free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/center/short-tail-n7-center-crt-audit-public.json \
  --pratt free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/pratt/pratt-certificates.json \
  --historical-audit free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/historical/n7-short-tail-orbits-crt-audit.json \
  --historical-center free-boundary-grid-twin-peaks/certificates/r031/n7-orbit-public/historical/short-tail-n7-center-crt-audit-original.json \
  --output /Volumes/data/CodexTemp/T2/n7-orbit-public-replay/n7-short-tail-orbits-crt-audit-public.json
```

The command must print `"status": "passed"`. It checks all 36 receipt
triples, source and declared-source hashes, nine r029 layer histograms, the
49-node Pratt proof, direct and Garner CRT reconstruction, all ten orbit
aggregates, and normalized agreement with the historical r031 output.

Before replay, verify the packaged whitelist and hashes in its directory:

```sh
shasum -a 256 -c SHA256SUMS
```
