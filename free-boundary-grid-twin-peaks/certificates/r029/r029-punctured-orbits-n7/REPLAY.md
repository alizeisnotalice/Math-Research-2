# Portable selected-hole replay

Set `A022_R028_CERT_DIR` to r028's certificate directory.  The following
replays each of the nine new orbit representatives; replace the four residue
paths by the corresponding `residues/v*/p0.json`--`p3.json` directories.

```sh
python3 -B verify_selected_hole.py --vertex 7 \
  "$A022_R028_CERT_DIR/r028-punctured-layer/n7-center-crt/pratt-certificates.json" \
  residues/v7/p0.json residues/v7/p1.json residues/v7/p2.json residues/v7/p3.json \
  > verify-v7.json
cmp verified/v7-verified.json verify-v7.json
```

The same verifier is run for `v=0,8,14,15,16,21,22,23`; `v=24` uses r028's
central residues.  `verified/n7-orbits-verified.json` records the independently
checked D4 partition and totals.  It is a frozen result record; its original
orchestrator's historical-binary checks are intentionally not represented as
portable binary replay.
