# Portable CRT replay

Set `A022_R028_CERT_DIR` and `A022_R029_CERT_DIR` to the corresponding
certificate roots and use a writable external replay directory.

```sh
python3 -B audit_short_tail_crt.py --raw raw \
  --pratt "$A022_R028_CERT_DIR/r028-punctured-layer/n7-center-crt/pratt-certificates.json" \
  --reference "$A022_R029_CERT_DIR/r029-punctured-orbits-n7/verified/v24-verified.json" \
  --r29-source "$A022_R029_CERT_DIR/r029-short-tail-mod/short_tail_mod.cpp" \
  --output "$A022_REPLAY_DIR/short-tail-crt-replay.json"
```

The public output differs from the historical audit JSON only in declared
path-string metadata (`certificate_path`, reference `path`, and three filename
fields per raw record).  All hashes, Pratt nodes, CRT values, bounds, and
integer constraints are compared exactly in the recorded replay.
