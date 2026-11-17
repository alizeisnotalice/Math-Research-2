# Portable replay

Run from this directory.  Set `A022_R028_CERT_DIR` to r028's certificate
directory and `A022_R027_CERT_DIR` to r027's certificate directory.

```sh
mkdir -p "$A022_REPLAY_DIR"
clang++ -std=c++17 -O3 -DNDEBUG tail_full_dp.cpp -o "$A022_REPLAY_DIR/tail_full_dp"
"$A022_REPLAY_DIR/tail_full_dp" 3 all > "$A022_REPLAY_DIR/n3.json"
"$A022_REPLAY_DIR/tail_full_dp" 4 all > "$A022_REPLAY_DIR/n4.json"
for v in 0 1 2 6 7 12; do "$A022_REPLAY_DIR/tail_full_dp" 5 "$v" > "$A022_REPLAY_DIR/v$v.json"; done
cmp tail-full-n3.json "$A022_REPLAY_DIR/n3.json"
cmp tail-full-n4.json "$A022_REPLAY_DIR/n4.json"
python3 -B audit_tail_full_dp.py "$A022_REPLAY_DIR/n3.json" --reference-n 3 --reference-dir "$A022_R028_CERT_DIR/r028-punctured-transport"
python3 -B audit_tail_full_dp.py "$A022_REPLAY_DIR/n4.json" --reference-n 4 --reference-dir "$A022_R028_CERT_DIR/r028-punctured-transport"
for v in 0 1 2 6 7; do cp "$A022_REPLAY_DIR/v$v.json" "$A022_REPLAY_DIR/tail-full-n5-v$v.json"; done
cp "$A022_REPLAY_DIR/v12.json" "$A022_REPLAY_DIR/tail-full-n5-center.json"
python3 -B summarize_n5_d4.py --table "$A022_R027_CERT_DIR/r027-global-normalizer/punctured-table-n5.json" --input-dir "$A022_REPLAY_DIR" --output "$A022_REPLAY_DIR/summary.json" > "$A022_REPLAY_DIR/summary.stdout"
cmp n5-tail-summary.json "$A022_REPLAY_DIR/summary.json"
```

The conclusion is an exact finite n=5 fact, not an all-n contraction theorem.
