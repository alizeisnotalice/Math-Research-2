# Portable exact cross-check

Run from this directory. The two arguments point to the corresponding
certificate directories from rounds 27 and 28.

```sh
python3 -B punctured_h_exact_check.py \
  --r027 "$A022_R027_CERT_DIR/r027-global-normalizer" \
  --r028 "$A022_R028_CERT_DIR/r028-punctured-n6" > verify-crosscheck.json
cmp punctured-h-crosscheck.json verify-crosscheck.json
```

The proof note establishes exact finite recurrences and CRT-width conditions.
The computation independently checks all D4 representatives only through
n=6. The proposed lower bound mentioned in the note remains a conjecture;
neither this finite computation nor the CRT plan proves it.
