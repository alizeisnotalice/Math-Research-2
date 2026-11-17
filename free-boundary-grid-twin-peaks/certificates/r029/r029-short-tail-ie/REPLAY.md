# Portable replay

Set `A022_R028_CERT_DIR` to r028's certificate directory and run from this
directory:

```sh
python3 -B short_tail_ie_check.py \
  --tail-dir "$A022_R028_CERT_DIR/r028-punctured-transport" \
  --output verify-short-tail-ie.json
cmp short-tail-ie-crosscheck.json verify-short-tail-ie.json
```

The exact formula controls only `h_1,h_2,h_3`.  In particular, a zero-prefix
term must not be read as asserting separate divisibility of
`C_H(I)/(M)_(t+1)`; the proof records the n=3 counterexample.
