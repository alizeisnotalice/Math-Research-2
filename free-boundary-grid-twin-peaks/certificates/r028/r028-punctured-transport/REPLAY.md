# Portable finite tail-DP replay

Run from this directory. `A022_R027_CERT_DIR` is the r027 certificate root.

```sh
python3 -B punctured_tail_dp.py 3 --output verify-tail-n3.json
python3 -B punctured_tail_dp.py 4 --output verify-tail-n4.json
cmp tail-n3.json verify-tail-n3.json
cmp tail-n4.json verify-tail-n4.json
python3 -B audit_punctured_tail_dp.py \
  --r027 "$A022_R027_CERT_DIR/r027-global-normalizer" > verify-tail-audit.json
cmp tail-dp-audit.json verify-tail-audit.json
```

The prefix recurrence and fibre identities are exact. The global condition
`T[t+1] <= 4*T[t]/5` is an explicitly unproved sufficient condition, not an
equivalent formulation or a result of the finite n=3,4 tables.
