# Replay

The audit needs Python 3 standard library only and reads shipped raw JSON; it
does not compile or run the production enumerator.

```sh
OUT=/path/to/precreated/r035-audit-output
test -d "$OUT"
python3 scripts/audit_spatial_crt.py --root . --output "$OUT/n7-spatial-crt-audit.json"
```

The optional output is deliberately outside this frozen certificate.
