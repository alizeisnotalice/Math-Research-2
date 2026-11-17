# Exact necessary-cone replay

Run from `r027-l64-recovery`:

```sh
python3 -B cone/audit_cone.py canonical/n24-move64-canonical-eps-1-2.json > cone/verify.json
cmp cone/l64-cone.json cone/verify.json
```

The script uses `Fraction` throughout. Its positive slack means that this
stored plug-in pair does not violate that necessary cone; it is not a
certificate of sampling accuracy, mixing, or an asymptotic statement.
