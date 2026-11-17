# Portable replay

Run from this directory.  No executable, remote address, or local source tree
is needed for this verification.

```sh
python3 -B verify_n7_crt.py pratt-certificates.json \
  residues/n7-v24-p0.json residues/n7-v24-p1.json \
  residues/n7-v24-p2.json residues/n7-v24-p3.json \
  > verify-n7-center.json
cmp n7-center-crt-verified.json verify-n7-center.json
```

The final `cmp` is exact JSON equality, including the input hashes and the
independently regenerated histogram.  It does not rerun the memory-intensive
modular production calculation.
