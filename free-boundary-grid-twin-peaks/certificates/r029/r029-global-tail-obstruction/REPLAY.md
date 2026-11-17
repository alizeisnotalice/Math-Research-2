# Portable replay

Run from this directory.  Set `A022_R028_CERT_DIR` to r028's certificate
directory; only its small `r028-punctured-transport/tail-n3.json` input is
needed for the finite terminal-bucket check.

```sh
python3 -B check_obstruction_maps_n3.py > verify-n3-map.json
cmp n3-map-check.json verify-n3-map.json
python3 -B verify_terminal_tail_n3.py \
  "$A022_R028_CERT_DIR/r028-punctured-transport/tail-n3.json" \
  verify-terminal-tail-n3.json
cmp terminal-tail-n3.json verify-terminal-tail-n3.json
python3 -B verify_kappa_4_5_threshold.py > verify-kappa.json
cmp kappa-4_5-threshold.json verify-kappa.json
```

The proof in `global-tail-obstruction.md` is the general result.  The two
finite scripts only audit its n=3 maps and explicit n=100 endpoint arithmetic.
