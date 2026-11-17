# Portable finite modular replay

This package archives 18 D4-representative cases for n=3..6, with four prime
residues per case. It excludes the executable and all duplicated stdout.
`A022_R027_CERT_DIR` and `A022_R028_CERT_DIR` must point to the final
certificate directories of the named rounds.

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread punctured_layer_mod.cpp -o /path/to/punctured_layer_mod
/path/to/punctured_layer_mod --boundary-test /path/to/boundary-test.json
/path/to/punctured_layer_mod 6 14 0 1 /path/to/n6-v14-p0.json
python3 -B verify_residues.py --runs runs \
  --table-root "$A022_R027_CERT_DIR/r027-global-normalizer" \
  --table-n6 "$A022_R028_CERT_DIR/r028-punctured-n6/punctured-table-n6-assert.json" \
  --output verify-crosscheck.json
cmp finite-crosscheck.json verify-crosscheck.json
python3 -B resource_model.py > verify-resource-model.json
cmp resource-model.json verify-resource-model.json
```

The n=7 center data here is a guarded resource model only. A production n=7
result requires four completed residues and an independent CRT verification;
none is asserted by this package.
