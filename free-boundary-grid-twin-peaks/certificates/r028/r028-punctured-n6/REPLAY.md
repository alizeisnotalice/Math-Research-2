# Portable finite replay

Build outside this certificate directory and run the asserted and optimized
variants. `A022_R027_CERT_DIR` must name the final r027 certificate directory
that contains `r027-global-normalizer`.

```sh
clang++ -std=c++17 -O3 -DNDEBUG punctured_normalizer_dp.cpp -o /path/to/punctured_dp
python3 -B verify_r27_n3_n5.py /path/to/punctured_dp "$A022_R027_CERT_DIR/r027-global-normalizer"
/path/to/punctured_dp 6 > verify-n6.json
cmp punctured-table-n6.json verify-n6.json
clang++ -std=c++17 -O3 punctured_normalizer_dp.cpp -o /path/to/punctured_dp_assert
/path/to/punctured_dp_assert 6 > verify-n6-assert.json
cmp punctured-table-n6-assert.json verify-n6-assert.json
cmp verify-n6.json verify-n6-assert.json
python3 -B oracle/subset_birth_oracle.py \
  --reference-dir "$A022_R027_CERT_DIR/r027-global-normalizer" \
  --output verify-subset-oracle.json
python3 -B oracle/verify_oracle_output.py \
  oracle/subset-birth-oracle.json verify-subset-oracle.json
```

This is an exact finite n=6 table. It does not establish a uniform punctured
normalizer comparison or an asymptotic statement.
