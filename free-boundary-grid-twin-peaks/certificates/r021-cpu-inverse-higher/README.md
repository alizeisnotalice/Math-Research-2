# r021 higher inverse fibres

This archive evaluates exact terminal inverse-fibre observables `m2` through `m5` for the least-numbered canonical shortest path of a uniform one-peak word. The target is `A2_distance_d / A1 = E_1[m_d]`. Its n=7 and n=11 entries are finite SMC diagnostics; they do not estimate the original K=2 law or prove mixing or an asymptotic statement.

## Source provenance

`running-source.cpp` is the source used for the n=11 production batch. Its SHA-256 is `38388ec36db5a11431d8c574e7f1754c70a5d9600714b4896f46f7dbdc449e54`, matching `pre-n11-source-binary.sha256`. The run binary had SHA-256 `3459c8d8296d120a88c4295ead8d413c2e231eb71b1557cc1175d389ce344dcf`.

`comment-fixed-source.cpp` has SHA-256 `c2f3d58c7348bcc102def2bbeff9176c5c2fd65c72a9c516128686a784af4f10`. It differs at exactly one comment line and is not the literal run source. See `SOURCE_AND_TABLE_RECEIPT.md`.

## Reproduction

The independent r020 3x3 exact oracle is vendored as `reference-r020-exact_inverse_fibre_n3.py`.

```sh
python3 make_n3_higher_table.py
python3 verify_n7_higher.py
python3 summarize_higher.py n11-p65536-r256-l16.csv --output recomputed-n11.json
cmp recomputed-n11.json n11-strict-summary.json
```

The n=7 verifier uses `reference-r018-n7-distance-histogram.json`, checking all 32 rows, metadata and seed formula, zero dead counts, proposal counts, and stored ratio identities. The n=11 Decimal summary checks all 256 rows and is a diagnostic, not a confidence interval.

The n=11 record has `P=65536`, `R=256`, `L=16`, `T=12`; its resource report gives 878.36 wall seconds, 210,796,544-byte maximum RSS, and zero swaps. The computation is finished; no production binary is archived here.
