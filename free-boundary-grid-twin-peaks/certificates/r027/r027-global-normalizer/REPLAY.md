# Portable replay

Run from this directory. The first three commands recompute the exact tables;
the n=3 invocation includes the independent 9-factorial enumeration. The
last command exhaustively checks the rank-one deletion and promotion fibres
for n=2 and n=3.

```sh
python3 -B punctured_normalizer_table.py 3 --output verify-n3.json
python3 -B punctured_normalizer_table.py 4 --output verify-n4.json
python3 -B punctured_normalizer_table.py 5 --output verify-n5.json
cmp punctured-table-n3.json verify-n3.json
cmp punctured-table-n4.json verify-n4.json
cmp punctured-table-n5.json verify-n5.json
python3 -B verify_rank_one_deletion_n3.py > verify-rank-one.json
cmp rank_one_deletion_n2_n3.json verify-rank-one.json
```

The deletion and stable-promotion formulae are exact identities on finite
connected simple graphs. The table only evaluates n=3,4,5. It localizes the
remaining normalizer problem to a punctured-grid root-transport inequality;
it neither proves that inequality nor resolves the twin-peaks problem.
