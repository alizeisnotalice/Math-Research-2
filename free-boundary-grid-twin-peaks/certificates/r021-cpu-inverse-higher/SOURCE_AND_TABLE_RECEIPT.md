# Source and 3x3-table receipt

The production receipt records source SHA-256 `38388ec36db5a11431d8c574e7f1754c70a5d9600714b4896f46f7dbdc449e54` and binary SHA-256 `3459c8d8296d120a88c4295ead8d413c2e231eb71b1557cc1175d389ce344dcf`. Applying the sole reverse comment edit to the later source yields `running-source.cpp` with exactly the production source hash.

```diff
-// Exact fibre size of the canonical greedy shortest-path lift, for d=2 or 3.
+// Exact fibre size of the canonical greedy shortest-path lift, for d=2 through 5.
```

No algorithmic line, compiler flag, or run CSV changed with that edit.

The retained `n3-python-higher-fibres.csv` is the independent Python oracle: 29,568 one-peak words plus one header, SHA-256 `e71cdfdf599957a32ecbcadf375b0b9a49f640f01b48998397596f75755560be`. Its sums are `(m2,m3,m4,m5)=(63552,64512,15744,0)`.

At production time, the C++ table made by the running binary had that same SHA-256 and compared byte-for-byte with the Python oracle. After the comment-only correction, the recompiled C++ n=3 table again compared byte-for-byte with the oracle. The redundant C++ tables and binaries are omitted; this receipt preserves the two-table equality evidence while retaining one complete table.
