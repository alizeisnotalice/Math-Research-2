# r022 index and integer bounds

Let `P<=2^22`, `G=ceil(P/256)`, `n<=32`, and `m=2n-1<=63`.

The source prefix is addressed at `particle*1024+k`.  With
`particle<=2^22-1` and `k<=m-1<=62`, this is strictly below `2^32`; the
existing global P guard is evaluated before `run()` allocates any buffer.
The same bound already covers the terminal rank table.

The new output address is `group*m+(d-1)`.  Since `G<=16384`, it is less than
`16384*63=1,032,192`, far below a Metal `uint`.  The buffer has at most
`1,032,192` UInt64 entries, i.e. 8,257,536 bytes.  Each stored entry is a sum
of at most 256 boolean values, so it lies in `[0,256]`; an across-group raw
sum lies in `[0,P]`, therefore below `2^22` and far below UInt64 overflow.

For `n>=2`, `2n-1<=n^2-1`.  Consequently all requested positions are present
in the actual terminal prefix (which stores its first `N-1` values); no
inferred final vertex or fabricated high word is used.

The equality edge `n=2`, where `2n-1=N-1=3`, was also run at P=17 and fully
replayed by the independent prefix/group verifier; see
`results/n2-p17-geodesic-receipt.json`.

The actual invalid input `n=2,P=4,194,305` was rejected with exit code 133
before output-file creation; see `results/invalid-particle-receipt.json`.
