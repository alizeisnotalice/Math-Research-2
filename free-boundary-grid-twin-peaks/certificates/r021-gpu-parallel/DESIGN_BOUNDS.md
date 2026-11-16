# r021 rank-parallel terminal inverse fibres

`baseline.swift` is a terminal-only optimization of the frozen r020 source
whose SHA-256 is `1f593425af5fe53f98cc3b01109d6498519ad2d555c2ba14332d715c07f26b48`.
The count, CDF, parent selection, prefix update, random-number and move
kernels are unchanged.

At terminal rank the new `inverseRanks` kernel writes the inverse final-word
map `rank[particle,vertex]`.  `inverseParticle` then handles one particle per
thread.  A DP predecessor test examines the at-most-four grid neighbours of a
vertex and their final ranks, instead of scanning every earlier tau vertex.
`inverseReduce` reproduces the original 256-particle group totals.  The
terminal-only rank table costs `P*256*2` bytes; the two `UInt64` particle
buffers cost `P*8` bytes each (0.5 MiB each at P=65536).

For `N<=256`, rank entries fit `UInt16`; the DP is restricted to d=2,3 and
has at most four masks.  A word has at most four endpoint candidates, hence
`m2<=4(N-1)<=1020` and `m3<=4(N-1)(N-2)<=259080`.  At terminal rank each
SMC selection count is in `{0,1}`.  Thus an inverse group is at most
`256*259080=66324480`; the existing `P^2N<=2^64-1` entry guard also makes the
host total safe (the largest bound over allowed square grids is below
`6.96e13`).

The public entry point also requires `P<=2^24`, before `run()` allocates
buffers.  This makes every `particle*256+vertex` address fit a Metal `uint`:
the largest valid address is `(2^24-1)*256+255=2^32-1`.  It also bounds the
old prefix buffers, terminal rank table and all group indices.  Actual
Metal-buffer allocation can of course impose a much smaller practical P on a
particular unified-memory device; the arithmetic bound is not a promise that
an 8 GiB terminal-rank buffer will allocate.

## Actual checks

The optimized source compiled with Swift 6.4.2 on Apple M6 / Metal 4.
For n=3, P=17, L=0/4/16; n=8, P=257, L=4; and n=16, P=1, L=0, all old
mathematical output fields and each m2/m3 GPU group were identical to r020.
An independent Python terminal DP also reproduced every group.  The n=16
case exercises the final omitted-prefix vertex.

For an n=12, P=4096, L=16 terminal run with identical seed, all mathematical
fields were again token-identical.  Under the same concurrent-GPU load, old
and new recorded GPU times were 5.9925626667 and 0.4483013333 seconds; wall
times were 6.3014999628 and 0.8572199345 seconds.  This is a local benchmark,
not a probability calculation or a standalone hardware comparison.

All test processes exited.  The r020 independent host checker is at
`verify_inverse.py`.
