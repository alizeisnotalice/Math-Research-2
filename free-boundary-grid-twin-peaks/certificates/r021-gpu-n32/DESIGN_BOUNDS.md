# r021 n<=32 ordered-prefix extension

This directory is an isolated extension of the frozen n<=16 r021 source
`6aa1f829dbe97c141781c7440c0f99ef62fa966cbe2c665fab7432a48b2ec9f8`.
It does not modify the production source or directory.

## Representation and bounds

The existing `W=32` occupancy and external-boundary words already provide
1024 bits.  This extension changes only the prefix and terminal-rank stride
from 256 to 1024, admits `n<=32`, and changes the pre-allocation particle cap
to `P<=2^22`.  Hence every Metal `uint` index
`particle*1024+vertex` is at most `2^32-1`.

`UInt16` holds both vertices and terminal ranks because `N<=1024`.  With
coordinates centred as `2x-(n-1)`, a root has squared norm at most 1922.
For a 256-particle terminal group,

| quantity | maximum at n=32 |
| --- | ---: |
| root F group sum | 984064 |
| absolute root G group sum | 492032 |
| distance-squared group sum | 492032 |
| inverse m2 per word | 4092 |
| inverse m3 per word | 4182024 |
| inverse m3 group sum | 1070598144 |

At `P=2^22`, full host totals are at most 16122904576 for F and
17540679991296 for m3, all below `UInt64.max`.  The arithmetic address cap
does not assert that the hardware can allocate the corresponding buffers:
each of the two prefix buffers and the rank buffer alone would be 8 GiB.

## Independent trace audit

`verify_n32_trace.py` reconstructs every word from the exact graph, every
allowed set, stratified parent, random draw, move, all 32 occupancy/boundary
words, and terminal root quantities.  It compares real high words rather than
padding them.  `verify_n32_inverse.py` independently recomputes terminal
d=2,3 inverse fibres; on K=2 traces it instead verifies the required inactive
null fields and zero group buffers.

Actual full traces were executed for n=24 and n=32, each at K=2 and one peak,
with P=1 and positive move counts.  n<=16 same-seed n=3 and n=16 regression
traces retained exact old mathematical fields and inverse group totals.  These
are implementation checks, not estimates of either conditional law.

The stronger P=17 audit is in [MULTIPARTICLE_AUDIT.md](MULTIPARTICLE_AUDIT.md).
It includes complete n=24/n=32 K=2 and one-peak traces, non-dead Double
underflow checks, and three destructive-copy rejection tests.
