# r023 exact integer and address bounds

For `n<=32`, `N=n²<=1024`, `d<=2n-2<=62`, and the existing tau=2 lift
invariant gives `1<=L<=N-3<=1021`.  Thus an individual histogram bin is at
most `256*1021=261376`, while a group L² sum is at most
`256*1021²=266864896`.  Across `P<=2²²` particles the L² total is below
`P*1021²<2^64`; UInt64 is exact throughout.

There are `G=ceil(P/256)<=16384` groups and `2n-1<=63` bins, so the new
histogram has at most `1,032,192` UInt64 cells (8,257,536 bytes).  Its address
`group*(2n-1)+d` is below 2³².  Terminal-prefix access is
`particle*1024+rank`; the inherited pre-allocation `P<=2²²` bound keeps this
Metal uint index below 2³².  The kernel reads only ranks `0,...,N-2`, never
fabricates an omitted terminal vertex.

The host checks each group `sum_d hist[d]=sum L`, `hist[0]=hist[1]=0`, and the
displayed L² group bound.  Therefore the JSON raw values support exact integer
tail numerators and exact first/second birth-rank moments before any floating
normalization.
