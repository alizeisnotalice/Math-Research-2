# Large-side-only row-orbit reduction

Fix the original labelled bipartition `A+B` with `|A|=a` and `|B|=b`.
Its edge matrix has `a` rows and `b` columns. For the original small-side
families `a=4,5`, enumerate instead the `b` rows of the transposed matrix.
Equivalently, sort the neighbourhood masks of vertices in **B only**. This is
the action of `S_B`, not an action of `S_A x S_B`.

Vertex relabelling preserves every global prefix count `h_j`; hence it
preserves `Rj7`. Each `S_B` orbit has exactly one nondecreasing tuple of `b`
masks in `[0,2^a)`. If its mask multiplicities are `c_1,...,c_t`, its coverage
weight is `b!/(c_1! ... c_t!)`. Thus every original labelled edge matrix
occurs exactly once after weights, including empty and disconnected graphs,
and the total weight is `(2^a)^b=2^(ab)`.

For original `a=4`, sort the transposed `b=7` rows: there are
`C(16+7-1,7)=170544` representatives, weighted total `2^28`.
For original `a=5`, sort the transposed `b=6` rows: there are
`C(32+6-1,6)=2324784` representatives, weighted total `2^30`.
Their total is `2495328`, compared with `11716640+10424128=22140768` for
sorting the original small-side rows: a factor `8.873...` reduction relative
to that prior row-orbit plan.

No claim is made that both sides have been quotiented simultaneously. The
rank intervals are ordinary one-side multiset ranks, so `[start,stop)` shards
remain disjoint and checkpointable.
