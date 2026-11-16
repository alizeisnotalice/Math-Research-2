# Exact `UInt64` range and index audit

Let `N=n^2<=256`, fix a one-peak terminal output `sigma`, and fix `d` in
`{2,3}`.  For any canonical endpoint candidate `v`, `u=sigma[0]` and the
tau-word are fixed.  A preimage differs only by the positions and order of
the `d-1` labelled deleted internal points among the `N-1` slots after `u`.
Consequently its count is at most `(N-1)_(d-1)`.  Birth tests only discard
interleavings.  There are at most four candidates because `v` must be a
neighbor of the last internal point.  Hence

| quantity | exact upper bound at N=256 |
| --- | ---: |
| one-candidate d=2 DP | 255 |
| `m2(sigma)` | `4*255 = 1020` |
| one-candidate d=3 DP | `255*254 = 64770` |
| `m3(sigma)` | `4*255*254 = 259080` |
| one terminal workgroup's m3 sum | `256*259080 = 66324480` |

At rank `r=N`, `a_i` counts allowed unoccupied vertices of a prefix with
`N-1` occupied distinct vertices.  There is at most one.  Thus the final
kernel multiplication is by 0 or 1, rather than by a potentially large SMC
weight.  A workgroup is launched by one GPU thread, which serially visits at
most 256 particles; it is not a 256-thread terminal reduction.

The host then reduces **all** `ceil(P/256)` group values, so its raw total is
bounded by `P*259080`, not by one group.  The archived entry-point precondition requires `P^2*N <= 2^64-1`.
Since n>=2 and hence N>=4, this already implies P<=2147483647, safely within UInt32.  For fixed `N`, take
`P_max=floor(sqrt((2^64-1)/N))`.  The resulting upper bound
`P_max*4*(N-1)*(N-2)` increases over the permitted square-grid range and is
largest at `N=256`:

```
P_max = 268435455
P_max * m3_max = 268435455 * 259080 = 69546257681400 < 2^64.
```

Thus every nonnegative addition in `inv2`, `inverseTerm`, and the final host
reduction is safely within `UInt64`; the `P<=UInt32.max` condition also makes
the `uint` Metal particle-count argument non-truncating.

Index audit of the actual kernel:

* `inverseTerm` invokes only `inv2(...,2)` and `inv2(...,3)`, so `plen<=2`,
  `full<=3`, and `cur[4]`, `nxt[4]` cover every mask.
* `r` ranges from 1 through `N-1`; for nonzero reachable states,
  `ti=r-popcount(mask)` is the number of placed tau points.
* `tv2` is called only when `ti<T`.  Its mapped original-word index is at most
  `N-1`.
* `sv2` reads `pref[id*256+k]` only when `k+1<N`; for N=256 that means
  `k<=254`.  The sole `k=255` case scans the 8 occupancy words to recover the
  omitted vertex, and no prefix-array access occurs.

This audits the fixed d=2,3 kernel only.  It is not a claim about a generic
larger-d extension.
