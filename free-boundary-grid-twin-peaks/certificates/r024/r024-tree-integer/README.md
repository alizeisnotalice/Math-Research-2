# Independent integer verification of a bounded-degree auxiliary counterexample

This directory does not use the continuous-label integral messages in
`a022-tree-top2-scaling-r024`.  It counts rank permutations directly.

For a rooted tree, a state `(m,r,k,h)` records a subtree of size `m`, its
root rank `r`, the number `k` of already determined non-root peaks, and the
fact `h` that a child root has higher label.  When a child state of size `c`
is interleaved, choosing `s` child vertices above the old root has exactly

```
C(r-1+s,s) C(m-r+c-s,c-s)
```

interleavings.  The child root becomes a peak precisely when its rank is at
most `s` and it had no higher child.  This yields `A1` and `A2` with integers.
`A1` is separately checked by the rooted-tree hook formula.

For `Btop2`, fix each ordered pair of non-adjacent vertices as the two highest
peaks.  Delete them.  In every remaining component, vertices adjacent to one
of the deleted roots are marked and can never peak; the same interleaving DP,
with all unmarked peaks rejected, counts its valid relative orders.  The
component orders are then interleaved by the exact multinomial coefficient.

`run_verify.py` exhaustively checks paths and stars through eight vertices,
then checks the 30-vertex tree formed by joining roots of two height-three
complete binary trees.  It obtains

\[
 A_2=4570935802581417984000000,\qquad
 B_{\rm top2}=1382099471396359372800000,
\]

and therefore

\[
 A_2/B_{\rm top2}=838907595/253657849>3=\Delta.
\]

This is a finite counterexample to the proposed all-bounded-degree inequality
`A2 <= Delta Btop2`.  It says nothing directly about free-boundary grids or
the original twin-peaks asymptotic question.
