# Independent review of r041 joint geometry kernel

## Final verdict: pass after endpoint clarification

The former comma notation was ambiguous.  The staged proof now explicitly
assigns `+infinity` only at the two feasible endpoint pairs
`(theta,u)=(0,0),(1,0)`, while at `u=0,0<theta<1` it retains the finite value
\([2\theta(1-\theta)]^{-1}\).  This is correct.  On the closed feasible set
`u<=theta<=1-u`, `min(kappa,K)` is continuous: approaching either endpoint
forces the denominator to zero and the numerator to one.

The conditional root-swap retains the root geometry, so no orientation is
lost.  The exact identity \(\kappa=q(q-1)/(2ab)\), including `q=2` and
singleton fibres, and the resulting size-biased `(D,E)` marginal are correct.
The UI discussion is explicitly conditional.  The interior/endpoint product
conditions (16)--(17) are valid sufficient conditions: kappa is bounded on a
fixed interior interval and `Z_n>=1`, while (17) controls the remaining
normalized product.  Nothing in the package claims these conditions are
proved or that the original theorem is resolved.
