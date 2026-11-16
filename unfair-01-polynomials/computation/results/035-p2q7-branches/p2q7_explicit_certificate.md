# Explicit four-offset certificate: `T={0,1}`, `ell=4`, `j=7`

This is one branch of the finite propagation ledger, not a statement about
the five branches still open under the same limited rules.

Here `p=2`, `q=7`, and
\[
\operatorname{supp}G=\{0,1,7,8\}.
\]
The selected even index is `2ell=8`, so the first-negative equations force
positive `F` coefficients at offsets
\[
f_{r-2},\qquad f_{r-28},\qquad f_{r-7}.
\]
Multiplication by the indicated positive `G` coefficients gives four
positive `H` coefficients below `r`:
\[
\begin{array}{c|c}
\text{factor product}&\text{positive coefficient of }H\\ \hline
f_{r-2}g_0&H_{r-2}\\
f_{r-7}g_1&H_{r-6}\\
f_{r-7}g_0&H_{r-7}\\
f_{r-28}g_7&H_{r-21}.
\end{array}
\]
The first pair of indices has separation `4=2*2`; the second has separation
`14=7*2`.  Hence two distinct positive terms occur at the same coefficient
of `C=H(x^2)H(x^7)`:
\[
2(r-7)+7(r-6)=2(r-21)+7(r-2).
\]
All four displayed `H` indices are strictly below the first negative index
`r`, so cancellation by `f_r` is impossible.  Thus this branch contradicts
binary `C`.
