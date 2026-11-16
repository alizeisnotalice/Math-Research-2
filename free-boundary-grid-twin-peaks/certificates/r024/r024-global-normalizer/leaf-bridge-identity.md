# A leaf-bridge identity and the remaining bridge problem

Let \(H\) have \(n\) vertices, fix \(u\in V(H)\), and form \(G\) by
adding a leaf \(y\) joined only to \(u\).  Write a rank order of \(H\) as
\(\pi\), let \(K(\pi)\) be its number of peaks, let \(p(\pi)\) be the
rank position of \(u\), and write \(I_u(\pi)\) for the indicator that
\(u\) is a peak.  Smaller rank means a larger label.

## Exact insertion identity

Each labelling of \(G\) is uniquely a labelling \(\pi\) of \(H\), plus a
rank insertion position \(r\in\{1,\ldots,n+1\}\) for \(y\).  If
\(r>p\), the new leaf is not a peak and no old peak changes.  If
\(r\le p\), the leaf is a peak and it destroys the peak at \(u\) exactly
when \(u\) had been a peak.  Therefore

\[
\begin{aligned}
A_2(G)={}&(n+1)\#\{\pi:K=2,I_u=1\}\\
&+\sum_{\pi:K=2,I_u=0}(n+1-p(\pi))\\
&+\sum_{\pi:K=1,I_u=0}p(\pi). \tag{L}
\end{aligned}
\]

This is a direct bijective identity, with no probabilistic approximation.

## Why it is the relevant stress test

The last term in (L) contains a first moment of a designated vertex's
rank under one-peak orders.  On a connected graph it can be of order \(n\)
for selected roots/vertices.  Merely adding one bridge edge does not make a
pointwise delay estimate valid.

The compensating terms can also be displayed exactly.  Let
\(\mathcal B_H\) be the set of `B`-orders of \(H\), and use `#` and sums
only over rank orders of \(H\).  Splitting according to whether the new
leaf has rank 1, rank 2, or a later rank gives

\[
\begin{aligned}
B(G)={}&(n-1)\#\{\pi\in\mathcal B_H:I_u=1\}\\
&+\sum_{\pi\in\mathcal B_H:I_u=0}(n+1-p(\pi))\\
&+2\#\{K=1,I_u=0\}\\
&+2\#\{K=2,\ u\text{ is a peak but not the rank-1 peak}\}. \tag{B}
\end{aligned}
\]

The first two lines are the cases where the new leaf occurs after the old
top two ranks.  The last two lines are the two choices in which the leaf has
rank 1 or 2.  Formula (B) was independently checked by direct permutation
enumeration on a path and two non-isomorphic four-vertex examples.

Thus (L) versus (B) identifies the exact missing comparisons.  A
counterexample driven specifically by a large one-peak rank moment needs
the \(\sum_{K=1,I_u=0}p\) term to dominate the terms in (B).  This is not a
necessary mechanism for every bridge counterexample: either of the first
two terms in (L) could also be anomalously large relative to (B).  A
positive global charge must control all three terms of (L), not merely its
last term.  This is the exact point where the present unbalanced-bridge
attempt is incomplete.

For orientation only, the targeted family consisting of a star with one
leaf extended by one extra edge has exact ratios
\(A_2/B=1.3333,1.4211,\ldots,2.1713\) for 2 through 15 original star
leaves.  It is not bounded-degree and is not evidence for a limit.  Its
role is solely to check the insertion identity against a family with large
one-peak rank moments.

## Averaging over every possible attachment vertex

There is a useful exact aggregate identity.  For each \(u\in V(H)\), let
\(G_u\) be obtained by attaching one leaf at \(u\).  Set
\(t=\binom{n+1}{2}\), and let

\[
 T_2(H)=\sum_{\pi:K(\pi)=2}\tau(\pi),
\]

where \(\tau\) is the rank of the lower of the two peaks.  Summing (L) and
(B) over all attachment vertices gives

\[
\boxed{\ \sum_u A_2(G_u)=(t-1)A_1(H)+(t+1)A_2(H)+T_2(H)\ }, \tag{A}
\]
\[
\boxed{\ \sum_u B(G_u)=(t-1)B(H)+2A_2(H)+2(n-1)A_1(H)\ }. \tag{C}
\]

For (A), a one-peak order contributes the rank sum
\(2+\cdots+n=t-1\).  In a two-peak order, summing its three terms in (L)
over all vertices gives \(t+1+\tau\).  For (C), a `B`-order contributes
\(t-1\) after its attachment vertices are summed; the last two terms of
(B) contribute respectively \(2(n-1)A_1\) and \(2A_2\).

Using the second-birth fibre representation, the established marked-sum
identity is

\[
 T_2(H)=\frac12\sum_{\sigma\in B(H)}\bigl(L(\sigma)^2+3L(\sigma)\bigr).
\]

Thus uniform attachment displays a second \(L\)-moment in an exact
identity, while the denominator has only first-order terms.  This fact alone
does not say that the second moment needs independent control.  The next
comparison records precisely what a hypothetical uniform normalizer bound
would imply.

## What a uniform normalizer bound actually yields

Suppose, only for this paragraph, that every one-leaf attachment obeys
\(A_2(G_u)\le C'B(G_u)\), where \(C'\) is independent of \(n\).  Put
\(R=A_2(H)/B(H)=\mathbb E_B L\), \(\alpha=A_1(H)/B(H)\), and
\(M_2=\mathbb E_B L^2\).  Substitution of (A), (C), and the marked-sum
identity gives the exact implied upper bound

\[
\begin{aligned}
 M_2\le{}&2C'\bigl(t-1+2R\bigr)-\bigl(2(t+1)+3\bigr)R\\
 &+2\alpha\bigl(2C'(n-1)-(t-1)\bigr). \tag{U}
\end{aligned}
\]

This should not be overstated.  Independently of any attachment argument,
\(1\le L\le n-2\) gives
\[
 M_2\le(n-2)R. \tag{T}
\]
If the proposed normalizer bound also gives \(R\le C'\), (T) is \(O(n)\),
whereas the leading part of (U) is
\(2t(C'-R)\), generally \(O(n^2)\).  Therefore (U) is weaker than the
trivial bound unless the gap \(C'-R\) is at most order \(1/n\), or the
negative \(\alpha\)-term makes a separately verified improvement.  No such
near-saturation or \(\alpha\) estimate is proved here.

This neither preserves a bounded maximum degree nor supplies an estimate
for square grids.  It is only an exact reduction that prevents treating the
appearance of \(L^2\) as a new necessary obstruction.
