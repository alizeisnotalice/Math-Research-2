# Compensating flat exceptional sources in the fractional time budget

Assume the setting of R193. Thus \(X\) has standard Gaussian density \(g\), \(v\ge0\) is bounded continuous, and for all \(s>0\)
\[
 \int g(x)p_{v(x)+s}(a-x)dx=p_{H+s}(a),\qquad H>0.
\]
Let \(\mu\) be a probability supported on a compact Ahlfors \(d\)-regular subset \(E\) of \(\{v=0\}\), \(0<d<1\), and let \(Q=1-d\), \(j=Q/2\), \(r(x)=\operatorname{dist}(x,E)\). Write \(h_b=P_b\mu\), and set
\[
 F(b,x)=\int_0^\infty(h_b-h_{b+t})(x)t^{-1-\kappa}dt,
 \qquad 0<\kappa<1.
\]

**Theorem.** Fix \(K>0\) and choose \(\kappa>0\) sufficiently small for the R193 positivity bound for this \(K\) and regularity data. Partition a punctured collar \(0<r<\delta\) into measurable sets \(G,B\), with \(r\le K\sqrt v\) on \(G\). Suppose
\[
 \int_B \min\{1,v(x)/r(x)^2\}\,r(x)^{-Q-2\kappa}\,dx<\infty. \tag{1}
\]
For the exceptional set define its baseline heat occupancy by
\[
 J_B(s)=\int_B g(x)h_s(x)dx.
\]
One sufficient condition on it is that, for some \(\alpha>\kappa\),
\[
 J_B(s)\le J_*+C s^\alpha\quad(0<s<s_0),\qquad
 J_*:=\limsup_{s\downarrow0}J_B(s).                            \tag{2}
\]
Then
\[
 \int_G v(x)^{-j-\kappa}dx<\infty.                             \tag{3}
\]
Condition (2) can also be replaced by a uniform lower bound along some sequence \(s_n\downarrow0\) for \(\int_0^\infty(J_B(s_n)-J_B(s_n+t))t^{-1-\kappa}dt\).

**Compensation estimate.** The Ahlfors upper bound and Gaussian estimates, summed over annuli of width \(\sqrt b\), give
\[
 |\partial_b^\ell h_b(x)|\le C_\ell b^{-j-\ell}e^{-r(x)^2/(8b)}
 \quad (b>0,\ \ell=0,1,2).
\]
For \(r>0\), optimizing in \(b\) gives respectively \(Cr^{-Q},Cr^{-Q-2},Cr^{-Q-4}\). Split the defining difference at \(t=r^2\). For smaller \(t\), use the age derivative and for larger \(t\) use the size bound. This gives \(|F(b,x)|\le C_\kappa r^{-Q-2\kappa}\) for every \(b>0\). Applying the same split to the first age derivative of the difference, with the second derivative on the short interval, gives
\[
 |\partial_b F(b,x)|\le C_\kappa r^{-Q-2\kappa-2}.
\]
The split estimates justify differentiation by dominated convergence at every \(b>0\). Combining these bounds,
\[
 |F(v+s,x)-F(s,x)|
 \le C_\kappa\min\{1,v/r^2\}\,r^{-Q-2\kappa}.                \tag{4}
\]
In particular the signed integral over \(B\) differs by a uniformly bounded amount from the baseline Marchaud difference of \(J_B\), by (1) and Fubini. For fixed \(s>0\), short increments in the Fubini exchange are bounded using the age derivative and long ones by the kernel size.

Here is why (2) suffices. Reduce \(\alpha\) if needed to at most one, preserving \(\alpha>\kappa\). Choose \(s_n\) maximizing \(J_B(s)-(C+1)s^\alpha\) on \([\epsilon_n,\eta]\), where \(0<\eta<s_0\), \(\epsilon_n\downarrow0\). This function is continuous for \(s>0\), has upper limit \(J_*\) at zero, and is bounded above by \(J_*-s^\alpha\). The maximum values tend to \(J_*\), so \(s_n\downarrow0\) after passing to a subsequence. For \(0<t\le\eta-s_n\), maximality gives
\[
 J_B(s_n+t)-J_B(s_n)
 \le (C+1)((s_n+t)^\alpha-s_n^\alpha)\le(C+1)t^\alpha.
\]
For larger \(t\), the boundedness \(|J_B|\le\|g\|_\infty\) suffices. The baseline differences are therefore uniformly bounded below along this sequence.

**Output estimate and the good portion.** Exactness and the same Fubini argument give
\[
 \int g(x)F(v(x)+s,x)dx
 =\int d\mu(a)\int_0^\infty
       [p_{H+s}(a)-p_{H+s+t}(a)]t^{-1-\kappa}dt.
\]
This is uniformly bounded for small positive \(s\) because \(H>0\). On \(G\) sufficiently near \(E\), R193 bounds \(F(v+s,x)\) below by \(c(v+s)^{-j-\kappa}\); the smoothing preserves \(r\le K\sqrt{v+s}\) and \(v+s\) is uniformly small near \(E\). Any remaining part of \(G\) with distance bounded below has uniformly bounded integral of \(|F|\), as does the complement of the collar. Its contribution to (3) is also finite by \(v\ge r^2/K^2\) on \(G\). The regular set \(E\) is Lebesgue null. Subtract the uniform lower bound on the integral over \(B\), use the uniform bound on the full output integral, and apply Fatou along the sequence from (2). The Gaussian density is bounded below on the compact collar. This proves (3).

## Finite-interval patching and a geometric extension

Suppose the flat phase is \(B=D\cap\{0<r<\delta\}\), where \(D\) is a finite union of intervals, allowing rays. It may cut through \(E\), and the geometry need not be symmetric. The baseline occupancy extends to
\[
 J_B(s)=\int_E g(a)1_D(a)\mu(da)+O(s^{d/2}+s).
\]
Indeed the replacement of \(B\) by \(D\) changes its baseline occupancy by \(O(e^{-\delta^2/(Cs)})\). For each \(a\in E\), the difference between \(P_s(g1_D)(a)\) and \(1_D(a)P_sg(a)\) is bounded by \(C e^{-\operatorname{dist}(a,\partial D)^2/(Cs)}\), by the Gaussian crossing bound. Ahlfors upper regularity at the finitely many boundary points (also for off-set centers) and a dyadic annulus sum bound its \(\mu\) average by \(Cs^{d/2}\). Finally \(P_sg=p_{1+s}=g+O(s)\) uniformly. More generally the same argument gives \(O(s^{\beta/2}+s)\) under \(\mu(\{a:\operatorname{dist}(a,\partial D)\le u\})\le C u^\beta\), \(\beta>0\).

For example, suppose \(v\le C r^{2+\theta}\) on \(B\) for some \(\theta>0\), and \(v\ge K^{-2}r^2\) on \(G\). Then (1) follows from the upper tube volume \(|\{r<u\}|\le C u^Q\) if \(2\kappa<\theta\). Take also \(2\kappa<d\) for finite-interval patching, or \(2\kappa<\beta\) for the geometric extension. The conclusion can exclude near-quadratic opening on \(G\) even when \(\int_B r^{-Q-2\kappa}dx=\infty\). For instance if on infinitely many disjoint geometric shells \(u/a<r<u\) the good phase contains a subset of length \(c u^Q\) with \(v\le C u^p\), where \(p(j+\kappa)\ge Q\), then each such shell consumes a positive amount of (3), a contradiction. A fixed good spatial region containing a nontrivial portion of \(E\) strictly inside it supplies the good tube lower bound by the same ball-packing estimate as R193.

Another example takes reflection symmetry of \(E,\mu\) about zero and lets \(B\) select exactly one point of almost every reflected pair in the collar. Then \(J_B(s)\) equals half the full collar occupancy for every \(s>0\), even for a rough selection. The full collar has occupancy \(\int_E p_{1+s}d\mu+O(e^{-\delta^2/(Cs)})\), so (2) holds. The flatness assumption above gives (1), and the original distance-weighted cost of \(B\) diverges since it has half the volume of every shell.

## Pointwise limitation

For every off-set \(x\), as \(b\downarrow0\) the fractional test tends to
\[
 F(0,x)=-\frac{2^\kappa\Gamma(\kappa+1/2)}{\sqrt\pi}
             \int |x-a|^{-1-2\kappa}\mu(da)
       \asymp -r^{-Q-2\kappa}.
\]
The identity is obtained by the substitution \(u=|x-a|^2/(2t)\), after dominated convergence using the split estimate above. The two mass bounds yield the stated order by a nearest-point ball and annulus sum. The derivative bound gives \(F(b,x)-F(0,x)=O(b r^{-Q-2\kappa-2})\); hence when \(b/r^2\) is sufficiently small the negative singularity persists. The gain in (4) thus requires a baseline occupancy estimate and cannot come from a pointwise smaller negative kernel on very flat bad sources.

## Local calculation

Run `/Users/zzh/.venvs/main/bin/python check_compensation.py`. Python 3.13.18 and SymPy 1.14.0 verify the exact zero-age change of variables with zero residual. At 60 decimal digits, mpmath 1.3.0 checks the zero-age constant and the first age correction for a point source at unit distance and \(\kappa=1/4\). The normalized first age correction is \(1.0000150207\) at \(b=2^{-18}\), consistent with the linear vanishing of the compensated error. A 64-atom reflection-symmetric Cantor approximation gives left-half baseline occupancy equal to half the full-line occupancy at all four tested ages. The quadratures illustrate the kernel calculation; the compact Ahlfors-set theorem follows from the preceding analytic estimates.
