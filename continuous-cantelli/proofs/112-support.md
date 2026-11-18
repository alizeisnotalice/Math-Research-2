# Opposite support bounds force a strict inversion and reflected moment budgets

Let \(f\ge0\) be bounded and continuous with exact normalized identity \(\int e^{zx-z^2f(x)/2}\gamma(dx)=1\), and suppose \(f(0)=0\). Write \(g=p_1\) and \(Q_t=p_t(x)dx-\int_{f\le t}P_{t-f(x)}\delta_x\,\gamma(dx)\). No positivity of \(Q_t\) is assumed. Suppose on the right \(H_+=\liminf_{x\to+\infty}f(x)>1\). Fix \(0<t\le H_+-1\), and let
\[
 B=\sup\{x:f(x)<t\},\qquad A=\inf\{x:f(x)>t\}.
\]
Here \(0<B<\infty\), by continuity at 0 and the right tail gap; \(A\) may be \(-\infty\). The upper support bound below is at the strict past endpoint \(B\), not at the inclusive sublevel endpoint.

## Strict support and no instantaneous outside layer

For every \(s<t\), the one-sided resolvent support theorem of R087 gives \(\operatorname{supp}Q_s\subset(-\infty,\max\{x:f(x)\le s\}]\subset(-\infty,B]\). As \(s\uparrow t\), weak convergence gives
\[
\operatorname{supp}Q_{t-}\subset(-\infty,B],\qquad
Q_{t-}=p_t(y)dy-\int_{f(x)<t}P_{t-f(x)}\delta_x\gamma(dx).
\]
The future transform for this strict-past remainder is
\[
\int e^{\lambda y}Q_{t-}(dy)=\int_{f(x)\ge t}e^{\lambda x-(f(x)-t)\lambda^2/2}\gamma(dx),\qquad \lambda>0.
\]
The left side is at most \(2e^{\lambda B}\) in absolute value. The right side is positive. If \(\gamma\{f=t,x>B\}>0\), some set of positive measure at distance \(d>0\) beyond \(B\) contributes at least \(c e^{\lambda(B+d)}\), a contradiction. Hence the zero-age layer at time \(t\) has no measure beyond \(B\). Subtracting it shows \(\operatorname{supp}Q_t\subset(-\infty,B]\) as well. This applies also at \(t=H_+-1\), because only earlier times were used, and it permits a measure-zero record birth at \(t\) to change the inclusive geometric hull.

Round 58 gives \(\operatorname{supp}Q_t\subset[A,\infty)\) from the strict future support. Thus \(Q_t\) is supported on \([A,B]\). It is absolutely continuous: the strictly positive-age part has a density and the zero-age layer is a restriction of \(\gamma\). Since its total mass is \(m_t=\gamma(f>t)>0\), if \(A\) is finite then \(A<B\). In particular there exist strict inversion points \(x<y\) with \(f(x)>t>f(y)\) whenever \(A\) is finite. The inequality \(A<B\) uses the strict future and strict past simultaneously and cannot be attained by a level atom.

## Reflected moment hierarchy

Write \(R_B(y)=p_t(y)-p_t(2B-y)\) on \(y<B\). For any strict-past mean \(x\le B\), the Gaussian kernel at \(y<B\) is at least its value at \(2B-y>B\), because \(|y-x|\le|2B-y-x|\). At the latter point the strict-past potential equals \(p_t(2B-y)\) by the strict-past support. This is a pointwise equality there, because all means are separated and the off-diagonal Gaussian bound gives continuity. Hence \(Q_{t-}\le R_B(y)dy\) on \(y<B\), and subtracting the positive level layer gives the same upper bound on \(Q_t\). Since \(B>0\), \(R_B\ge0\) on that half-line. It follows for every polynomial \(P\) nonnegative on \([A,B]\) (or on \((−\infty,B]\) when \(A=-\infty\)) that
\[
\int_{f>t}\big[\exp(-\tfrac12(f(x)-t)\partial^2)P\big](x)\,\gamma(dx)
 \ \le\ \int_A^B P(y)R_B(y)\,dy. \tag{1}
\]
The exponential on a polynomial is a finite sum. Its left side equals \(\int P\,dQ_t\) by differentiating the exact future transform at zero (bounded \(f\) and Gaussian integrability justify this). Equivalently the reflected moment sequence minus the future moment sequence is the moment sequence of a nonnegative measure supported in the corridor. This is an upper hierarchy in the signed regime; no positivity of \(Q_t\) is inferred.

For finite \(A\), put \(d=B-A\). Since \(R_B(y)\le K(B-y)\) for \(y\le B\), with \(K=2/(\sqrt{2\pi e}\,t)\), tests \(P=1\), \(P(y)=y-A\), and \(P(y)=(y-A)(B-y)\) give respectively
\[
\begin{split}
 m_t &\le 2\Phi(B/\sqrt t)-\Phi(A/\sqrt t)-\Phi((2B-A)/\sqrt t)\le Kd^2/2,\\
 \int_{f>t}(x-A)\,d\gamma &\le Kd^3/6,\\
 \int_{f>t}[(x-A)(B-x)+(f(x)-t)]\,d\gamma &\le Kd^4/12.
\end{split}
\]
These inequalities give quantitative necessary corridor widths. In particular a jump inward of the future endpoint cannot close the overlap with strict past. A point birth beyond \(B\) does not count toward the room needed at that time.

There is a related plateau constraint. On the part of \(\{f=t\}\) lying to the left of \(A\), where \(Q_t=0\), one has \(Q_{t-}(dy)=g(y)dy\) restricted to that level. Hence \(g(y)\le R_B(y)\) almost everywhere on that portion. This excludes a positive-length disappearing plateau in any sufficiently small left collar of \(B\), where \(R_B\) tends to zero but \(g\) stays positive.

## Global consequences and local dying maxima

If positive-side lower early sublevels are intervals for all times \(s<H_+-1\), then necessarily
\[
 \sup_{x\le0}f(x)>H_+-1.
\]
Indeed otherwise, at \(t=H_+-1\), \(A\ge0\). Earlier connectivity makes the positive strict past an interval up to \(B\), so no strict inversion can exist on the positive side; this contradicts the crossed endpoints. By the one-sided version of the R087 birth estimate, this connectivity follows if every positive record minimum below \(H_+-1\) satisfies its one-flank logarithmic regularity condition. The estimate needs only the high right tail to give right exterior vanishing. No left tail regularity or interior geometry is required. In the continuous single-well case with finite unequal tail limits, it excludes a difference of at least one between those two limits, without analytic matching at the critical radius.

Separately, if a point \(a\) with \(f(a)=t>0\) lies strictly outside the convex hull of \(\{f>t\}\), then a neighborhood of \(a\) contains no future points and \(Q_t=0\) there. The local past output equals \(p_t\), with variance deficit \(v(x)=t-f(x)\) and full input density \(g\) in that neighborhood. A dying separated future component at an outward endpoint jump has precisely such a local maximum. The usual contact obstructions localize: all omitted input centers are separated and their positive-age Gaussian potentials are holomorphic on a small disc about \(a\); zero-age inputs outside contribute no local mass. For example a nonconstant real-analytic germ at \(a\) is impossible by the R023 core Fourier argument, which uses only analyticity of the core density at the contact. Likewise a one-sided comparable deficit \(c|u|^\alpha\le v(a+u)\le C|u|^\alpha\) with \(1<\alpha<2\) is excluded by the R068 opposite-wing monotonicity argument. If a flank is \(o(u^2)\) even in rescaled measure, it saturates \(g(a)\) in the corresponding one-sided small output interval; any separated positive strict-past band adds strictly positive density. Therefore \(p_t(a)>g(a)\) when strict past is nonempty. For \(t>1\) it is nonempty by \(\mathbb E f=1\), giving \(|a|>\sqrt{t\log t/(t-1)}\) at any such dying flat maximum. These conclusions require an exterior gap in the strict future and do not constrain arbitrary interior maxima between two surviving future components.

## Local verification

`check_constraints.py` verifies the three polynomial budgets and the backward-heat quadratic term exactly using SymPy. It checks the reflected Gaussian CDF formula and its quadratic upper bound at three widths at 50 and 75 digits using mpmath. The support and comparison proofs above are analytic.
