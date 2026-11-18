# Slow infinite-order flatness and the first-hitting heat potential

On the positive side take \(v(u)=m_0e^{-d[(\kappa+\log(r_0/u))^p-\kappa^p]}\) for \(0<u\le r_0\), with \(p>1\), and a fixed positive fraction of the same variance on the negative side. The first hitting radius of \(m_j=m_0 2^{-j}\) is
\[
r_j=r_0e^{\kappa-L_j},\qquad
L_j=(\kappa^p+j\log 2/d)^{1/p}.
\]
This shape decays faster than every power and slower than \(e^{-d_1/u^\gamma}\) for every fixed \(d_1,\gamma>0\). Its relative gap obeys
\[
1-r_j/r_{j-1}\sim\frac{\log 2}{dp L_{j-1}^{p-1}}\asymp j^{-(p-1)/p}.
\]
The gaps are not summable, as required for the contradiction in the exact-map first-hitting argument.

For \(R=r_{j-1}\), put \(q_R=R^2/(2m_j)\) and \(A_R=dp L_{j-1}^{p-1}-2\). Write \(g=p_1\) and \(q_0=(2\pi)^{-1/2}\). Use
\[
H_v(u)=\frac{\sqrt v\,q_0 e^{-q/2}}q W(q),\qquad
q=u^2/v,\qquad
W(q)=\int_0^\infty e^{-t}(1+2t/q)^{-3/2}dt.
\]
The formula follows in \(H_v(u)=\tfrac12\int_0^v p_s(u)ds\) by setting \(s=v/(1+2t/q)\). In the positive half-core integral \(S_R^+(0)\), set
\[
u=R e^{-x},\qquad x=\frac{2t}{A_R q_R}.
\]
The log ratio of \(q(u)\) to \(q_R\) is \(d[(L_{j-1}+x)^p-L_{j-1}^p]-2x\), evaluated with \(\log1p\) and \(\operatorname{expm1}\) for small \(x\). It is at least \(A_Rx\), so the scaled integrand has an integrable bound proportional to \(e^{-t}\). Also \(q_R\to\infty\). The ratio to the endpoint Laplace term converges to one, giving
\[
S_R^+(0)\sim \frac{2g(a+R)q_0 R\sqrt{2m_j}}{A_R q_R^2}\,e^{-q_R/2}.
\]
This yields \(-\log(S_R^+(0)/g(a))/q_R\to1/2\). The other side has a weakly larger boundary action when its variance fraction is at most one.

At \(y=r_j\), set \(u=y+\sqrt{m_j}t\). The variance ratio tends to one for fixed \(t\), and it is bounded by two everywhere in the potential core. The hitting point and the outer core boundary separate by an unbounded number of thermal widths, because \((R-y)/\sqrt{m_j}\to\infty\). With \(H_2(t)\) as an integrable envelope, dominated convergence gives
\[
\frac{2S_{r_{j-1}}(r_j)}{m_jg(a+r_j)}\longrightarrow1.
\]
The integrable envelope also controls the quadrature cutoff of 24 thermal widths. If \(A_R\ge 1/\log4\), then \(q(u)/(8\log(4R/u))\) takes its minimum at \(R\). Indeed its log derivative in \(x\) is \(dp(L_{j-1}+x)^{p-1}-2-1/(\log4+x)\ge0\). Thus the endpoint action yields an explicit Cauchy coefficient cutoff \(q_R/(8\log4)\) on this model.

The two-profile smoke calculation uses 48-digit tanh–sinh and independent double adaptive integration on the same scaled variables. It compares endpoint and near-hitting normalizations across eight levels up to \(j=300\). Round-36 coordinated CPU measurements may reuse the already active 12-worker high-order leakage-jet pool. The separate slow-flat program then remains a tested auxiliary computation.
