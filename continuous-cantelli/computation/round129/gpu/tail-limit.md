# An exact obstruction for fixed-sign regularly decaying wings

Let \(F_f(t)=E_\gamma e^{tX-t^2 f(X)/2}\), with \(0\le f\le H\). On each wing assume
\[
f(\sigma r)-1=a_\sigma (R/r)^{p_\sigma}(1-R/r)\qquad(r\ge R),
\]
where \(\sigma\) is \(+1\) or \(-1\), \(p_\sigma>0\), and \(a_\sigma\ne0\). The scanned profiles have \(f\ge1/2\). More generally the conclusions below hold whenever \(\inf_{[-K,K]}f>0\) for each finite \(K\).

For \(T\to\infty\), on the favorable wing put \(p=p_\sigma\), \(a=a_\sigma\). Then
\[
\begin{aligned}
p<2&:\quad T^{p-2}\log F_f(\sigma T)\longrightarrow -aR^p/2,\\
p=2&:\quad F_f(\sigma T)\longrightarrow e^{-aR^2/2},\\
p>2&:\quad T^{p-2}(F_f(\sigma T)-1)\longrightarrow -aR^p/2.
\end{aligned}
\]
The limiting signs in the sampled noncompact family depend only on the chosen favorable tail coefficient. Every nonzero coefficient in this family obstructs exact normalization, even for \(p>2\) where the residual tends to zero.

For the proof, complete the Gaussian square and put \(X=\sigma(T+Z)\), with \(Z\sim N(0,1)\) under the shifted Gaussian reference. Fix \(\epsilon\in(0,1/2)\). The contribution outside \(|Z|\le\epsilon T\) is exponentially small in \(T^2\). To see this, restrict first to \(|X|\le K\), where \(f\) has a positive minimum \(m_K\), and use the original integrand to bound the contribution by \(e^{TK-m_KT^2/2}\). On \(|X|>K\) choose \(\sup|f-1|\le\epsilon^2/2\) and bound the shifted Gaussian tail there by \(2e^{-\epsilon^2T^2/2}\), multiplied by \(e^{\epsilon^2T^2/4}\). On \(|Z|\le\epsilon T\), write
\[
T^p(f(X)-1)=aR^p(1+Z/T)^{-p}(1-R/(T+Z)).
\]
For \(p<2\) its upper and lower bounds uniformly approach \(aR^p\) when \(\epsilon\downarrow0\), giving the limit of the scaled logarithm. For \(p=2\), convergence in probability and boundedness on this window yield the stated limit of the integral. For \(p>2\), expand the last exponential to first order uniformly on the same window and take the Gaussian expectation. The same reasoning applies under reflection with the other prescribed coefficient.

The compact core can be arbitrary subject to the local positive lower bound. The conclusion is specific to wings with a nonzero one-term asymptotic deviation. Other tail geometries and compensation mechanisms remain open.
