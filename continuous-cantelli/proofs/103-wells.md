# Global Gaussian tests for separated variance wells

Let \(X,G\) be independent standard Gaussians and assume
\[
 Y=X+\sqrt{v(X)}G\sim N(0,M),\qquad 0\le v\le M .
\]
This covers the shifted conditional variance used after subtracting a common minimum. Write \(p_s\) for the centered Gaussian density of variance \(s\). Convolving the output identity with \(p_\epsilon\) gives, for every \(\epsilon>0\) and \(y\in\mathbb R\),
\[
 \int_{\mathbb R}p_1(x)p_{\epsilon+v(x)}(y-x)\,dx=p_{\epsilon+M}(y).
\]
Also \(\int v\,d\gamma=M-1\).

## Sharp constraints from interval caps

For a piecewise constant profile of allowed variance ranges \(l(x)\le v(x)\le u(x)\), integrating the pointwise minimum and maximum of \(p_{\epsilon+v}(y-x)\) gives necessary lower and upper inequalities at every \((\epsilon,y)\). The formulas are valid over the whole real line and allow any arrangement of wells and barriers. In particular the unknown background still contributes its lower envelope with \([l,u]=[0,M]\).

For a fixed displacement \(z=y-x\), the function \(T\mapsto p_T(z)\) increases up to \(T=z^2\) and then decreases. Thus on \(T\in[T_1,T_2]\) its minimum is at an endpoint. If \(T_1<T_2\), the minimum uses \(T_2\) inside and \(T_1\) outside the radius
\[
 c(T_1,T_2)
 =\sqrt{T_1T_2\log(T_2/T_1)/(T_2-T_1)}.
\]
The maximum uses \(T=\min(T_2,\max(T_1,z^2))\). A constant variance bound gives identical endpoints. These envelopes are sharp for each individual Gaussian test over measurable profiles with the pointwise bounds: choose a minimizing or maximizing endpoint or variance at each input.

For each piece and each endpoint variance, the finite integrals are explicit:
\[
 \int_a^b p_1(x)p_T(y-x)\,dx
 =p_{1+T}(y)\left[
 \Phi\left(\frac{b-y/(1+T)}{\sqrt{T/(1+T)}}\right)
 -\Phi\left(\frac{a-y/(1+T)}{\sqrt{T/(1+T)}}\right)
 \right].
\]
Split each input interval at \(y\pm c(T_1,T_2)\). This yields the whole-line lower test using only Gaussian CDFs.

If disjoint wells \(I_i\) have upper caps \(u_i\) and all other caps are \(M\), then the mean condition requires
\[
 \sum_i (M-u_i)\gamma(I_i)\le1 .
\]
When this is strict and no lower caps are imposed, the target mean is compatible with the caps alone. Smooth transitions at the endpoints can preserve it, by interpolating between a profile close to the lower choice and one close to the upper cap. Finitely many extra positive barriers can be enforced on sufficiently small neighborhoods without consuming the strict slack.

## An exact two-well certificate

Take \(\epsilon=1/4\), \(M=9/4\), \(y=0\). Prescribe \(v\le1/64\) on the two wells
\[
 [-5/8,-1/16]\quad\text{and}\quad[1/16,5/8]
\]
and retain the global bound \(v\le M\) everywhere. The lower-envelope integral divided by \(p_{\epsilon+M}(0)\) is strictly greater than \(5/4\). With only one of these wells prescribed, the lower-envelope ratio for this same Gaussian test is below \(97/100\). The cap mean budget is below \(95/100\).

The rational certificate works with the ratio integrand
\((2\pi)^{-1/2}\sqrt{(M+\epsilon)/T}\exp[-(1+1/T)x^2/2]\). Both endpoint branches are decreasing on the positive half-line at \(y=0\), so rational left and right rectangle sums enclose the integral over \([0,3]\) without having to approximate the switching point. Integer square roots and alternating Taylor bounds enclose the endpoint values. The omitted upper tail is bounded by Mills applied to the narrower endpoint branch. The strict excess remains for any admissible variance profile with both wells. The individual inequalities concern this test and its envelope.

## The sharp envelope with a fixed variance mean

The mean can be incorporated into a given Gaussian test, including a nonnegative finite mixture of such tests. Let \(K(x,v)\) denote the sum of the test kernels, and write \(\mu=M-1\). With pointwise caps as above, the exact infimum over measurable profiles with \(\int v\,d\gamma=\mu\) is
\[
 \sup_{\lambda\in\mathbb R}
 \left\{\lambda\mu+
  \int_{\mathbb R}\min_{l(x)\le v\le u(x)}(K(x,v)-\lambda v)\,d\gamma(x)\right\}
\]
for \(\int l\,d\gamma<\mu<\int u\,d\gamma\). The inequality in the direction needed for exclusion holds simply by subtracting and adding \(\lambda v\). For equality, the displayed concave function tends to minus infinity at both ends because the mean is interior. At a maximizing \(\lambda\), write \(v_-(x),v_+(x)\) for the smallest and largest pointwise minimizers, which are measurable by continuity in \(v\) and measurability of the piecewise caps. The left and right derivatives of the pointwise minimum are \(-v_-(x)\) and \(-v_+(x)\), respectively. This follows by comparing minimizers at neighboring multipliers and taking compact subsequential limits. The slopes are bounded by \(M\), so dominated convergence shows
\[
 \int v_-\,d\gamma\le\mu\le\int v_+\,d\gamma .
\]
Choose \(v_+\) on a left half-line and \(v_-\) elsewhere. The Gaussian measure has no atoms, so the mean varies continuously with this division and can be set to \(\mu\). The pointwise infimum remains attained almost everywhere. This establishes exactness over the stated measurable class. Endpoint means force \(v=l\) or \(v=u\) almost everywhere.

For a single kernel the one-dimensional tilted minimum is tractable: \(\partial_T p_T(z)=p_T(z)(z^2-T)/(2T^2)\) and the second derivative changes sign only where \(z^2/T=3\pm\sqrt6\). Thus the tilted objective can be minimized by endpoints and derivative roots on finitely many monotonicity intervals.
