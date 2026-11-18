# Periodic stationary variance graphs and measurable minima

Let \(X\) be Haar on \(\mathbb T_L=\mathbb R/L\mathbb Z\), \(G\sim N(0,1)\) independent, and \(w\ge0\) belong to \(L^\infty(\mathbb T_L)\). The law of \(X+\sqrt{w(X)}G\pmod L\) is Haar if and only if
\[
\int_{\mathbb T_L}e^{ik_nx}e^{-k_n^2w(x)/2}dx=0,
\qquad k_n=2\pi n/L,\quad n\ne0.
\tag{1}
\]
This follows from the conditional Gaussian characteristic function and Fourier uniqueness. Subtracting the essential infimum of \(w\) preserves (1). Put \(v=w-\operatorname{ess\,inf} w\), and assume the graph is nonconstant.

## Microscopic normalization

In the marked tail parametrization, if \(t^3(f(t+z/\theta'(t))-1)\) has a bounded nonnegative limiting profile \(B\) and \(\theta'(t)\sim\lambda t^{3/2}\), the mark cutoff is \(H=t^2(f-1)=O(t^{-1})\). Writing \(p=tV\) for the rescaled exponential mark, the real shift \(-\theta'(t)V/t\) tends to zero and the imaginary shift has limiting variance \(\lambda^2p\). The corresponding spatial noise width is of order \(t^{-3/2}\) in the original coordinate. The truncated-exponential weight tends to a constant on \(0\le p\le B(z)\). If passage to the limiting marked identity is justified, a nonzero periodic mode gives
\[
\int_0^L e^{ik_nz}\int_0^{B(z)}e^{\lambda^2k_n^2p/2}dp\,dz=0,
\]
which is (1) with \(w(z)=\lambda^2(\operatorname{ess\,sup}B-B(z))\). The inner integral is \(2(e^{\lambda^2 k_n^2B/2}-1)/(\lambda^2 k_n^2)\), the constant disappears, and the common exponential factor may be removed. The variance complement accounts for the sign. This verifies the limiting coefficient separately from any justification of extending R47 Fourier control to frequencies of order \(t^{3/2}\). The results below apply to all graphs satisfying (1).

## Uniform measurable exclusion of flat intervals

For \(s>0\), let \(P_s(y)=\sum_{\ell\in\mathbb Z}p_s(y+\ell L)\) be the wrapped heat density with respect to Lebesgue measure. Since \(v\) is nonconstant and nonnegative, some band \(E=\{x:v_0\le v(x)\le V_0\}\) with \(v_0>0\) has positive measure. Its conditional wrapped densities are bounded below by
\[
a=(2\pi V_0)^{-1/2}e^{-L^2/(8v_0)}>0,
\]
by choosing an argument in \([-L/2,L/2]\) for one wrap. This band alone contributes at least \(\kappa\) times Haar to each Borel target, where \(0<\kappa<1\) can be chosen with \(\kappa\le |E|a\). In particular the minimum set \(\{v=0\}\) has measure zero, since the zero sources already contribute its full Haar mass there.

Let \(I_r(a_0)\) be the circular interval of radius \(r\) centered at \(a_0\). Choose \(\epsilon>0\) small enough to make the right side below less than 1. For small \(r\), the set \(I_r(a_0)\cap\{v\le\epsilon^4r^2\}\) is disjoint from \(E\) and its output remains within \(I_{(1+\epsilon)r}(a_0)\) whenever \(|G|\le1/\epsilon\). Adding the fixed bulk contribution gives
\[
\frac{|I_r(a_0)\cap\{v\le\epsilon^4r^2\}|}{2r}
\le\frac{(1-\kappa)(1+\epsilon)}{\mathbb P(|G|\le1/\epsilon)}=1-\delta<1.
\tag{2}
\]
This bound is uniform in the center. It excludes even a sequence of essentially flat intervals with moving centers and requires no pointwise representative of \(v\).

## One-sided quadratic averages are impossible

Adding \(\zeta>0\) to the variance preserves stationarity. With unnormalized Lebesgue input the resulting continuous wrapped density is 1 at every target. At any fixed target \(a_0\), Fatou as \(\zeta\downarrow0\) and one wrap give
\[
\int_{\{0<|u|<L/4,\ \eta\le v(a_0+u)/u^2\le C\}}\frac{du}{|u|}<\infty
\qquad(0<\eta<C<\infty),\tag{3}
\]
because on this set the Gaussian density at \(-u\) is at least \((2\pi C)^{-1/2}e^{-1/(2\eta)}/|u|\). Consequently intermediate ratios in any fixed such range occupy \(o(r)\) length in \([r,2r]\), and similarly on the other side. View \([r,2r]\) as an interval in (2) with radius \(r/2\). On at least a fraction \(\delta\), its variance exceeds \(\epsilon^4 r^2/4\), hence \(v(a_0+u)/u^2>\epsilon^4/16\). Removing intermediate ratios using (3), for every fixed \(C>0\) we obtain
\[
\liminf_{r\downarrow0}\frac{|\{u\in[r,2r]:v(a_0+u)>C u^2\}|}{r}\ge\delta.\tag{4}
\]
Reflection gives the same bound on the other side. In particular \(r^{-3}\int_r^{2r}v(a_0\pm u)du\to\infty\), since its limit inferior is at least \(\delta C\) for every \(C\). Thus a nonconstant stationary graph can have neither an essential one-sided quadratic upper contact nor bounded one-sided quadratic averages along any sequence of shrinking annuli. No ordinary continuity or differentiability is required.

## Finitely many isolated minima with monotone wings

Suppose the essential minima have finitely many locations \(a_i\), with disjoint fixed small neighborhoods, a positive essential variance gap outside, and nondecreasing outward variance on each local wing. Monotone representatives suffice. The infimum approached among these wings is assumed to be zero. The one-sided average divergence then forces \(v(a_i\pm u)/u^2\to\infty\) along each wing, since its value at \(u\) bounds the average over \([u/2,u]\).

All local negative variance moments must be finite at those points. For each \(j\ge1\), differentiate the smoothed wrapped output at \(a_i\) to order \(2j\). The principal wrap has sign \((-1)^j\) and magnitude at least \(c_j(v+\zeta)^{-j-1/2}\) in a sufficiently small source neighborhood, because \(|u|/\sqrt{v+\zeta}\) is uniformly small and \(He_{2j}(0)=(-1)^j(2j-1)!!\). The principal derivative outside that neighborhood and the other wrap derivatives are uniformly bounded in \(\zeta\) by the Gaussian bound \(C_j|u-\ell L|^{-2j-1}\) at separated arguments, summable over \(\ell\ne0\). The output derivative is zero. Monotone convergence as \(\zeta\downarrow0\) gives \(\int v^{-j-1/2}<\infty\) locally for every \(j\).

Let \(b(r)\) be the smallest outward wing value at radius \(r\), chosen between the one-sided limits of the monotone representatives. For small \(r\) the variance outside all radius-\(r\) minimum neighborhoods is essentially at least \(b(r)\). Choose positive integers \(n_j=2^j+O(1)\) for which all phases \(2\pi n_j(a_i-a_1)/L\) lie in a fixed small neighborhood of zero modulo \(2\pi\). These simultaneous returns have bounded gaps, since the orbit closure of the phase vector is compact and finitely many translates of the open neighborhood of zero by orbit points cover it. Put \(k_j=2\pi n_j/L\), \(r_j=c/k_j\) for small fixed \(c>0\). The real phase of \(e^{ik_j(x-a_1)}\) has a positive lower bound throughout the union of radius-\(r_j\) neighborhoods. At least one wing of length \(r_{j+1}\) has variance at most \(b(r_{j+1})\). All possibly negative contribution to the real Fourier identity lies outside the radius-\(r_j\) neighborhoods. Thus
\[
c_0 r_{j+1}e^{-k_j^2b(r_{j+1})/2}\le L e^{-k_j^2b(r_j)/2}.
\]
It follows that \(b(r_j)-b(r_{j+1})=O(k_j^{-2}\log k_j)\); telescoping gives \(b(r_j)=O(r_j^2\log(1/r_j))\), using \(b(r_j)\to0\). At least one wing then has a negative variance moment over length \(r_j\) of size at least \(r_j b(r_j)^{-p}\), which diverges for \(p>1/2\). This contradicts the finite sum of the corresponding wing integrals. In particular a continuous nonconstant stationary graph cannot have finitely many monotone pieces. A positive-length minimum plateau is already excluded, and otherwise finitely many monotone pieces give finitely many isolated minima with the required positive gap and monotone wings.

## Local verification

The companion script uses SymPy to check the marked exponential antiderivative and the Gaussian Fourier multiplier. It calibrates wrapped Fourier coefficients at three variances and three frequencies with mpmath at 50 and 75 digits. The analytic bounds above prove the rigidity statements.
