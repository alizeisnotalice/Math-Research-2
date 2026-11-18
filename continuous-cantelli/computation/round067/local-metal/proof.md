# Localization and Fourier coefficients of the periodic heat model

Write $g$ for a positive integer and $a_j=2\pi j/g$ for its valley locations. Let $Q>0$ be continuous and $2\pi$-periodic. Write $p>0,d>0$, put $\rho_g(\theta)=2|\sin(g\theta/2)|/g$, and set $v=d\rho_g^p Q$. Uniform input phase followed by Gaussian noise of variance $v$, taken modulo $2\pi$, produces
$$
F_n=\frac1{2\pi}\int_0^{2\pi} e^{in\theta-n^2v(\theta)/2}\,d\theta,\qquad n\ge1.
$$
Adding a common positive variance $v_*$ multiplies each coefficient by $e^{-n^2v_*/2}$. It preserves all vanishing questions.

## Valley scaling

Partition the circle into cells $a_j+[-\pi/g,\pi/g]$. Let $s_n=(2/(dn^2))^{1/p}$, and substitute $\theta=a_j+s_n z$. Because $\rho_g(a_j+s_n z)/s_n\to|z|$ and $Q(a_j+s_n z)\to Q_j=Q(a_j)$, the rescaled real heat weight on each cell tends to $e^{-Q_j|z|^p}$. Extend each weight by zero off the rescaled cell. The sine bound
$$
\rho_g(a_j+u)\ge 2|u|/\pi,\qquad |u|\le\pi/g,
$$
dominates the extended integrands by $e^{-c|z|^p}$ for some $c>0$, so the convergence is in $L^1(\mathbb R)$.

For $p<2$, $ns_n\to0$. Along a sequence with fixed residue $n\bmod g$, dominated convergence gives
$$
n^{2/p}F_n\longrightarrow
\frac{\Gamma(1+1/p)}{\pi}\left(\frac2d\right)^{1/p}
     \sum_{j=0}^{g-1}e^{2\pi i n j/g} Q_j^{-1/p}.
$$
For $p=2$, $ns_n=\sqrt{2/d}$. The Gaussian Fourier integral instead gives
$$
n F_n\longrightarrow
\frac1{\sqrt{2\pi d}}\sum_{j=0}^{g-1}
      e^{2\pi i n j/g}Q_j^{-1/2}\exp\left[-\frac1{2d Q_j}\right].
$$
For $p>2$, $ns_n\to\infty$. Replacing the rescaled heat weight by its $L^1$ limit costs $o(1)$ uniformly in the Fourier phase. The Riemann–Lebesgue lemma for each limiting weight gives $n^{2/p}F_n\to0$. This argument does not quantify the rate for superquadratic valleys.

## Vanishing and positivity

If $Q$ is $2\pi/g$-periodic, the heat profile has that period. The cell contributions differ only by their roots of unity. Thus $F_n=0$ exactly for $g\nmid n$.

None of the prescribed profiles yields a uniform output. For $p\le2$, restrict to $g\mid n$. The corresponding limiting sum above is strictly positive, so some nonzero Fourier coefficient survives. For $p>2$, use the continuous density after adding independent heat noise of variance $\varepsilon$. At a valley $a_j$, on $|u|\le K\sqrt\varepsilon$ the variance $v(a_j+u)=o(u^2)$ is $o(\varepsilon)$ uniformly for fixed $K$. The centered term of the periodized Gaussian kernel on this arc converges, after rescaling, to the standard Gaussian integral over $[-K,K]$. Its density contribution tends to $ (2\pi)^{-1}\mathbb P(|Z|\le K)$. A distinct closed interval where $v>0$ adds a fixed positive amount, uniformly for small $\varepsilon$, by strict positivity and continuity of the heat kernel on that interval. Choose $K$ large. The smoothed density at the valley exceeds $1/(2\pi)$ for small $\varepsilon$, incompatible with the preserved constant density of a uniform output. This is a statement about the prescribed circle model.

## Small symmetry breaking

The computation uses $Q=Q_0+\eta B$, with $Q_0$ of period $2\pi/g$, $Q_0\ge0.6$, $|B|\le1$, and $0\le\eta\le1/4$. Let $s(\theta)=dn^2\rho_g(\theta)^p/2$. For $g\nmid n$, symmetry gives the stable identity
$$
\frac{F_n(\eta)}{\eta}=
\frac1{2\pi}\int e^{in\theta-sQ_0}
       \frac{e^{-s\eta B}-1}{\eta}\,d\theta .
$$
The derivative at $\eta=0$ is $L_n=-(2\pi)^{-1}\int e^{in\theta-sQ_0}sB\,d\theta$. The second derivative in $\eta$ is bounded in modulus by $\sup_{s\ge0}s^2e^{-0.35s}=4e^{-2}/0.35^2$. Taylor's theorem gives $|F_n(\eta)/\eta-L_n| \le 2e^{-2}\eta/0.35^2<(9/4)\eta$. The real-arithmetic estimate is uniform in $n$.

## Discarded tails and polynomial subtraction

The local integration endpoint is the smaller of $\pi/g$ and
$$
\frac{\pi}{2}\left(\frac2{dn^2}\right)^{1/p}
       \left(\frac{37}{0.35}\right)^{1/p}.
$$
By the sine inequality, on any excluded region $s\ge37/0.35$. The unit-normalized raw integrand there is bounded by $e^{-37}$. The scaled difference and $L_n$ are bounded by $s e^{-0.35s}\le(37/0.35)e^{-37}$ after the mean-value theorem. The integration measure has total mass at most one. These omitted contributions are below $10^{-16}$ and $10^{-14}$, respectively.

When $|-s\eta B|<0.1$, the GPU replaces $e^{x}-1$ by its degree-six polynomial. The Taylor error divided by $\eta$, multiplied by $e^{-sQ_0}$, is at most
$$
\frac{10}{9}\frac{10^{-6}}{7!}s e^{-0.6s}.
$$
Indeed $e^{|x|}\le10/9$ and $|x|^7/\eta\le10^{-6}s$. The supremum of this expression is below $1.4\times10^{-10}$. The Fraction calculation verifies the numerical real-arithmetic bounds. Binary32 quadrature and rounding are checked separately against independent numerical implementations.
