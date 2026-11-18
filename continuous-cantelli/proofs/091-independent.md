# Stefan uniqueness, the first opposite fork, and flat germs

Let $f$ be an exact normalized graph. Suppose the true early exit on $0<t<T$ is already established as in rounds 83–84. Let $g=p_1$. The surviving probability density is $q(t,x)$ and the stopped law has density $g1_{I(t)}$ on the expanding bounded interval. The density of the law stopped at time $t$ or earlier is $H=q+g1_I$.

## Uniqueness and symmetry

The pointwise identity $q=(H-g)_+$ uses only support and positivity. For two candidate early exits, set $K=H_1-H_2$, $Q=q_1-q_2$, and $W=\int_0^t Q(s)\,ds$. Itô's formula gives $W_{xx}=2K$ in distributions and $W_t=Q$. The monotone one-Lipschitz map $h\mapsto(h-g)_+$ gives $KQ\ge Q^2$.

On every compact early slab the two differences are in $L^2(dt\,dx)$, including time zero, since $\|p_t\|_2^2=O(t^{-1/2})$ and both absorbed profiles are bounded on a fixed compact interval. With spatial Fourier cutoff $|\xi|\le N$,
$$
\frac12\|\partial_xW_N(t)\|_2^2
=-2\int_0^t\langle K_N,Q_N\rangle\,ds.
$$
Indeed the spatially truncated functions have the required derivatives in $L^2$ and are absolutely continuous in time. The right side tends to the same untruncated pairing. Its nonnegativity forces $Q=0$ almost everywhere, then $K=0$. Weak continuity of the stopped laws and $g>0$ identify the endpoints at every early time. Reflection is another candidate, so the early boundary is symmetric. Its limiting radius is denoted by $b$; thus $f(x)=f(-x)$ on $[0,b]$ and their common value at $b$ is $T$.

For a finite grid, conditional on the current Brownian position $y$, the earlier position has a centered-Gaussian translate with mean proportional to $y$. Multiplication by an even centered interval indicator and Gaussian convolution preserve radial nonincrease. The bridge probability of survival to the observation time has that property by passing to dense grids with strict margins on $[\varepsilon,t]$. At $T$ one can condition on $B_T$ directly and use finite grids strictly before $T$. This avoids requiring pointwise convergence of density ratios at the terminal time.

## Endpoint bound and separated births

The survivor at $T$, equivalently the strict-past remainder $Q_{T-}$, is supported on $[-b,b]$. Its path never exceeded $b$ before $T$, so by reflection its density is bounded by $p_T(y)-p_T(2b-y)$ for $y<b$. In particular its one-sided Laplace transform $M(z)=\int e^{zy}Q_{T-}(dy)$ satisfies
$$
M(z)\le Ce^{bz}/z^2,\qquad
\limsup_{z\to\infty}z^2e^{-bz}M(z)\le (2b/T)p_T(b).
$$
The precise constant follows from the linear upper cap at $y=b$ and the change of variables $b-y=s/z$. The future transform for the strict-past remainder includes the whole layer $f=T$:
$$
M(z)=\int_{f\ge T}e^{zx-z^2(f(x)-T)/2}\gamma(dx).
$$
For any region where $x\ge b+D$ and $T\le f\le T+h$, its Gaussian mass is bounded by $C z^{-2}e^{-Dz+hz^2/2}$. Choosing $z=D/h$ for fixed $D>0$ gives $C h^2 D^{-2}e^{-D^2/(2h)}$. Consequently at a remote point $x_0>b$ with $f(x_0)=T$, the same one-flank modulus criterion as in round 87 excludes the separated birth. Reflection proves the left statement. This uses the endpoint survivor directly and does not require continuing the Stefan solution or left/right resolvent support past $T$. If the early no-birth modulus was imposed only at levels strictly below $T$, an additional hypothesis at remote $T$-level minima is required here.

Taking $D=0$ and $z=\sqrt{2/h}$ gives $\gamma\{x\ge b,T\le f\le T+h\}=O(h)$. For a nondecreasing high wing, this gives a pointwise lower bound $f(b+u)-T\ge c_0u$. If the outward derivative at $b$ exists, no local monotonicity is needed for this latter conclusion. For a finite derivative $c$, positivity of the future integral and scaling $x=b+s/z^2$ give $\liminf z^2e^{-bz}M(z)\ge2g(b)/c$, with infinity understood if $c=0$. Thus necessarily $c\ge Tg(b)/(b p_T(b))>0$.

## First fork and derivative-free direction test

Suppose the two tail limits are finite, $H_+>H_->1$, and $T=H_--1$. Assume there are no separated $T$-level births on the high wing, for example by the endpoint estimate and the stated modulus. Then $f(x)>T$ for every $x>b$, giving a strict gap on each compact subinterval away from $b$. For $z>0$ use the exact odd real transform
$$
0=2\int_0^\infty g(x)\sinh(zx)
\left[e^{-z^2 f(x)/2}-e^{-z^2 f(-x)/2}\right]dx.
$$
If $f(x)\ge f(-x)$ on $(b,b+\varepsilon)$, all positive contributions lie in an intermediate compact interval on the high wing with $f\ge T+\delta$. Their upper bound is $C e^{-(T+\delta)z^2/2+O(z)}$. Farther out the bracket is negative by the unequal limits. On $x\in[z,z+1]$, completing the Gaussian square gives a negative contribution of magnitude at least $c e^{-Tz^2/2-o(z^2)}$. It cannot be canceled by the compact positive part. Thus every right neighborhood of $b$ contains $x$ with $f(x)<f(-x)$, opposite to the far-tail ordering. The earliest asymmetric radius is $b$ and the infimum asymmetric level is $T$. This assertion gives arbitrarily near points of the indicated sign, not a fixed-sign interval.

## Quasianalytic and Gevrey thresholds

If the two branches are smooth through $b$ in a common reflection- and subtraction-stable quasianalytic germ class, the difference has every Taylor coefficient zero at $b$ by equality to its left. Taylor-map injectivity makes it zero on a neighborhood of $b$, contradicting the immediate fork. Neither convergence of the Taylor series nor analytic continuation of it is needed. Punctured one-sided analyticity with all derivatives extending to $b$ is insufficient: a nonzero smooth flat germ is analytic away from the joining point. A uniform quasianalytic germ bound at the joining point is essential.

Under the stronger low-wing tail bound $f(-x)\le H_-+C/x$ eventually, the far negative term has magnitude at least $c e^{-Tz^2/2-Cz/2}$ up to fixed factors. In particular $C=0$ applies to a nondecreasing low wing. If the high wing has $f(b+u)-T\ge c_0u$ and $\Delta(u)=(f(-b-u)-f(b+u))_+$ were eventually at most $e^{-\kappa/u}$, then the positive near contribution, after normalization by $e^{Tz^2/2}$, would be bounded by
$$
C' z^2 e^{bz}\int_0^\varepsilon e^{-(cz^2/2-z)u-\kappa/u}\,du,\qquad 0<c<c_0.
$$
The integral is at most $\varepsilon e^{-2\sqrt{(cz^2/2-z)\kappa}}$. This cannot offset the far contribution if $\sqrt{2c\kappa}>b+C/2$. Hence
$$
\liminf_{u\downarrow0}u\log(1/\Delta(u))
\le (b+C/2)^2/(2c_0).
$$
Under equal left jets, a uniform Gevrey germ bound $|D^nh|\le M R^n(n!)^s$ for the difference gives, by Taylor remainder optimization, flat decay $e^{-c u^{-1/(s-1)}}$ for $1<s<2$. The preceding inequality rules it out. For $s=2$, the same optimization gives $e^{-\kappa/u}$ for any $\kappa<1/R$. Mere existence of the unequal tail limits gives only an $e^{-o(z^2)}$ far lower bound after normalization, sufficient for the qualitative fork and quasianalytic exclusion. The quantitative Gevrey argument needs the stated low-tail control or a comparable bound.

## Product form of the fork–tail tradeoff

The pointwise slope can be avoided in the tradeoff itself. With $h(u)=f(b+u)-T>0$ and the same positive difference $\Delta$, put
$$
\kappa=\liminf_{u\downarrow0}h(u)\log(1/\Delta(u)).
$$
For a nondecreasing low wing let $\theta=\limsup_{z\to\infty}z(H_--f(-z))$. For any $k<\kappa$ with $k\ge0$, one eventually has $\Delta(u)\le e^{-k/h(u)}$. The normalized positive near integral is then at most
$$
C_\varepsilon z^2\exp[(b+\varepsilon)z-z\sqrt{2k}]
$$
by minimizing $z^2h/2+k/h$ over $h>0$, without a derivative or lower bound on $h/u$. Along a subsequence realizing the limsup, use $x\in[z-1,z]$ and monotonicity of the low wing to bound the normalized negative far term below by $c\exp[z^2(H_--f(-z))/2]$. The intervening compact positive terms are exponentially smaller on the quadratic scale. Let first $z\to\infty$, then $\varepsilon\downarrow0$ and $k\uparrow\kappa$. This gives
$$
\sqrt{2\kappa}+\theta/2\le b.
$$
If $\kappa=0$, the same comparison with the elementary bound $\Delta\le O(1)$ gives $\theta/2\le b$. If either quantity is positive infinite, fix arbitrarily large finite subthresholds to obtain a contradiction. For a nonmonotone low wing the argument instead uses the limsup of $(2/z)$ times the logarithm of its actual normalized Gaussian-window contribution. An eventual bound $f(-x)\le H_-+C/x$ ensures that quantity is at least $-C$.
