# Last negative time of the fixed fat Cantor source

Let \(s=7/8\), \(\varepsilon=1/8\), \(I=[40001/20000,40007/20000]\), and \(\ell=|I|=3/10000\). At stage \(n\), the fixed Cantor construction removes the central relative fraction \(2^{-n}/256\) from every surviving interval. Both endpoints of \(I\) survive, so \(\operatorname{diam}K=\ell\). The removed gaps have total length at most \(\ell/256\). We treat both \(\eta=c1_K\gamma\) with \(c=15/16\) in the mixed model and \(c=1\) in the deterministic graph of round 154.

## Bridge identity and interval envelopes

For \(0<u<\varepsilon\), the only applied source is at time \(s\), and \(Q_{s+u}(dy)=[p_{s+u}(y)-(P_u\eta)(y)]dy\). Completing the Gaussian square gives
\[
 p_s(x)p_u(y-x)=p_{s+u}(y)p_v(x-m),\qquad
 v=\frac{su}{s+u},\quad m=\frac{s}{s+u}y .
\]
Thus the supremum over \(y\) of the applied-source/free-density ratio equals
\[
 F_c(v)=\sup_{m\in\mathbb R}\int_K c\frac{g(x)}{p_s(x)}p_v(x-m)\,dx .
\]
In this parametrization \(c g/p_s=c\sqrt{s}\exp(x^2/14)\), which increases across \(I\). Denote its endpoints by \(w_-\) and \(w_+\). Put \(\xi=u/\ell^2\), \(\theta=v/\ell^2=s\xi/(s+\xi\ell^2)\) and \(M(\theta)=\operatorname{erf}(1/\sqrt{8\theta})\). The Gaussian mass of \(I\) is maximized when its mean is the center \(a\) of \(I\), by differentiation or symmetric monotone rearrangement. Upper-bound by the full interval and \(w_+\). For the lower bound use the same center, replace the weight by \(w_-\), and subtract at most \(\ell/256\) times the maximum Gaussian kernel. This gives
\[
 w_-\left[M(\theta)-\frac{1}{256\sqrt{2\pi\theta}}\right]
 \le F_c(v)\le w_+M(\theta) .
\]

A strict left ratio above one gives a negative density at the center-bridge target \(y=(1+u/s)a\), and, by spatial continuity, negativity against a nonzero nonnegative smooth test on an open interval around it. A right ratio below one gives nonnegativity everywhere, with strictly positive density. For \(v>0\), the function \(m\mapsto P_v[(c g/p_s)1_K](m)\) is continuous, nonconstant and vanishes at infinity, so its positive supremum is attained. The heat semigroup and strictly positive Gaussian kernel make \(F_c(v)\) strictly decreasing as \(v\) increases. Continuity in \(v>0\) follows from \(L^1\) continuity of the Gaussian kernels. Since \(v\) also increases with \(u\), the two strict signs enclose a unique critical time \(t_*=s+u_*\). At the critical time \(Q_{t_*}\ge0\), with density zero at at least one spatial point. For \(s<t<t_*\) it has strictly negative mass on an open spatial interval. The initial slice at \(s\) is negative on \(K\), since \(w_->1\).

## Exact certificates

The certified open enclosures are

| Source | \(u_*/\ell^2\) | \(t_*=s+u_*\) |
|---|---|---|
| Mixed, \(c=15/16\) | \((0.1139039,0.1166260)\) | \((0.875000010251351,0.87500001049634)\) |
| Deterministic, \(c=1\) | \((0.1472570,0.1500079)\) | \((0.87500001325313,0.875000013500711)\) |

At the left and right endpoints the mixed margins, measured above one for the lower ratio and below one for the upper ratio, exceed \(1.08476\times10^{-8}\) and \(3.91850\times10^{-8}\). The corresponding deterministic margins exceed \(2.65808\times10^{-8}\) and \(8.27109\times10^{-8}\). These signs are rational checks. The verifier bounds \(\pi\) by Machin alternating series and uses rational integer-square-root enclosures and Taylor bounds for both endpoint exponential weights. For the Gaussian mass it uses odd and even Taylor sums of
\[
 M(\theta)=\frac{1}{\sqrt{2\pi\theta}}
 \sum_{n\ge0}\frac{(-1)^n}{(8\theta)^n(2n+1)n!} .
\]
The sum comes from integrating the Taylor formula for \(e^{-x^2/(2\theta)}\) over \([-1/2,1/2]\). The integrand's Taylor remainder sign gives the enclosures directly. An independent 78-decimal-digit check compared the Gaussian integral with `erf` and found agreement to over 67 decimals, and the geometric envelopes fell within the rational bounds.

## Later times and the Root completion

In the mixed three-layer construction the time-1 source is \(\gamma-\eta-P_\varepsilon\eta\), so \(Q_1=\eta\). For \(1<t<1+\varepsilon\), \(Q_t=P_{t-1}\eta\), and thereafter \(Q_t=0\). All of these remainders are nonnegative.

For the deterministic graph take the residual starting measure \(\mu=\gamma-P_\varepsilon\eta\) at time 1 and its regular Root completion to \(\lambda=\gamma-\eta\), as in round 154. The two measures have the same mass and first moment. Their absolute-value-potential gap is strictly positive at every finite spatial point:
\[
 \int|x-y|\,(\lambda-\mu)(dy)
 =\int\left[\mathbb E|x-y-\sqrt{\varepsilon}G|-|x-y|\right]\eta(dy)>0 .
\]
A nondegenerate Gaussian crosses both sides of the kink of the absolute value, giving the strict inequality. The regular Root obstacle has no initial contact, so there is no time-1 stop. Consequently \(Q_1=\mu\). For times beyond 1 the remainder is the nonnegative survivor law of the Root embedding started from \(\mu\): the heat evolution of \(\mu\) less the heat propagation of the completed stops, by the strong Markov property. This establishes that the pre-1 transition detected above is the last negative time in each completed model.

## Reproduction

Run `bash run_both.sh` with `certify_both.py` in the configured MacBook project directory. The calculation uses Python 3.12.13 and mpmath 1.3.0. Three bounded processes independently check both brackets at rational Taylor orders 20, 28 and 36. The batch ran 10:02:54.377572–10:02:54.400646 UTC with 0.065562 cumulative worker CPU seconds. Every worker and its supervisor exited zero, and a follow-up process query found none of their PIDs. All 16 hashes in the MacBook transfer ledger matched. Full exact endpoints and decimal displays are in `both/36.json` and `both-summary.json`.
