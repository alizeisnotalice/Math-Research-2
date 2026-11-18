# Weighted Gaussian center areas on a positive tail

Let \(f:\mathbb R\to[0,M]\) be continuous and satisfy
\[
 F_f(t)=\int e^{tx-t^2f(x)/2}\gamma(dx)=1,\qquad t\in\mathbb R,
\]
and write \(\nu_t(dx)=e^{tx-t^2f(x)/2}\gamma(dx)\).
The bounds below use no tail derivative assumptions and do not require the global supremum to be attained.

## Exact area formula and estimates uniform over tail sectors

The zero set of \(f\) has measure zero. For example, on \(\{x\ge a,f(x)=0\}\), with \(a>0\), the identity gives an upper bound \(e^{-ta}\) on its Gaussian measure as \(t\to+\infty\). Reflection covers the negative half.
Where \(f(x)>0\), put
\[
 m(x)=x/f(x),\qquad \sigma(x)^2=1/f(x),\qquad
 d\rho(x)=f(x)^{-1/2}e^{x(m(x)-x)/2}dx.
\]
Completing the square and applying Tonelli shows that for any Borel interval \(J\) and source set \(E\),
\[
 \int_J \nu_t(E)\,dt
 =\int_E \mathbb P\{N(m(x),\sigma(x)^2)\in J\}\,d\rho(x).
                                                        \tag{1}
\]
In particular the right side for \(E=\mathbb R\) equals \(|J|\).
For \(x\ge R>0\), one has \(\sigma^2=m/x\le m/R\). Set \(p=\mathbb P(|G|\le1)>1/3\).
For \(0\le a<b\),
\[
 \rho\{x\ge R,\ m(x)\in[a,b]\}
 \le p^{-1}\bigl[(b-a)+2\sqrt{b/R}\bigr].                \tag{2}
\]
Indeed, for such \(x\) the normal in (1) lies in
\([a-\sqrt{b/R},b+\sqrt{b/R}]\) with probability at least \(p\).

A quantitative area correspondence follows. For \(RT\ge1\), \(J=[a,b)\subset[T,2T]\), and any source set \(E\subseteq[R,\infty)\),
\[
 \left|\rho\{x\in E:m(x)\in J\}-\int_J\nu_t(E)dt\right|
 \le 640\sqrt{T/R}.                                    \tag{3}
\]
Here \(R\) and \(E\) may vary with \(T\).
To prove (3), let \(w=\sqrt{T/R}\), \(\kappa=\sqrt{RT}\). Integrate the absolute difference of the indicator and the Gaussian window probability in (1).

- For \(T/2\le m\le4T\) the variance bound is \(4w^2\). The classification error is at most \(\sum_{y\in\{a,b\}}\exp(-(m-y)^2/(8w^2))\). On each strip \(nw\le|m-y|<(n+1)w\) there are at most two center intervals of length \(w\), each of \(\rho\)-mass at most \(5w/p\) by (2). Thus this part is at most \(20w(1+\sqrt{2\pi})/p<80w/p\).
- For \(m<T/2\), the crossing probability is at most \(e^{-RT/4}\). Bound (2) gives at most \(w(\kappa/2+\sqrt2)e^{-\kappa^2/4}/p<2w/p\).
- For \(2^jT\le m<2^{j+1}T\), \(j\ge2\), the probability of reaching \(J\) is at most \(e^{-2^j\kappa^2/8}\). The strip contribution is at most \(w(\kappa2^j+2\sqrt{2^{j+1}})e^{-2^j\kappa^2/8}/p\). Use \(\kappa e^{-\kappa^2 2^j/16}<1\) and \(\sum_{j\ge2}2^j e^{-2^j/16}\le 2\int_2^\infty e^{-s/16}ds<32\) to bound the sum by \(128w/p\).

Adding proves (3) since \(p>1/3\).
In particular, if \(\nu_t(x<R)\to0\), the tail center measure on intervals of length comparable to \(T\) is asymptotically Lebesgue. For a non-vanishing compact contribution, (3) measures the residual area
\[
 |J|-\int_J\nu_t(x<R)dt
\]
up to the same error. This retains the compact contribution rather than assuming it vanishes.

## Low and high sources in the area formula

For \(0<\eta<1\), the source subset \(x\ge R\), \(m\in[T,2T]\), \(f\le\eta\) has
\(\omega(x)=d\rho/dx\ge\eta^{-1/2}e^{(1-\eta)RT/2}\).
Its ordinary Lebesgue measure is at most
\[
 { \sqrt\eta \over p}\left(T+2\sqrt{2T/R}\right)e^{-(1-\eta)RT/2}.
                                                        \tag{4}
\]
This bound does not require that the low source subset carry little \(\rho\)-mass.

One may also sum the geometry over the whole positive tail. Suppose \(R\ge\sqrt\eta\), and set
\(q=e^{-(1-\eta)R/(2\sqrt\eta)}\).
On \(x\ge R,f(x)\le\eta\) we have \(m\ge R/\eta\).
Use center slabs of length \(\eta^{-1/2}\) starting at \(R/\eta\).
For slab \(j\ge0\), (2) bounds its mass by \((2j+5)/(p\sqrt\eta)\).
The increasing exponential lower bound on \(\omega\) then gives
\[
 |\{x\ge R:f(x)\le\eta\}|
 \le {e^{-(1-\eta)R^2/(2\eta)}\over p}
       \left({5\over 1-q}+{2q\over(1-q)^2}\right).
                                                        \tag{5}
\]
For example this is at most \(42e^{-(1-\eta)R^2/(2\eta)}\) if
\(R\ge\max(\sqrt\eta,2\sqrt\eta\log2/(1-\eta))\).
In particular an entire positive tail cannot approach a limit \(<1\), including a limit of zero. Individual low wells remain possible under this necessary condition.

At a high source \(f\ge1+\delta\), \(\omega(x)\le(1+\delta)^{-1/2}
e^{-\delta x^2/[2(1+\delta)]}\).
Also \(m\ge a>0\) implies \(x\ge (1+\delta)a\). Integration gives
\[
 \rho\{f\ge1+\delta,\ m\ge a\}
 \le {e^{-\delta(1+\delta)a^2/2}\over
            \delta a\sqrt{1+\delta}}.                 \tag{6}
\]
Thus high sources cannot provide order-\(T\) residual area on center intervals in \([T,2T]\). This conclusion still allows a fixed compact set or low wells to supply area.

For a more geometric valley bound, let \(f\ge1+\delta\) on \([u,v]\), take \(\eta>0\), and let
\(J=[a,b)\subset[T,2T]\cap[(1+\eta)u,v/(1+\delta)]\).
A source with \(m\in J\) and \(f<1+\delta\) must satisfy \(x<u\), since \(x\in[u,v]\) is high by hypothesis and \(x>v\) implies \(f=x/m>1+\delta\). It then has \(f\le1/(1+\eta)\). Consequently,
\[
 \begin{split}
 &\rho\{R\le x<u,\ m\in J,\ f\le1/(1+\eta)\}\\
 &\quad\ge |J|-\int_J\nu_t(x<R)dt
     -640\sqrt{T/R}
     -{e^{-\delta(1+\delta)a^2/2}\over\delta a\sqrt{1+\delta}}.
 \end{split}                                           \tag{7}
\]
This is a cumulative lower bound for left-side low wells on a positive plateau.

## Precise area deficit when the positive tail stays at or above one

Assume \(f\ge1\) on \(x\ge R>0\). Put
\[
 B(x)=1-f(x)^{-1/2}e^{-x^2(1-1/f(x))/2},\qquad x\ge R,
 \qquad c_t=\nu_t(x<R).
\]
Now \(\sigma\le1\). The strip expansion used for (3), with unit strips and no exterior split, gives \(O(1)\) discrepancy in (1), uniformly for all half-open center intervals. Since \(m\le x\), changing from the window \([a,b)\) in \(m\) to the same window in \(x\) changes its \(\rho\)-mass by at most \(2/a+2/b\) for \(b\ge a\ge R\): at either boundary \(y\), crossing sources have \(x\ge y\) and \(m\le y\), so their total mass is at most
\[
 \int_y^\infty e^{-x(x-y)/2}dx\le 2/y.
\]
It follows that
\[
 \int_a^b B(x)\,dx=\int_a^b c_t\,dt+O(1)              \tag{8}
\]
uniformly in the interval. For \(f>0\) everywhere, \(c_t\) is integrable at \(+\infty\), by the split at any fixed negative location and positivity on a compact interval. Consequently \(B\) is integrable. If \(f(x)=cx^2+o(x^2)\) at zero, with \(c>0\) and no positive zeros elsewhere, the compact asymptotic is \(c_t\sim A/t\), \(A=e^{1/(2c)}/\sqrt c\). This follows as in round 25 by the substitution \(s=tx\) at the origin and dominated convergence. In this case (8) gives \(\int_R^T B(x)dx\sim A\log T\).

If \(c_t\to0\), the discrepancy in (8) is in fact \(o(1)\) as \(a\to\infty\), uniformly in \(b\ge a\). To see this, the signed error for a single boundary \(y\) is
\[
 L_y=\int_{x\ge R}\left[{\bf1}_{m\ge y}
                 -\Phi((m-y)/\sigma)\right]d\rho(x).
\]
For fixed \(H>1\), the part \(|m-y|>H\) is uniformly bounded by a vanishing Gaussian strip tail as \(H\to\infty\), using \(\sigma\le1\) and (2) with enlargement by 1 instead of \(\sqrt{b/R}\). On \(|m-y|\le H\), sources \(x\ge y+H+1\) have \(\rho\)-mass tending to zero as \(y\to\infty\) by \(\omega(x)\le e^{-x(x-y-H)/2}\). The remaining sources lie in \([y-H,y+H+1]\).

The real identity and the eventual sign give
\[
 \int_{x\ge R}(1-e^{-y^2(f(x)-1)/2})p_1(x-y)dx\le c_y.
\]
On this fixed window and for every fixed \(\varepsilon>0\), the set where \(y^2(f-1)>\varepsilon\) thus has ordinary measure \(o(1)\). Its \(\rho\)-measure is also \(o(1)\) since \(\omega\le1\). On its complement one has \(|m-x|=O(\varepsilon/y)\), \(\sigma=1+O(\varepsilon/y^2)\), and \(\omega=1+O_H(\varepsilon)\). After discarding \(O(\varepsilon/y)\)-width neighborhoods of \(y-H,y,y+H\), the remaining indicator comparisons are uniform. Hence, for fixed \(H\), the truncated \(L_y\) differs by at most \(o(1)+O_H(\varepsilon)\) from the standard \((m,\sigma,\omega)=(x,1,1)\) truncated integral on \([y-H,y+H]\), which vanishes by oddness. Let \(\varepsilon\downarrow0\), then \(H\uparrow\infty\). This proves \(L_y\to0\). Taking differences of the two boundary errors and using the \(2/a+2/b\) spatial-boundary estimate proves the refinement to (8).

For a quadratic origin and eventually \(f\ge1\), the refinement gives the exact dyadic area law
\[
 \int_T^{2T} B(x)\,dx \longrightarrow
 {e^{1/(2c)}\over\sqrt c}\log2.                         \tag{9}
\]
It allows arbitrary continuous high excursions on the positive tail. If a smooth asymptotic \(f(x)-1\sim Kx^{-3}\) happens to hold as well, the left side is asymptotic to \((K/2)\log2\), forcing \(K=2A\).

## Local checks

The companion Python script checks the square completion and shell constants in SymPy 1.14.0 and Fraction. mpmath 1.3.0 at 80 decimal digits uses \(f=1\) to check (3) via the exact truncated-normal antiderivative at six combinations of the spatial cutoff and parameter window. Largest observed absolute error divided by \(\sqrt{T/R}\) is about 0.56419. The proven bound \(640\) rests on the deterministic estimates rather than this calibration. Python 3.13.18 ran the script to completion.
