# Pointwise and finite-moment inverse estimates

Let $X,G$ be independent standard Gaussian variables, let $a\ge0$ be bounded, and set $m=E_\gamma a$, $b=a-m$, $A=\|b\|_\infty$. Exactness means $X+\sqrt{a(X)}G\sim N(0,1+m)$. A spatial modulus for $b$ is denoted by $\omega$. All estimates below use only values of $b$, not derivatives.

## Pointwise heat inversion

Let $P_\tau$ be Gaussian convolution of variance $\tau$. For any $\tau>A$, exactness gives
$$
P_\tau(b\gamma)(y)=
-\frac12\partial_y^2
\int_0^1(1-s)\int b(x)^2 p_{\tau+s b(x)}(y-x)\,\gamma(dx)\,ds .
$$
In fact, cancellation of the common noise factor in the characteristic functions yields
$$
0=\int e^{i\xi x}
\left[-\frac{\xi^2 b(x)}{2}+
\frac{\xi^4 b(x)^2}{4}\int_0^1(1-s)e^{-s\xi^2 b(x)/2}\,ds\right]\gamma(dx).
$$
Multiply by the smooth heat factor after dividing by $-\xi^2/2$, then invert. The apparent division at $\xi=0$ extends by $E_\gamma b=0$. All heat variances in the inverse expression have a positive common lower bound.

For $\tau\ge2A$ and $q\in[\tau/2,3\tau/2]$, putting $t=z^2/\tau$ gives
$$
|p_q''(z)|\le \frac4{\tau}(1+2t)e^{-t/12}p_{2\tau}(z)
\le \frac{48}{\tau}p_{2\tau}(z).
$$
It follows that $|P_\tau(b\gamma)(y)|\le(12/\tau)P_{2\tau}(b^2\gamma)(y)\le(12 A^2/\tau)p_{1+2\tau}(y)$. Completing the square at $y=(1+\tau)x$ gives the bound for every $x$,
$$
 |b(x)|\le \frac{12 A^2}{\tau}e^{\tau x^2/2}
 +E\omega\!\left(\sqrt{\frac{\tau}{1+\tau}}\,|G|\right),
 \qquad \tau\ge2A>0 .
$$
For $0<\varepsilon\le1$, take $\tau=48A/\varepsilon$. At every point with $|x|\le\sqrt{2\log2/\tau}$ and $|b(x)|\ge\varepsilon A$, the Gaussian average of the modulus is at least $\varepsilon A/2$. Using the actual modulus bounded by $2A$ and cutting off a Gaussian tail shows $\omega(C_\varepsilon\sqrt A)\ge\varepsilon A/4$. This is a pointwise restriction in a window of radius proportional to $A^{-1/2}$, without any lower bound on $\|b\|_2/A$.

## Finite moment residual and Hermite inverse

This part allows a non-exact graph. Let $h_n=\mathrm{He}_n/\sqrt{n!}$ be the orthonormal Gaussian Hermite basis, $\Pi_N$ the projection onto orders at most $N$, and $\sigma=\|b\|_2$. Define
$$
R_j=[t^j]\left(E_\gamma e^{tX-t^2/2+t^2b(X)/2}-1\right),\qquad
E_N=\left(\sum_{n=0}^N n!|R_{n+2}|^2\right)^{1/2}.
$$
These are normalized moment residuals for the Gaussian output target. Coefficient extraction gives
$$
\langle b,h_n\rangle
=2\sqrt{n!}R_{n+2}-
\sum_{\substack{j\ge2\\2j\le n+2}}
\frac{2^{1-j}}{j!}
\sqrt{\frac{n!}{(n-2j+2)!}}\,
\langle b^j,h_{n-2j+2}\rangle .
$$
For each $j$ the map from the coefficients of $b^j$ to those of order at most $N$ is a coordinate shift followed by a diagonal factor bounded by $N^{j-1}$. Taking the Euclidean norm before summing over $j$ avoids an extra dimension factor. Since $\|b^j\|_2\le A^{j-1}\sigma$, writing $x=AN/2$ gives
$$
\|\Pi_N b\|_2\le 2E_N+\eta(x)\sigma,\qquad
\eta(x)=\frac{e^x-1-x}{x},\quad\eta(0)=0 .
$$
For the stationary Ornstein–Uhlenbeck pair $(X,X_h)$ of correlation $e^{-h}$, each Hermite mode contributes exactly $2(1-e^{-nh})|\langle b,h_n\rangle|^2$ to the squared increment. Its spatial displacement has law $\sqrt{2(1-e^{-h})}\,G$. Hence with
$$
J_h=E\omega\!\left(\sqrt{2(1-e^{-h})}|G|\right)^2
$$
one obtains the finite witness
$$
2E_N\ \ge\
\left[\sqrt{\left(\sigma^2-\frac{J_h}{2(1-e^{-(N+1)h})}\right)_+}
-\eta(AN/2)\sigma\right]_+ .
$$
For exact profiles with $0<A\le1$, take $N=\lfloor1/A\rfloor$ and $h=(\log2)/N$. Since $\eta(1/2)<0.3$, at least $0.91\sigma^2$ lies at Hermite orders above $N$ and $J_h\ge0.91\sigma^2$. Thus nonzero residuals $b$ of exact near-constant profiles cannot have $b/\sigma$ in an $L^2(\gamma)$-precompact family as $A\to0$. The conditions involve relative shape or a near-maximal point in the expanding window. They do not by themselves settle global rigidity for arbitrary continuous profiles.
