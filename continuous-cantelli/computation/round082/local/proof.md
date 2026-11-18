# Global Jensen sign and explicit compensation examples

Let $G\sim N(0,1)$, $s=t^2/2$, $B=1+\beta$, and $y=\sqrt B\,x$. Use probabilists' Hermite polynomials. Take $\beta>0$, $A\in(0,1]$, $0\le r\le1/12$, and $|g|\le\sqrt r/2$, and set
$$
b(x)=-A e^{-\beta x^2/2}
\left(\mathrm{He}_2(y)+g\,\mathrm{He}_3(y)+r\,\mathrm{He}_4(y)\right).
$$
This is a continuous bounded perturbation with $b\le1$. To check the upper bound, the polynomial after including the minus sign is
$$
P(y)=1-3r+3gy+(6r-1)y^2-gy^3-ry^4.
$$
When $r>0$, the upper bound $|g||y|^3\le r y^4/2+g^2 y^2/(2r)$ gives $P(y)\le1-3r+3gy-(3/8)y^2\le1-3r+6g^2\le1$. Here $6r\le1/2$, $g^2/(2r)\le1/8$, and the quadratic maximum is explicit. When $r=0$, $g=0$ and $P=1-y^2$. Multiplying by $Ae^{-\beta x^2/2}$ preserves $b\le A$. A positive additive constant makes the actual variance nonnegative.

Write $m(t)=\mathbb E b(G+t)$. Completing the Gaussian square changes $y$ to a unit Gaussian with mean $\xi=t/\sqrt B$, and contributes the factor $B^{-1/2}e^{-\beta t^2/(2B)}$. The Hermite moments give
$$
m(t)=-A B^{-1/2} e^{-\beta t^2/(2B)}\,\xi^2(1+g\xi+r\xi^2).
$$
If $r>0$, the last polynomial is at least $1-g^2/(4r)\ge15/16$. For $r=0$ it is identically one. Thus $m(0)=0$ and $m(t)<0$ on the entire punctured real line for every member of this family.

Put $J(t)=\mathbb E[e^{s b(G+t)}-1-s b(G+t)]$. It is nonnegative. The full Gaussian-shift equation requires
$$
\mathbb E e^{s b(G+t)}-1=s m(t)+J(t)=0.
$$
The numerical scalar is $D(t)=\log[J(t)/(-s m(t))]$. In particular, the sign of $D(t)$ determines the sign of this residual at nonzero $t$. The Jensen sign of $m(t)$ alone leaves both undercompensation and overcompensation possible.

## A strict tail threshold

Write $b_+=\sup_x b(x)>0$. The positive set of $b$ is contained in a compact interval because $P$ has negative even leading term. On $b(x)\le0$, the Jensen integrand is bounded by $s^2 b(x)^2/2$, and
$$
\mathbb E b(G+t)^2
\le C(1+|t|^8)e^{-\beta t^2/(1+2\beta)}
$$
by completing the Gaussian square again. This tail decays faster than the absolute linear term, since $\beta/(1+2\beta)>\beta/[2(1+\beta)]$. On the positive compact set the Jensen integrand is at most $e^{s b_+}$, giving an integral bounded by $C e^{-(1-b_+)t^2/2+O(|t|)}$. Thus if $b_+<1/(1+\beta)$, then $J/(-s m)\to0$ on both tails.

If $b_+>1/(1+\beta)$, a fixed compact interval has $b(x)$ strictly above a number $c>1/(1+\beta)$. For large $|t|$, the Jensen integrand there is at least $\tfrac12 e^{s c}$. The Gaussian mass of the interval gives a lower bound $c' e^{-(1-c)t^2/2-O(|t|)}$. Therefore $J/(-s m)\to+\infty$ on both tails. In the even family $g=0$ the exact peak is $b_+=A(1-3r)$. At equality $b_+=1/(1+\beta)$ the same ratio tends to zero for these even profiles: the positive branch's unique nondegenerate maximum is at zero, so after completing the square the maximum neighborhood contributes at most $C |t|^{-1}e^{-\beta t^2/(2B)}$. The rest of the positive set has a fixed peak gap and is exponentially smaller. The absolute negative linear term is at least $c\, t^4 e^{-\beta t^2/(2B)}$. The scan treats near-critical finite shifts separately from these limits.

## Two exact sign examples for all real shifts

Take $\beta=8,g=r=0$. Write $A$ for the amplitude. A direct square completion gives
$$
\mathbb E b(G+t)^2=
\frac{A^2}{\sqrt{17}}e^{-8t^2/17}(c_0+c_2t^2+c_4t^4),
\quad
c_0=\frac{226}{289},\quad
c_2=\left(\frac{54}{17}-2\right)\frac9{17^2},\quad
c_4=\frac{81}{17^4}.
$$
For $A=1/2$, the small-shift limit is
$$
\lim_{t\to0}\frac{J(t)}{-s m(t)}
=\frac{6102}{2312\sqrt{17}}<1.
$$
The strict supercritical tail bound gives $J/(-s m)>1$ eventually. Continuity supplies a nonzero balanced shift. The existence argument uses these analytic strict signs. The reported location from quadrature is numerical.

For $A=1/64$, use $J/s^2\le\tfrac12e^{sA}\mathbb E b(G+t)^2$. The ratio is bounded by
$$
\frac{27 A}{4\sqrt{17}} e^{-\kappa t^2}(c_0+c_2t^2+c_4t^4),
\qquad \kappa=4/153-A/2>0.
$$
For $x\ge0$, $xe^{-\kappa x}\le1/(e\kappa)$ and $x^2e^{-\kappa x}\le4/(e\kappa)^2$. Use $\sqrt{17}>4$ and $e>8/3$. The rational upper bound is below $1/10$, so the residual is strictly negative at every nonzero $t$. Both examples satisfy the same strict global Jensen sign on $m(t)$.

## Analytic error budget for the computational model

The stable integrand is $b^2\phi_2(sb)$, with
$$
\phi_2(z)=\int_0^1(1-u)e^{uz}\,du=\frac{e^z-1-z}{z^2}.
$$
Extending through zero gives $\phi_2(0)=1/2$. For $|z|\le0.3$, the degree-six polynomial's absolute error is at most $(10/7)(0.3)^7/9!<9\times10^{-10}$. This follows by the Taylor remainder in the integral. Its relative error is below $3\times10^{-9}$, using $\phi_2(z)\ge0.35$.

For the scan box $\beta\in[1/64,512]$, its polynomial Gaussian gives $|b|<2000$. Indeed $\beta/(2B)\ge1/130$, $|g|<1/6$, and the monomial suprema after weighting are bounded respectively by $6,65,1000,130^2$ for powers one through four. On $|G|>24$, the omitted $J/s^2$ is bounded by $2000^2 e^{72-288}$, via $|t|\le12$, $b_+\le1$, and the Gaussian tail. Divide by the explicit $-m(t)/s$, which is at least $2\cdot2^{-10}(15/16)e^{-72}/12000$ in the scan box. The omitted amount in the ratio is below $10^{-48}$. These are real-arithmetic truncation bounds; GPU quadrature and rounding are verified empirically on selected inputs.
