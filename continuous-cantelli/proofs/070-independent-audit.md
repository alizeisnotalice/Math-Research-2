# Independent check of the double-log flux coefficients

Use the flux identity from round 63, with
$$
F(d)=d^{-1/2}e^{-1/(2d)},\quad F(0)=0,\qquad
J(s,d)=-F(d+e^{-2s})+1_{s\ge0}F(d).
$$
Let $L$ denote the heat logarithmic scale and $X=\log L$. On either side write the critical profile as
$$
q_\pm(t)=\frac{1}{2d_\pm(t)}
=Y+\gamma\log Y+b+r_\pm(Y),\qquad Y=\log t,\qquad r_\pm(Y)=o(1),
$$
with $\gamma=5/2$ and $b=\log(4/\sqrt\pi)$. Set $C=\sqrt{\pi/2}$, $a(Y)=\gamma\log Y+b$, and let $g_0$ be the Gaussian source density at the spatial contact point.

## Frozen prepotential

All fixed $d$ derivatives of $J(s,d)$ have an integrable exponential majorant for bounded nonnegative $d$. At $d=0$ the derivatives of $F(d)$ vanish. Consequently the coefficient of $d^n$, $n\ge1$, in the Taylor expansion of $\int J(s,d)\,ds$ is
$$
-\frac{1}{2n!}\int_0^\infty F^{(n)}(y)\frac{dy}{y}
=-C(2n-1)!! .
$$
Indeed, the first expression uses $y=e^{-2s}$. Integrating by parts $n$ times gives $n!\int F(y)y^{-n-1}\,dy$ in the numerator. The boundary terms vanish by infinite flatness at zero and polynomial decay at infinity. Substituting $y=z^{-2}$ gives the Gaussian $2n$th moment.

For $|s|\le L/2$, replace $d_\pm(L-s)$ by the common smooth model value $d_0(L)=[2(X+a(X))]^{-1}$. The difference from $d_0(L-s)$ is $o(X^{-2})$ uniformly since $r_\pm=o(1)$. The change of the model between $L-s$ and $L$ is $O(|s|/(LX^2))$. The integrable $d$ derivative of $J$ controls both replacements. The removed region $|s|>L/2$ is exponentially small, and $g_\pm(L-s)$ can be replaced by $g_0$ with the same accuracy. This calculation differentiates only the common model and the heat kernel. The normalized right side of the flux identity is therefore
$$
\frac{A(\infty)-A(L)}{g_0 C}
=2d_0(L)+6d_0(L)^2+o(X^{-2})
=\frac1X+\frac{-a(X)+3/2}{X^2}+o(X^{-2}).
$$

## Tail without derivatives of the residual

The substitution $Y=\log t$ in the left side of the flux identity gives
$$
F(d_\pm(t))\,dt
=\sqrt2e^{-b}Y^{-2}e^{-r_\pm(Y)}
 \left[1+\frac{a(Y)}{2Y}+o(Y^{-1})\right]\,dY.
$$
Here $\sqrt2e^{-b}=C/2$. The square-root expansion is uniform at infinity and needs only $r_\pm=o(1)$. The source density replacement in the tail is exponentially small in $L$. Write
$$
U(Y)=1-\tfrac12(e^{-r_+(Y)}+e^{-r_-(Y)}),\qquad
V(X)=\int_X^\infty U(Y)Y^{-2}\,dY.
$$
The normalized left side becomes
$$
\frac1X -V(X)+\int_X^\infty \frac{a(Y)}{2Y^3}\,dY
-\int_X^\infty \frac{a(Y)}{2Y}\frac{U(Y)}{Y^2}\,dY
+o(X^{-2}).
$$
Initially bounded $U$ gives $V(X)=O(\log X/X^2)$ by comparison with the right side. Now integrate the cross term by parts using $V'=-U/X^2$ almost everywhere. With $B_0(Y)=a(Y)/(2Y)$ it equals $B_0(X)V(X)+\int_X^\infty B_0'(Y)V(Y)\,dY$, which is $O(\log^2X/X^3)$. Only this integrated primitive and the explicit multiplier are differentiated. The remaining integral is
$$
\int_X^\infty \frac{a(Y)}{2Y^3}\,dY
=\frac{\gamma\log X/4+b/4+\gamma/8}{X^2}.
$$
Hence the necessary integrated relation under the mere $o(1)$ residual assumption is
$$
V(X)=\frac{(5\gamma/4)\log X+5b/4+\gamma/8-3/2+o(1)}{X^2}.
$$

If the two pointwise expansions $r_\pm(X)=[A_\pm\log X+B_\pm]/X+o(X^{-1})$ are additionally assumed, then $1-e^{-r_\pm}=r_\pm+O(r_\pm^2)$. With mean coefficients $\overline A,\overline B$ this gives $V(X)=[(\overline A/2)\log X+\overline A/4+\overline B/2+o(1)]/X^2$. Comparison yields $\overline A=5\gamma/2$ and $\overline B=5b/2-\gamma-3$.

The odd spatial heat test from round 63 forces the two pointwise coefficients to agree. If $A_+\ne A_-$, then $d_+-d_-$ has a nonzero leading coefficient at scale $(\log X)/X^3$. Otherwise a difference $B_+\ne B_-$ has scale $X^{-3}$. Both scales are larger than the discarded odd-kernel bound $O(e^{-c\sqrt X})$ on $z\ge X^{1/4}$. Within $z\le X^{1/4}$ the interpolated variances tend uniformly to 1 and scaled differences have the required bounds without taking derivatives of the remainders. The nonzero moment for fixed offset $s>0$ is $2\Phi(s)-1$. Thus the two coefficients are equal separately, giving $A_\pm=25/4$ and $B_\pm=(5/2)b-11/2$.

The integrated relation alone does not establish existence of the two pointwise expansions.
