# Outer monotone levels of the signed intermediate heat mass

Let \(g=p_1\) be the standard Gaussian density. Suppose \(f\) is bounded, continuous, and nonnegative, with the exact identity
\[
 \int e^{zx-z^2 f(x)/2}g(x)dx=1,\qquad z\in\mathbb C.
\]
Let \(H=\sup f\) and
\[
 Q_t=P_t\delta_0-\int_{f(x)\le t}P_{t-f(x)}\delta_x\,g(x)dx.
\]
Throughout we use the support theorem from round 58, with no positivity assumption on \(Q_t\). Consider a level \(0<t<H\) with a single right outer boundary:
\[
 \{f>t\}=(a,\infty).
\]
The true past source is entirely to the left of \(a\). The main matching result below assumes \(f\) is nondecreasing on both sides near \(a\). The preliminary and phase statements do not require this monotonicity.

## 1. Constraints on the outer side and exclusion of finite plateaus

The time primitive \(W_t=\int_t^H Q_s ds\) has a continuous density, given by
\[
 W_t(y)=\int_0^1p_s(y)ds-\int_0^t p_s(y)ds
        +\int_{f(x)\le t} K_{t-f(x)}(y-x)g(x)dx,
 \qquad K_b(z)=\int_0^b p_s(z)ds.                            \tag{1}
\]
Here \(\int_0^H Q_sds=\int_0^1P_s\delta_0 ds\) follows from the full identity as in round 2. For \(s\ge t\), every future input lies to the right of \(a\), so \(W_t=0\) on \((-\infty,a)\) by round 58. In a nonconstant exact graph this forces \(t>1\). Indeed \(t<1\) makes the first difference in (1) positive everywhere; for \(t=1\), some positive-measure set \(f<1\) supplies a strictly positive integrated kernel. Also \(a>0\): the equality of the past density to \(p_t\) on the left gives
\[
 \Phi(a)-\Phi(a/\sqrt t)
 =\int_{x<a,\ f(x)<t} \bar\Phi((a-x)/\sqrt{t-f(x)})g(x)dx>0 .    \tag{2}
\]
The positive sign follows because the past has positive mass at variances \(t-f>0\), as \(E_\gamma f=1<t\). The equality includes zero-age past sources deterministically on the left.

There is no finite plateau \([l,a]\) of positive length at this level if \(f\ge t\) on \((l,\infty)\). Define \(Q_{t-}\) using \(f<t\) in the past integral. All its past means are at most \(l\); hence it is a real analytic density on \((l,\infty)\) by off-diagonal Gaussian estimates, even when past ages tend to zero. On \((l,a)\), \(Q_t=0\) and \(Q_{t-}=Q_t+g1_{\{f=t\}}\), so \(Q_{t-}=g\) there and therefore everywhere on \((l,\infty)\). This gives \(Q_t=g1_{(a,\infty)}\,dx\). But its Laplace transform for \(s>0\) is strictly larger than the exact future formula
\[
 \int_{a}^{\infty} e^{-sx-s^2(f(x)-t)/2}g(x)dx ,
\]
since \(f(x)>t\) on the ray. This is a contradiction.

There is also a quantitative outward thickness condition. For every sufficiently small fixed \(r>0\), let
\[
 R_r(h)=|\{u\in(0,r):0<f(a+u)-t\le h\}|.
\]
Then
\[
 \liminf_{h\downarrow0} h\log(1/R_r(h))=0.                     \tag{3}
\]
The analytic density \(Q_t\) on \((a,\infty)\) is not identically zero, since its mass is \(\gamma(f>t)>0\). Hence its support starts precisely at \(a\), by the identity theorem on that connected ray. If the liminf in (3) were positive, \(R_r(h)\le C e^{-c/h}\) for small \(h\). The shifted future Laplace transform,
\[
 L_t(s)=e^{as}\int e^{-sy}Q_t(dy)
       =\int_0^\infty e^{-su-s^2(f(a+u)-t)/2}g(a+u)du,           \tag{4}
\]
would then satisfy \(L_t(s)=O(e^{-\delta s})\) for some \(\delta>0\): outside \(r\) use \(e^{-sr}\); inside \(r\), the distribution-function integral with \(k=s^2/2\) is bounded by \(C k\int_0^{h_0}e^{-kh-c/h}dh+O(e^{-kh_0})\). The inequality \(kh+c/h\ge2\sqrt{kc}\) proves the exponential bound. Applying the round-58 support argument to the shifted transform \(e^{(a+\delta')z}\int e^{-zy}Q_t(dy)\), \(0<\delta'<\delta\), forces a strict support gap to the right of \(a\), impossible. Its variation bound and quadrant argument are unchanged; the nonnegative real-axis upper bound may be any constant. In particular, a nondecreasing outward contact cannot have \(\liminf_{u\downarrow0}(f(a+u)-t)\log(1/u)>0\).

## 2. Phase and jump law on the past side

Let \(\nu^-_t=\int_{f\le t}P_{t-f(x)}\delta_xg(x)dx\) be the nonnegative past output. Its density equals \(p_t\) on \((-\infty,a)\). On \((a,\infty)\), it has a smooth density \(U(y)\), which is decreasing because every positive-variance mean is at most \(a\). Its right trace \(U_0\) is the integral of these positive-variance kernels at \(a\), by monotone convergence. Fatou from the left equality gives \(U_0\le p_t(a)\). Thus
\[
 J=p_t(a)-U_0\ge0
\]
is the right trace of the remaining density \(Q_t\). At the boundary there is no atom since the source is absolutely continuous and every zero-age mean is its original position. In particular the past density is bounded on the full line.

Put \(b_-(u)=t-f(a-u)\ge0\), \(u>0\), and \(q_-(u)=b_-(u)/u^2\). Fatou also gives finite logarithmic measure \(du/u\) near zero to each intermediate band \(\eta\le q_-\le C\), because the positive-variance kernel at \(a\) there is bounded below by \(c_{\eta,C}/u\). The flat-phase occupation converges weakly on the past half-line:
\[
 1_{\{q_-(ru)\le1\}}\ \rightharpoonup^\ast\ \theta,
            \qquad \theta=J/g(a),\quad u>0.                    \tag{5}
\]

For the proof, take any weak-star subsequential limit \(\chi\) of these bounded indicators. At scale \(r\) the past output density converges on bounded intervals in \(L^1\) to \(U_0+J1_{\{z<0\}}\), using its two traces. For a smooth compactly supported test in the output coordinate \(z\), split past source means at \(x=a-rL\). On the finite near rescaled window \(u\in(0,L)\), the local output pairing is with \(P_{b_-(ru)/r^2}\psi(-u)\), times the source density tending uniformly to \(g(a)\). On a fixed annulus away from \(u=0\), use \(q_-\) small to approximate \(\psi(-u)\), \(q_-\) large to get zero, and vanishing rescaled occupation of each intermediate band. The small interval around zero is controlled by the sup norm of the test. Thus the limiting near output is \(g(a)\chi(-z)1_{\{z<0\}}\,dz\).

For the far part, the moving-center Gaussian derivative estimate used for whole-line tangents gives an \(O(1/L)\) Lipschitz bound in the bounded output window: shift \(z\) toward the far mean by \(L/4\), use \((d/s)e^{-(L/4)d/(2s)}\le C/L\), and dominate by the already bounded full **past** output density at the shifted point. Far zero-age sources have no mass in the central windows. Pair with zero-integral tests to erase the far constant, and send \(L\to\infty\). The distribution
\((g(a)\chi(-z)-J)1_{\{z<0\}}\) annihilates these tests, so is a constant distribution and vanishes by its support. This proves (5) and its uniqueness. In particular \(0\le J\le g(a)\).

## 3. Matching of arbitrary monotone contacts

If \(b_-\) is nondecreasing near zero, then (5) forces \(\theta\in\{0,1\}\). Indeed, a sequence \(q_-(r)\le1/4\) forces the whole interval \([r/2,r]\) flat; it implies \(\theta=1\). A sequence \(q_-(r)\ge5\) forces \([r,2r]\) nonflat; it implies \(\theta=0\). If neither occurs, eventually \(q_-\in[1/4,5]\), contrary to finite logarithmic measure. More generally if \(\theta=1\), a sequence \(q_-(r)\ge\eta>0\) would on \([r,2r]\) either lie in an intermediate band or give nonflat occupation, contradicting (5). If \(\theta=0\), any sequence \(q_-(r)\le C\) likewise contradicts the absence of flat and intermediate occupation on \([r/2,r]\). Hence the pointwise contact is either \(q_-\to0\) or \(q_-\to\infty\), respectively.

By the right trace and finite signed measure,
\(sL_t(s)\to J\) as \(s\to\infty\), an ordinary Abelian limit. From (4),
\[
 sL_t(s)=\int_0^\infty e^{-\omega}g(a+\omega/s)
               e^{-s^2 b_+(\omega/s)/2}d\omega,\qquad
 b_+(u)=f(a+u)-t .
\]
If \(b_+\) is nondecreasing as well, \(J=g(a)\) forces \(b_+(u)/u^2\to0\): a sequence \(b_+(r)\ge\eta r^2\) supplies a fixed loss on \(\omega\in[1,2]\) at \(s=1/r\). If \(J=0\), a sequence \(b_+(r)\le C r^2\) gives a strictly positive contribution on \([1/2,1]\), so \(b_+(u)/u^2\to\infty\). The contact is therefore flat on both sides or steep on both sides. No regularity of either contact ratio is needed beforehand.

By (1), \(W_t\) is \(C^1\): \(K_b\) is uniformly Lipschitz with bounded derivatives away from the center, and the exceptional center is one input point of measure zero. On the right,
\[
 W_t''=2(g-Q_t).
\]
This is obtained by integrating the heat equation in time from \(t\) to \(H\), with \(Q_H=0\) and future sources restricted to \(\{f>t\}\). Together with the zero left support of \(W_t\) it gives
\[
 W_t(a+y)=(g(a)-J)y^2+o(y^2),\qquad y\downarrow0 .
\]
A flat level also requires \(p_t(a)>g(a)\), since \(p_t(a)-g(a)=U_0>0\). Thus \(a^2>t\log(t)/(t-1)\) at every flat right level. On an interval of strictly monotone outer levels, \(a(t)\) is continuous. Fatou makes the off-diagonal trace \(U_0(t)\) lower semicontinuous, so \(\theta(t)\in\{0,1\}\) is upper semicontinuous. The flat-level set is closed relative to this interval and has zero Lebesgue measure in time: at a flat crossing the inverse boundary \(a(t)\) has infinite difference quotients, whereas a monotone inverse is differentiable with finite derivative almost everywhere. The steep levels are consequently open and dense.

## 4. A further exclusion of comparable past powers between one and two

The steep class can still be restricted if
\[
 c_1u^p\le b_-(u)\le c_2u^p\quad(0<u<r_0),\qquad 1<p<2 .
\]
This entails \(J=0\) by (5), so \(U_0=p_t(a)\). Compare the value \(p_t(a-h)\) to \(U_0\) for \(h>0\). The exterior source part at distance at least \(r_0\) has difference \(O(h)\); the local part is
\[
 \Delta=\int_0^{r_0}g(a-u)p_{b_-(u)}(u)
             [1-e^{uh/b_-(u)-h^2/(2b_-(u))}]du .
\]
The part with \(1-e^{-h^2/(2b_-(u))}\) in brackets is positive and bounded below by \(c h^{2/p-1}\) by restricting to \(h^{2/p}\le u\le2h^{2/p}\). The cross term has magnitude bounded by
\[
 C h\int_0^{r_0}u^{1-3p/2}
                   e^{-(u-h)^2/(2b_-(u))}du .
\]
Split at \(h/2\). On the lower interval use \(e^{-c h^2/u^p}\); on the other interval discard the exponential. For \(p<4/3\) the bound is \(O(h)\), for \(p=4/3\) it is \(O(h\log(1/h))\), and for \(p>4/3\) it is \(O(h^{4/p-2}+h^{3-3p/2})\) using \(u=h^{2/p}z\) in the lower part. In all cases this is \(o(h^{2/p-1})\). Thus \(U_0-p_t(a-h)\) is at least \(c h^{2/p-1}+o(h^{2/p-1})\), contrary to the analytic left value \(O(h)\).

An ordinary linear crossing is locally compatible with the half-line equality: for constant source density and \(b_-(u)=cu\), \(u>0\),
\[
 \int_0^\infty p_{cu}(y+u)\,du=e^{-(y+|y|)/c}.
\]
The profile is 1 on the left and smaller on the right. The argument therefore gives local exclusions and restrictions for arbitrary monotone contacts while leaving general exact outer branches unresolved.

## 5. Localization to a single exterior interface and local verification

The trace and phase laws, matching conclusion, and comparable-power exclusion require only
\(\{f>t\}\subset(a,\infty)\) and a nonempty right neighborhood
\((a,a+\delta)\subset\{f>t\}\). Past sources at \(x\ge a+\delta\) give a smooth bounded term through the local interface; the other sources and the Fatou argument are unchanged after separating this term. The moving-center derivative estimate requires boundedness of the past density only near the interface. The left support and the proof that \(t>1\) still hold. Formula (2) must now subtract any leakage from right past sources, so it is used above only for the full future ray.

The plateau exclusion also localizes. If \(f=t\) on \((l,a)\) and \((a,a+\delta)\) belongs to the future, the preactivation density is analytic on \((l,a+\delta)\), because all strictly past means lie outside it. It equals \(g\) on the plateau, implying that \(Q_t=g\) on \((a,a+\delta)\). The difference of this local \(Q_t\) and the future formula in (4) gives
\[
 \int_0^\delta e^{-su}\bigl(1-e^{-s^2(f(a+u)-t)/2}\bigr)g(a+u)du
      =O(e^{-s\delta}).
\]
All terms beyond \(\delta\) have this bound by finite total variation. Any closed interval \(J\subset(0,\delta)\) has a strictly positive lower bound on \(f(a+u)-t\), and its left-side contribution is at least \(c_J e^{-s\sup J}\) for large \(s\), a contradiction. In the local setting, the thickness assertion (3) additionally requires the right density to be nonzero on the analytic neighbor, so that its support starts at \(a\).

The companion calculation uses Python 3.13.18, SymPy 1.14.0 and mpmath 1.3.0. SymPy verifies the positive cross-term exponent differences and the one-dimensional heat resolvent. At 70 decimal digits the straight crossing with \(c=1.3\) matches its exact profile at five points; the power model \(p=3/2\) has normalized cusp discrepancy \(1.6043550948482\) at \(h=2^{-36}\), approaching \(1.6047611497507\). The outward logarithmic-barrier model gives \(-s^{-1}\log L(s)\approx1.41215716\) at \(s=2048\), approaching the analytic saddle rate \(\sqrt2\). These are local model checks, and support and contact restrictions in the proof are analytic.
