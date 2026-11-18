# All-line microscopic tangent and periodic rigidity at the noise scale

## 1. Microscopic output is uniform

Let the source have bounded continuous density \(g\), and let the output of \(X+\sqrt{v(X)}G\), \(G\) an independent standard normal, have bounded continuous density \(p\). Assume \(v\ge0\). Let \(x_n\to x_0\), with \(g(x_0)>0\), and \(\ell_n\to0^+\). Suppose
\[
 a_n(u)=\frac{v(x_n+\ell_n u)}{\ell_n^2}
     \longrightarrow a(u)
 \quad\text{in measure on every bounded interval},
\]
where \(a\) is nonnegative and bounded on the full line. No uniform bound on \(a_n\) or power envelope on \(v\) is required. Then the unit Lebesgue source transported by \(u+\sqrt{a(u)}G\) has unit Lebesgue output. Equivalently,
\[
 U_a(z)=\int_{\mathbb R}\int_0^{a(u)} p_s(z-u)\,ds\,du
 \quad\text{is constant},                                   \tag{1}
\]
where \(p_s\) is the Gaussian density of variance \(s\).

Normalize the output measure in \(z=(Y-x_n)/\ell_n\) by dividing it by \(\ell_n\). It has density \(p(x_n+\ell_n z)\), and its source in the \(u\)-coordinate has density \(g(x_n+\ell_n u)\). Fix an output test supported in \([-R,R]\), and split sources at \(|u|=L>2R\). Denote the far density in \(|z|<L\) by \(F_{n,L}(z)\). Far sources with zero variance are absent there. The positive-variance far contribution is smooth there by the uniform Gaussian off-diagonal derivative bounds and the finite mass of the normalized source at every fixed \(n\). By positivity and equality of the full output, \(F_{n,L}(z)\le p(x_n+\ell_n z)\) pointwise in that interval.

Far density becomes locally constant in a manner uniform in both \(n\) and the variance. For \(h=L/4\), \(d=|u-z|\), and \(\xi=\operatorname{sign}(u-z)\), one has
\[
 \frac{p_s(z-u)}{p_s(z+\xi h-u)}
  =\exp[-h(d-h/2)/s]\le e^{-hd/(2s)},\qquad d\ge h,
\]
hence
\[
 |\partial_z p_s(z-u)|
  \le (d/s)e^{-hd/(2s)}p_s(z+\xi h-u)
  \le \frac{2}{e h}\,p_s(z+\xi h-u).
\]
Split the far source into its two signs and integrate. Since \(z\pm h\in(-L,L)\), this proves
\[
 \sup_{|z|\le R}|F_{n,L}'(z)|\le C\|p\|_\infty/L .             \tag{2}
\]
The bound does not require \(a_n\) small or large in the exterior. The large far *constant* is harmless.

Take a smooth compactly supported test \(\psi\) with zero integral. The full output pairing tends to zero by continuity of \(p\). The far pairing is bounded in absolute value by \(C_\psi/L\) using (2), after subtracting the constant \(F_{n,L}(0)\). On the finite near source interval, the pairing is
\[
 \int_{-L}^L g(x_n+\ell_n u)(P_{a_n(u)}\psi)(u)\,du ,
\]
where \(P_0\psi=\psi\). For continuous bounded \(\psi\), \(P_s\psi(u)\) is continuous in \(s\in[0,\infty)\), including at zero, and bounded by \(\|\psi\|_\infty\). Convergence in measure plus this uniform bound proves that the near pairing tends to
\(g(x_0)\int_{-L}^L (P_{a(u)}\psi)(u)\,du\).
Now let \(L\to\infty\). The Gaussian tail bound from boundedness of \(a\) and compact support of \(\psi\) makes this integral absolutely convergent. We obtain
\[
 \int_{\mathbb R} (P_{a(u)}\psi)(u)\,du=0
                       \quad\text{if }\int\psi=0.
\]
Thus the micro output is a constant times Lebesgue measure.

The constant is one. The heat equation, tested distributionally, gives
\[
 U_a''=2(D_a-dz),\qquad
 D_a(\psi)=\int (P_{a(u)}\psi)(u)\,du.
\]
Since \(0\le U_a\le\sup a\) by comparison to the constant-variance integrated kernel, a constant \(U_a''\) has to be zero. Bounded \(U_a\) with zero distributional second derivative is itself constant. This derivation uses the integrated heat potential without needing continuity of the micro output density or a positive lower variance.

A location consequence is
\[
 p(x_0)\ge g(x_0).
\]
Indeed, repeat the split with a nonnegative test of total mass one and drop the far mass, taking \(n\to\infty\) and then \(L\to\infty\). If \(g\) is positive Gaussian and \(v\) is continuous and nontrivial, choose a fixed nondegenerate source interval \(J\) separated from \(x_0\) with \(v\ge v_0>0\). This interval is disjoint from every finite near micro window once \(n\) is large. Retaining its positive smooth output contribution gives
\[
 p(x_0)\ge g(x_0)+\int_J p_{v(x)}(x_0-x)g(x)\,dx>g(x_0).
\]
In the Cantelli setting this extends the earlier center exclusion to bounded microscopic profiles at moving locations, without requiring an attained zero.

## 2. Rigidity for bounded profiles with a quadratic upper contact at a minimum

Let \(a\) be a bounded continuous micro profile satisfying \(D_a=dz\). Suppose its global minimum \(c\) is attained, translated to zero, and \(a(u)-c\le C u^2\) at zero on both sides. Then \(a\) is constant.

Put \(b=a-c\ge0\). The output for \(a\) is Gaussian convolution of the output \(D_b\) with variance \(c\). The measure \(D_b\) is translation bounded since \(b\) is bounded: for a fixed output interval, the contribution from distant source means is controlled by a uniform Gaussian tail. Thus it is a tempered measure. If \(c>0\), Gaussian convolution is injective on these measures, by taking Fourier transforms and dividing by the nonzero Gaussian multiplier on the compact support of an arbitrary local smooth test. Hence \(D_b=dz\), also trivially when \(c=0\). For a periodic profile one may instead prove this by multiplying every nonzero Fourier identity by \(e^{k^2c/2}\).

Regularize with variance \(\varepsilon>0\). At output zero,
\[
 \int_{\mathbb R}p_{b(u)+\varepsilon}(-u)\,du=1 .
\]
Fatou shows \(\int_{\{b>0\}}p_{b(u)}(-u)\,du\le1\). If \(q(u)=b(u)/u^2\), then for every \(0<\eta<C\) the subset where \(\eta\le q\le C\) has finite logarithmic measure \(du/|u|\) near zero, because there the unregularized kernel is at least \(c_{\eta,C}/|u|\). Using \(q\le C\), this implies
\(q(\sqrt{\varepsilon}\,w)\to0\) in measure on every fixed \(w\)-interval. For example on a fixed annulus, the rescaled measure of the set \(q\ge\eta\) is at most a constant times its logarithmic measure truncated to \(|u|\le O(\sqrt{\varepsilon})\), which tends to zero. A small interval about \(w=0\) is harmless.

Consequently the contribution of \(|u|\le R\sqrt{\varepsilon}\), after scaling, tends to
\[
 \int_{-R}^R
 \lim_{\varepsilon\downarrow0}^{\text{in measure}}
 p_{1+w^2 q(\sqrt{\varepsilon}w)}(-w)\,dw
 =\int_{-R}^R p_1(w)\,dw .
\]
The kernel is uniformly bounded by \((2\pi)^{-1/2}\). If \(b\) is nonzero, choose a fixed positive-variance interval separated from zero, making a fixed positive contribution in the limit. Send first \(\varepsilon\downarrow0\), then \(R\uparrow\infty\). This contradicts total density 1.

In particular every continuous periodic nonconstant \(C^2\) profile is excluded: it attains its minimum, and Taylor's upper bound there is quadratic. The same follows if the derivative is merely locally Lipschitz near one global minimizer. The conclusion applies to any phase shift of a periodic microscopic limit.

In fact a one-sided quadratic upper bound at a minimum suffices. Apply the phase argument of round 21 with constant source weight. Fatou gives finite logarithmic measure for every intermediate band of \(q\), and the upper-bound side is flat in rescaled measure. For a weak-star subsequential limit \(\chi\) of \({\bf1}_{q(\sqrt{\varepsilon}w)\le1}\), the differentiated smoothed output identity at \(y=\sqrt{\varepsilon}z\), after multiplication by \(\sqrt{\varepsilon}\), gives \((p_1*\chi)'(z)=0\). On a fixed annulus, small \(q\) gives the standard kernel, large \(q\) gives zero, and the intermediate occupation tends to zero; the whole derivative kernel has integrable tails \(O(|w|^{-2})\). Therefore \(\chi\) is constant by Gaussian injectivity and equals 1 because it equals 1 on the controlled half. The local density comparison then gives the same contradiction. A completely arbitrary rough periodic profile is outside this quadratic upper-contact argument.

## 3. Consequence for near-macroscopic peaks

In the contact notation use \(x_n=x_0+r_n\), \(\ell_n=r_n^{q/2}\), \(q>2\). Any value convergence of \(v(x_n+\ell_n u)/r_n^q\) in measure on every fixed microscopic source interval to a bounded profile gives (1), independently of the variance in the rest of the macroscopic core. In particular, a continuous \(C^2\) periodic profile of nonzero amplitude cannot appear as that full-line limit. A leading area order or a uniform envelope on the macroscopic core is useful for determining the value of the constant in (1), but is unnecessary for rigidity.

## Local verification

The companion script ran under Python 3.13.18. SymPy 1.14.0 verified the exact shifted Gaussian ratio exponent \(-dh/s+h^2/(2s)\) and the bound \(\sup_{w\ge0}w e^{-hw/2}=2/(eh)\). At 100 decimal digits, mpmath 1.3.0 calibrated both regularized minimum mechanisms. For \(b(u)=1-\cos u\), the two-sided radius-\(1/4\) contribution divided by \(\log(1/\varepsilon)\) approaches \(e^{-1}/\sqrt\pi\). For \(b(u)=(1-\cos u)^2\), the piece \(|u|\le16\sqrt\varepsilon\) tends to the standard Gaussian mass of \([-16,16]\), while the disjoint piece \(2\le u\le3\) tends to \(0.084438690492693\ldots\). These are independent model calibrations; the proof of the all-line tangent uses the exact output and uniform bound (2).
