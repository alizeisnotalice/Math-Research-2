# One-sided density jets at finite-order contacts

Write \(g=p_1\), with center \(a\). On the side \(\sigma=\pm1\) let
\[
v(a+\sigma z)=v_\sigma(z)=c_\sigma z^{n_\sigma}(1+b_\sigma z+d_\sigma z^2),\quad z>0,
\]
where \(n_\sigma>2\) is integer and \(c_\sigma>0\), on a sufficiently small fixed core. The positive shape factor can also be replaced by a one-sided real-analytic germ with constant term one. Let \(F(a+y)\) be the local output density for the random map \(X+\sqrt{v(X)}Z\), with independent standard Gaussians \(X,Z\). The remote source contribution, if added, is analytic in \(y\) near zero because its source distance from the observation is bounded away from zero.

For each finite Taylor truncation, the core density as \(z=|y|\downarrow0\) has common leakage coefficients plus the side-dependent diagonal coefficients
\[
\sum_{\ell\ge0}\frac{1}{2^\ell\ell!}
 \frac{d^{\,2\ell}}{dz^{\,2\ell}}
 \left[g(a+\sigma z)v_\sigma(z)^\ell\right].
\]
The expression gives finite asymptotic coefficients, with no convergence claim for its complete series. The first term is \(g(a+\sigma z)\); the first correction is \((g v_\sigma)''/2\) in the positive side coordinate.

## Common coefficients away from the diagonal

For source coordinate \(|u|\ge2z\), the difference \(y-u\) has magnitude at least \(|u|/2\). Every fixed output derivative of the Gaussian kernel has a uniform bound of the form \(C_Nv(a+u)^{-(N+1)/2}e^{-u^2/(16v(a+u))}\) after absorbing the Hermite polynomial into a weakened Gaussian. It is integrable at the contact because \(v(a+u)=O(|u|^{n_\sigma})\) with \(n_\sigma>2\). Taylor expand in \(y\) with this bound. Each omitted coefficient integral over \(|u|<2z\) tends to zero faster than every power of \(z\). For sources with \(|u|<2z\) and outside a fixed relative neighborhood of the matching diagonal \(u=\sigma z\), displacement is at least a fixed multiple of \(z\). Together with the local variance upper bound on each side, this also gives an error smaller than every power of \(z\). The resulting leakage jet is the same from the two sides.

## Conditional inverse for the diagonal

Use the positive side coordinate. Put \(\tau(z)=\sqrt{v_\sigma(z)}\), \(\lambda(z)=\tau(z)/z\). Conditional on fixed noise \(\zeta\), the matching branch is the inverse of \(s\mapsto s+\zeta\tau(s)\), absorbing the unimportant sign into the symmetric normal noise. For \(|\zeta|\le\lambda(z)^{-1/2}\), there is a unique local inverse inside a fixed relative neighborhood of \(z\). Indeed \(\tau(z+zw)/\tau(z)\) is uniformly analytic for small complex \(w\), and the scaled perturbation parameter is \(\zeta\lambda(z)\). A uniform local implicit-function radius after this rescaling allows finite Taylor remainder bounds in the scaled parameter. If \(s(z,\zeta)\) denotes the inverse, the conditional density there is \(g(a+\sigma s(z,\zeta))/(1+\zeta\tau'(s(z,\zeta)))\).

Lagrange inversion, applied to an antiderivative of \(g(a+\sigma s)\) and differentiated in \(z\) with \(\zeta\) fixed, gives the coefficient of \((-\zeta)^m\) in the conditional density as \(m!^{-1}(d/dz)^m[g(a+\sigma z)\tau(z)^m]\). Uniform finite Taylor remainders in \(\zeta\lambda(z)\) are integrable against the standard Gaussian. The rest of the noise range contributes a term smaller than every power of \(z\), as can also be seen directly from the diagonal kernel by \(u=\sigma(z+\tau(z)t)\) and bounded nearby variance ratios. Using the even moments of \(\zeta\) yields the displayed diagonal coefficients.

In particular let \(n=\min(n_+,n_-)\) and set \(c_\sigma^*=c_\sigma\) on sides whose order is \(n\), zero otherwise. The jump in the order-\((n-2)\) Peano derivative, right minus left in the common signed coordinate \(y\), is
\[
\frac{g(a)n!}{2}\big[c_+^*-(-1)^n c_-^*\big].
\]
The lower-order Peano derivatives agree. A nonzero jump is incompatible with an analytic target. If the two sides are smooth up to their boundary, the same result is a derivative-limit jump. To see this, one may differentiate the conditional inverse on the relative neighborhood with a smooth noise cutoff and increase the fixed Taylor remainder order. Derivatives of the cutoff and the omitted large-noise contribution are bounded by inverse powers of \(z\) times the same superpolynomial Gaussian tail.

For the next term on a side, divide the diagonal correction by \(g(a)z^{n_\sigma-2}\). Its constant is \(c_\sigma n_\sigma(n_\sigma-1)/2\), and its \(z\)-coefficient is
\[
\tfrac12 c_\sigma n_\sigma(n_\sigma+1)(b_\sigma-\sigma a)
        +45c_\sigma^2\,\mathbf1_{n_\sigma=3}.
\]
For cubics the additional term comes from \((g v_\sigma^2)^{(4)}/8\) at the same order. This term is absent at the displayed order for \(n_\sigma\ge4\).

## Computation

The script integrates the diagonal change of variables on \(|t|\le20\) and compares the correction to \(g(a+\sigma z)\). Its two numerical forms are 87-digit direct log-ratio integration and binary64 stable even pairing, which avoids subtracting two near-equal densities or cancelling odd terms in binary64. The common leakage cancellation is established by the argument above. Finite-scale observations calibrate the constants and the sign conversion between positive side coordinates and the signed output coordinate.
