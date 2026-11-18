# Round 79 MacBook scan of global call inequalities

Let (Z=X+\sqrt{v(X)}Y), (M=1+\mathbb E v(X)), and
let (G_M) have distribution (N(0,M)). The strike discrepancy is
[
 \Delta(k)=\mathbb E[(Z-k)_+]-\mathbb E[(G_M-k)_+].
]
Both convex orders hold exactly when the call functions agree.
This experiment uses continuous bounded periodic (v), with
reference variance (M) and exact fourth-moment matching.

## Matched profile and stable integration

For (r_omega=e^{-\omega^2/2}), set
[
 g(x)=\cos(kx)-r_k+\eta(\cos(lx)-r_l),\qquad
 v(x)=m+\lambda g(x),\qquad
 \lambda=2(k^2r_k+\eta l^2r_l)/\operatorname{Var}(g(X)).
]
Then (\mathbb E g(X)=0) and
(\operatorname{Cov}(X^2,g(X))=-(k^2r_k+\eta l^2r_l)).
Thus (3\operatorname{Var}(v(X))+6\operatorname{Cov}(X^2,v(X))=0),
which is exactly the fourth-moment condition. The law is symmetric.
The grid screens out standard deviations of (v(X)) below 0.004,
requires the global bound (v\le m+1) by a cosine amplitude
bound, and keeps (v) at least 0.35.

The tested grid has 57 frequencies (k\in[2.5,7]), 51 frequencies
(l\in[k+0.15,9]), 81 mixing coefficients (\eta\in[-2,2]),
and five logarithmic positive noise margins in ([0.35,6]).
Each accepted profile is evaluated at 81 standardized strikes
(y=k_{\rm strike}/\sqrt M\in[0,7]).
The reported call discrepancy is (\Delta(y\sqrt M)/\sqrt M).
To avoid subtracting large calls, the shader uses the identity
[
 \partial_v\mathbb E[(x+\sqrt vY-k_{\rm strike})_+]
 =\varphi((k_{\rm strike}-x)/\sqrt v)/(2\sqrt v)
]
and integrates between (m) and (v(x)) by four-point
Gauss–Legendre quadrature. A 512-panel midpoint rule on ([0,9])
and symmetry give the Gaussian outer integral. CPU validation
uses adaptive outer integration with a 16-node variance rule,
an independent explicit conditional call evaluation, and
60-decimal working-precision checks for selected small residuals.

## Exact sign and tail conclusions

For a symmetric law and (q>1),
[
 \mathbb E|Z|^q=2q(q-1)\int_0^\infty k^{q-2}
                   \mathbb E[(Z-k)_+]\,dk .
]
If the even moments through (2p) match, (Delta) is orthogonal
to (1,k^2,\dots,k^{2p-2}) on the half-line. A continuous nonzero
(Delta) with fewer than (p) sign changes would have a fixed
sign after multiplication by a polynomial in (k^2) of degree
at most (p-1), using one factor at each sign transition.
This contradicts the moment equalities. In particular the
equal-variance and exact-fourth-moment profile requires at least
two sign transitions.

A bounded nonconstant continuous periodic (v) has
[
 \lim_{k\to\infty}k^{-2}\log\mathbb E[(Z-k)_+]
 =-\frac1{2(1+\max v)}.
]
For the upper bound use the Gaussian mgf with (v) replaced
by its maximum. For the lower bound pick (v_0<\max v) and
a periodically recurrent interval of positive length where
(v\ge v_0). For large (k), a translate lies at bounded
distance from (k/(1+v_0)). Restrict (X) to that interval
and (Y) to an interval permitting (Z\ge k+1); its log
probability is at least (-k^2/[2(1+v_0)]+O(k)).
Taking (v_0\uparrow\max v) gives the stated rate.
Since the maximum exceeds (m), the discrepancy is eventually
positive. The accepted cosine frequencies are rational, so
their profiles are periodic and satisfy this conclusion.

## Large common noise and the first unmatched Hermite term

Fix the even perturbation (v-m) and add constant variance.
At an admissible reference (m_0), let (W_0) denote the
mixture output, (M_0=1+m_0), and
(\delta_6=\mathbb E W_0^6-15M_0^3).
For (M=1+m\to\infty), the standardized discrepancy has
the uniform expansion
[
 \frac{\Delta_m(y\sqrt M)}{\sqrt M}
 =\frac{\delta_6}{720M^3}H_4(y)\varphi(y)+O(M^{-4}).
]
Indeed the signed base law has moments through five zero.
Its characteristic function is
(\delta_6(i\xi)^6/720+O(|\xi|^8)). Adding the normal variance
and standardizing multiplies its rescaled transform by
(\exp[-(1-M_0/M)\xi^2/2]). Dividing the resulting transform
by (-\xi^2) gives the call transform, since the second strike
derivative is the density difference. The Fourier remainder is
bounded in (L^1) by (O(M^{-4})); the leading inverse transform
is (H_4\varphi). For a later first unmatched even moment
the same argument gives the corresponding Hermite shape.

The CPU calibration with (k=3,\eta=0) has raw sixth excess
5.40727393917. At (m=256), the call residual divided by
(\delta_6/(720M^3)) is 1.20396535 at (y=0)
and -0.74053979 at (y=1.4), versus Hermite limits
1.19682684 and -0.73641957.

## Numerical outcome and receipts

The main counted **35,492,580 finite strike tests** over
438,180 accepted parameter combinations. Every profile had
at least two numerical mesh sign transitions above the
relative/absolute threshold, and 351,250 had exactly two.
The leading seven-transition representative was independently
confirmed by both CPUs, with zeros at approximately
0.34914, 1.09977, 1.94025, 2.86388, 3.65255, 4.81246,
and 5.29061 standardized strike units. These are detected
transitions of a representative, and the mesh does not give
a global count of all zeros.

For the least sensitive selected profile after dividing by the
noise standard deviation, the strongest positive and negative
discrepancies refined to (1.99135\times10^{-9}) and
(-1.64599\times10^{-9}) at standardized strikes 0.88669 and
1.93052. At the strongest positive profile the discrepancy was
about 0.007506 at zero, while the strongest negative profile
refined to about -0.004385 at 1.06705. Detailed coefficients
and 60-decimal checks are in the CPU output.

On the Apple M5 Pro, production PID 49593 ran 04:39:15–04:39:18
UTC on 2026-11-18 and exited zero. Its native command buffers
totaled **2.319082 seconds**. All five whole-device GPU samples
were 100%; a separate earlier device probe also read 100%.
The post-run process query found neither production nor pilot
PID. This was a finite mesh run.

From the authenticated MacBook directory
`~/CodexRemoteMath/cc79-metal-nov18`:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc79-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o call_scan
sh run.sh smoke 5 5 9 3 21 128 1024
sh run.sh sync 57 51 81 5 81 512 8192
"$HOME/.venvs/math/bin/python" check_cpu.py sync.json sync.cpu.json
"$HOME/.venvs/math/bin/python" smoothing_check.py smoothing.remote.json
```

The MacBook ran macOS 27.3 build 26C5035r, Swift 6.4.2,
Python 3.12.13, SciPy 1.18.1, and mpmath 1.3.0.
The grid is deterministic. Cross-device hashes of selected
sources and outputs are in `remote.sha256`.
