# Gaussian smoothing of a fat Cantor set

Start from \([-1/2,1/2]\). At stage \(n\ge1\), remove the middle relative fraction \(c2^{-n}\) of each surviving interval, with \(0<c<2\). Let \(K_c\) be the limiting closed set and \(H_v=p_v*1_{K_c}\), where \(p_v\) is the Gaussian density of variance \(v\). This report uses the pure smoothing, with no Gaussian background factor. The reference target is
\[
\theta_{\rm ref}=\frac{p_{7/8}(2)}{p_1(2)}=\sqrt{8/7}\,e^{-2/7}=0.803363018465151\ldots .
\]
For each \(0<\theta<1\), there is a unique variance where \(\max_x H_v(x)=\theta\). The maximum decreases strictly with \(v\) by the Gaussian semigroup, tends to 1 at zero by the density theorem, and tends to 0 at infinity. Symmetry and \(H_v'(x)<0\) for \(x>1/2\) place every maximizing location in the searched half interval or its reflection.

## Center and split certificates

Use positive-half moments
\[
M_n(v)=\int_{K_c\cap[0,1/2]} y^{2n}e^{-y^2/(2v)}\,dy.
\]
The center curvature has the sign of \(M_1-vM_0\). For \(x>0\), after a positive factor is removed, \(H_v'(x)\) has the series
\[
\sum_{n\ge0}\left[M_{n+1}-(2n+1)vM_n\right]\frac{(x/v)^{2n+1}}{(2n+1)!}.
\]
Once \((2n+1)v\ge1/4\), that coefficient is nonpositive by support. Thus finitely many signed moment tests prove global strict unimodality. Convolution with another Gaussian preserves this property. More explicitly, if \(H_v'<0\) on the positive half-line, then for each \(h>0\)
\[
H_{v+h}'(x)=\int_0^\infty H_v'(y)[p_h(x-y)-p_h(x+y)]\,dy<0
\quad(x>0),
\]
and \(H_{v+h}''(0)<0\) follows by differentiating this strictly signed formula at zero.

The interval verifier encloses each \(M_n\) using a finite stage whose components have length \(\ell_N\), weighted by the retained fraction \(P_{>N}\). In each component the true indicator and weighted stage indicator have the same zeroth and first moments, so for a test function \(g\) the absolute error is bounded by \((1-P_{>N})\,|K_{c,N}|\ell_N^2\|g''\|_\infty/4\), with half the stage length for moments on the positive half. The bound uses Taylor's theorem and the signed total variation at most \(2(1-P_{>N})\ell_N\) per interval. The finite-stage quadrature uses the midpoint value plus its second derivative correction, with an integrated fourth-derivative remainder. The tail product is enclosed by a finite product and the sum bound on remaining fractions.

For each listed \(c\), a positive lower curvature sign and an upper strict-unimodality certificate enclose its first split as \(v\) decreases from large variance. The interval arithmetic also verifies strictly negative higher moment coefficients throughout each bracket, using the support bound for the remaining orders. Therefore the zero inside each bracket has a unique center maximum and the same local square-root split proved below.

| \(c\) | Certified split variance bracket | Approximate height at split |
| ---: | ---: | ---: |
| \(1/256\) | [0.019245260, 0.019245461] | 0.992553765 |
| \(1/64\) | [0.024510711, 0.024510912] | 0.972180928 |
| 0.125 | [0.042580450, 0.042580651] | 0.813554761 |
| 0.25 | [0.056493368, 0.056493569] | 0.668048970 |
| 0.5 | [0.081763036, 0.081763237] | 0.448549861 |
| 1 | [0.133279387, 0.133279588] | 0.186407408 |
| 1.375 | [0.174896842, 0.174897043] | 0.081605221 |

At \(c=1/256\), the target \(\theta_{\rm ref}\) is reached at a unique centered maximum, with certified variance in \([0.147530539,0.147530740]\). In contrast, at \(c=0.25\), the first gap alone gives \(H_v(0)<0.757396527\) for all \(v\), so the maximum reaching \(\theta_{\rm ref}\) is off-center. The latter computes the largest Gaussian mass of \([c/4,1/2]\cup[-1/2,-c/4]\) over \(v\). Its closed-form maximizing variance is \(((1/2)^2-(c/4)^2)/(2\log((1/2)/(c/4)))\).

## Small-removal asymptotics

For small \(c\), the first central gap has half-width \(a=c/4\). Denote the union of all later gaps by \(G\). Uniformly under small \(c\),
\[
|G\cap[-w,w]|\le Cc w^2 .
\]
A stage-\(n\) descendant gap can meet this window only when \(2^{-n}\lesssim w\), since its distance from zero is bounded below by a fixed multiple of the parent component length. Its length is \(O(c4^{-n})\) and at most \(C(1+2^n w)\) such gaps can meet the window. Sum over those \(n\). Consequently \(\int_G (1+y^2/v)e^{-y^2/(2v)}dy=O(cv)\) and \(\int_G y^2e^{-y^2/(2v)}dy=O(cv^2)\).

The curvature at the center therefore has the form
\[
H_v''(0)=\frac{1}{\sqrt{2\pi}v^{3/2}}
\left[-e^{-1/(8v)}+\frac{c}{2}\,e^{-c^2/(32v)}+O(cv)\right].
\]
A zero lies near \(v_0=[8\log(2/c)]^{-1}\) by evaluating at \(v_0/2\) and \(2v_0\). At any zero in this range, \(e^{-1/(8v)}\sim c/2\). For \(n\ge1\) and \((2n+1)v<1/4\), integration by parts over \([0,1/2]\) and its gaps gives the normalized outer boundary contribution \(-e^{-1/(8v)}\) to \([M_{n+1}-(2n+1)vM_n]/[v(1/2)^{2n+1}]\). The positive part from the central gap is \(O(c^3/v)\), and the later gaps give \(O(cv)\). To check these last bounds use \((2n+1)=O(1/v)\) and \((2y)^{2n}\le(2y)^2\) in the gap integral. All higher moment coefficients are therefore negative at a curvature zero in this range. The support condition handles the omitted high orders.

It follows that a curvature zero here has a unique center maximum and that larger variances have strictly negative curvature. At smaller variances \(v<a^2\) the curvature is positive by direct integration. For \(a^2\le v\le v_0/2\), the displayed curvature formula gives positivity since the central contribution is at least a fixed multiple of \(c\). Hence there is a unique center-to-split variance \(v_*(c)\) for sufficiently small \(c\), and
\[
v_*(c)\sim \frac{1}{8\log(2/c)},\qquad
1-H_{v_*(c)}(0)\sim c\sqrt{\frac{\log(2/c)}{\pi}}.
\]
The last claim follows by integrating the first gap and bounding the outer and later-gap mass deficits by \(O(c\sqrt{v_*})\).

At the zero, the negative coefficient for \(n=1\) implies \(H_{v_*}^{(4)}(0)<0\). The heat equation \(\partial_v H_v=\frac12H_v''\) gives the crossing derivative of the curvature. Apply the implicit function theorem to \(H_v'(x)/x\) as a function of \(x^2,v\), using uniqueness of the global center at \(v_*\) to localize nearby maxima. The two maximizers immediately below the split obey
\[
x_\pm(v)=\pm\sqrt{3(v_*-v)}+O((v_*-v)^{3/2}).
\]

## MacBook Metal and independent replay

A fixed scan covered 137 removals in \([1/512,1.375]\), 225 logarithmic variances in \([2^{-32},1/2]\) and 6,144 spatial sites per pair. Spatial candidates combine a uniform mesh and points generated by Cantor addresses adapted to the Gaussian width. The kernel integrates over weighted stage components via three-point Gauss–Legendre, clips by the eight-standard-deviation window and saves height, center curvature and sampled maximizer. The table in `thresholds.json` brackets the sampled level for each \(c\) and all targets 0.70, 0.71, …, 0.99 plus the reference. There were 4,247 sampled level brackets. At 1,921 lower-bracket grid variances, the sampled peak exceeds the center by more than \(10^{-5}\).

Independent CPU comparisons use normal CDF integrals at stage depths 14 and 16 and local refinement of selected sampled maxima. The depth check differs by at most \(9.164\times10^{-11}\) on the selected inputs. The largest discrepancy between the sampled GPU value and CPU at the same point was \(6.886\times10^{-5}\), occurring at a very small variance with \(c=1.375\). The full table is a grid scan. The exact shape statements use only the interval certificates and analytic inequalities above.

The production run evaluated 189,388,800 sites on Apple M5 Pro from 10:08:02.064 to 10:08:05.217 UTC on 2026-11-18. Metal command time was 2.629848875 s over 3.153625125 s elapsed. Device-wide utilization samples were 96%, 99% or 100%. The two CPU replays agree across local and remote hosts. Both smoke and main runs exited successfully, and their two PIDs were absent after collection.

### Reproduce

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc157-temp"
TMPDIR="$cache/cc157-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o cantor_scan
sh run.sh main 129 225 6144 64
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json thresholds.remote.json
```
