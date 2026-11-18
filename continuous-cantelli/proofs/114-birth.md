# Early separated wells and two rigorous birth obstructions

Use the exact normalized graph \(f\ge0\), standard Gaussian source \(g=p_1\), and signed remainder \(Q_t\). Suppose the two tails have \(\liminf f\ge1+T\) and \(f(0)=0\). For \(0<t<T\), round 87 gives signed support inside the convex hull of \(\{f\le t\}\). No sign of \(Q_t\) is assumed in a disconnected interval.

## 1. Sharper mass bound for a separated outer well

Suppose \(t_0\in(0,T)\) is a rightmost record birth: \(b=B(t_0-)<x_0=B(t_0)\), \(d=x_0-b>0\), where \(B(t)=\max\{f\le t\}\). Assume \(x_0\) is an isolated strict minimum locally and its two flanks are nondecreasing away from \(x_0\) on a fixed small neighborhood. There may be other wells globally and no regularity of the flank derivatives is assumed. If \(M(h)\) is the Gaussian mass of the small sublevel in this particular well at excess level \(h\), then
\[
 M(h)\le C h^{5/2}e^{-d^2/(2h)}\qquad(h\downarrow0).             \tag{1}
\]
This improves the \(h^2\) reflected-wall bound by using the exponentially narrow spatial localization on the newborn well.

Choose a fixed spatial interval \([x_0-r,x_0+r]\), \(r<d/4\), within the local monotone flanks, with both endpoint levels above \(t_0\). For small \(u>0\), the local sublevel at height \(t_0+u\) is an interval of mass \(M(u)\), and its width is \(\asymp M(u)\) by the positive smooth source density. The reflected round-87 bound gives
\[
 M(u)\le C u^2 e^{-(d-r)^2/(2u)};
\]
in particular both local endpoint shifts are exponentially small. Outside this neighborhood to the right, \(f>t_0\) with a uniform gap by compactness and the far-tail bound. Therefore for small \(u\) the actual rightmost endpoint is the outer endpoint of this well.

Write \(h_n=(h^{-1}+n)^{-1}\) and \(\delta_n=h_n+h_n^2\). Choose an observation \(y_n\) strictly to the right of the actual rightmost endpoint at time \(t_0+\delta_n\) with \(y_n-x_0\) exponentially small in \(1/h_n\). This is possible by the preceding endpoint bound. The strict-past heat balance at that observation is
\[
 P_{\delta_n}Q_{t_0-}(y_n)=
  \int_{t_0\le f(x)\le t_0+\delta_n}
        p_{t_0+\delta_n-f(x)}(y_n-x)\,\gamma(dx).                 \tag{2}
\]
It is a pointwise equality; age-zero sources are spatially separated there. For a newborn source with \(h_{n+1}<f(x)-t_0\le h_n\), its age is between \(h_n^2\) and \(2h_n^2\), and \(|y_n-x|=o(h_n)\) uniformly. Thus its kernel in (2) is at least \(c/h_n\). Old signed heat has its positive part bounded by \(C(b-x)\) for \(x\le b\), by the reflected Gaussian comparison used in round 87. Therefore
\[
 P_{\delta_n}Q_{t_0-}(y_n)\le C\delta_n^{3/2}
           e^{-(y_n-b)^2/(2\delta_n)}
 \le C' h_n^{3/2} e^{-d^2/(2h_n)}.
\]
Here \(y_n\ge x_0\) and \(h_n^{-1}-\delta_n^{-1}\) is bounded. Keeping only the indicated positive well sources in (2) gives
\[
 M(h_n)-M(h_{n+1})\le C h_n^{5/2}e^{-d^2/(2h_n)}.
\]
The limiting mass at zero excess is zero since the minimum is strict. Summing the geometric sequence in inverse time proves (1).

For either flank let \(m(r)=\max_{0\le u\le r}[f(x_0\pm u)-t_0]\). Its spatial interval has mass at least \(c r\) and is within the well sublevel at \(h=m(r)\). With \(L=\log(1/r)\), (1) yields the necessary bound
\[
 L-\frac{d^2}{2m(r)}-\frac52\log L\ge -C.                     \tag{3}
\]
Indeed first \(m(r)\gtrsim L^{-1}\). If \(m(r)\le d^2/L\), use \(\log(1/m)\ge\log L-C\); in the other case the left side tends to \(+\infty\). The reflected statements hold for a leftmost birth. Thus a monotone local well with a one-flank ceiling \(m(r)\le (d^2/2)/(L-\beta\log L+O(1))\) for \(\beta<5/2\) cannot arise as a positive-time outer birth. This remains valid when complicated components are already present on the other side. The argument does not bound the hull width of a well from its mass when both flanks can oscillate.

## 2. Extra wells are incompatible with earlier or equal continuation on both old collars

Suppose the early sublevels were connected up to a first separated birth time \(t_0\in(0,T)\), with unique zero at the origin. Round 84 gives a real Brownian first-exit interpretation before \(t_0\), and early enthalpy uniqueness identifies it with the symmetric canonical Gaussian first-exit boundary. Let \(b\) be its radius at \(t_0\), and \(F_*\) its increasing spatial inverse. Use only a short finite-time canonical continuation. Suppose the actual activation time on both old outward collars obeys
\[
 f(x)\le F_*(|x|)\qquad (b<|x|<b+\eta).                       \tag{4}
\]
No modulus assumption is needed on the newborn wells, whose minima at time \(t_0\) are separated from the old interval. Under (4) such a separated birth is impossible.

At time \(t_0-\) the actual and canonical heat remainders agree. For all sufficiently short times \(t_0+h\), each canonical newly activated mean is in the stated collars. The actual times there are at least \(t_0\) because all prior sublevels were inside \([-b,b]\). Condition (4) ensures that each canonical newly activated mean is activated in the actual graph no later, and consequently has an actual heat age at least its canonical age. Choose a fixed exterior point \(y\) far enough to the right of the union of both active spatial supports over the short slab, with \(|y-x|^2>h\) for every newly activated mean. This is possible by compactness from the high tails. For such \(y\), \(p_s(y-x)\) increases in \(0<s\le h\), as follows by differentiating in \(s\). The actual deducted past potential from the same spatial means is therefore at least the canonical one. There are additionally actual newborn well sources, with positive mass at some level strictly before \(t_0+h\) by continuity, all disjoint from those canonical means for small \(h\). They contribute a strictly positive heat density at \(y\).

Both remainders at \(y\) must be zero: the canonical one by the real exit and the actual one by signed exterior support. Subtracting their identities from their common initial remainder contradicts the strict extra past contribution. This is a two-sided ordering restriction for the old collars, not a condition on the old or newborn slopes. In particular, if the actual \(f\) is analytic through both old endpoints, the early agreement and the positive-time analyticity of the canonical inverse give equality on both outward collars, which is impossible if any separated birth occurs. The same implication holds for a common subtraction-stable quasianalytic class.

## Scope of spatial continuation

At fixed time the signed heat remainder is real analytic on each open spatial interval separated from zero-age source locations, namely away from the graph level interfaces. Exterior zero extends throughout that exterior component. A disconnected newly active island generally has its own zero-age endpoints; spatial analyticity alone does not continue this zero through those endpoints into its interior or into the gap on its other side. Bound (1) instead uses exterior equality arbitrarily close to the new edge, and Section 2 uses exact exterior equality plus the heat-age ordering of both old collars. Arbitrary continuous multiple-well early graphs without additional hypotheses remain open.

## Local calculation

Run `/Users/zzh/.venvs/main/bin/python check_multiwell.py`. Python 3.13.18 with SymPy 1.14.0 verifies the three exact identities for the heat-age derivative and adjacent bins. With 70 decimal digits, mpmath 1.3.0 checks the normalized half-line wall heat and the bin mass for the comparison profile at distance one. At \(h=2^{-24}\) these give \(0.3989422090650152\) and \(0.3934694306674684\), respectively, approaching \((2\pi)^{-1/2}\) and \(1-e^{-1/2}\). The JSON records these deterministic calculations. The exact bounds follow from the preceding proof.
