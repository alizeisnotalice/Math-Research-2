# Arbitrary coefficient contrast for subquadratic spatial scales on the middle-thirds set

Let E be an affine middle-thirds Cantor set of highest points with a fixed collar U containing no other zeros; write r=dist(x,E), d=log 2/log 3, Q=1−d. Suppose on the punctured collar c V(r)≤v(x)≤C V(r) with arbitrary constants 0<c≤C<∞. Assume V is positive continuous near 0, tends to zero, V(r)/r²→∞, and for some q,c₀>0
$$
 V(s)/V(t)\ge c_0(s/t)^q \quad(0<t\le s\le r_0).
$$
Then the exact Gaussian output identity is impossible. No radial equality, small coefficient oscillation, or regularity in the gap address is required. This covers V(r)=r^p for 0<p<2 and quadratic opening times any positive power of log(1/r), as well as iterated logarithms tending to infinity multiplying r². For each of the specified logarithmic profiles the displayed inequality follows from the logarithmic derivative approaching 2 with respect to log r. Round 68 used Fourier freezing under a more restrictive near-constant coefficient description; here the ternary semicircle phase controls the sign directly.

Use [0,1] coordinates on the hull; a fixed affine scale is absorbed in k and the constants. At frequency k=π3^m let t=t_k be a solution of k²V(t)=1 tending to zero (continuity suffices). Then kt→0. By the scaling bound, any r≤t in the collar has V(r)≤1/(c₀k²). Within each depth-m interval choose its fixed interior subcylinder [2/9,1/3] in rescaled coordinates. All these subcylinders have sin(kx)≥sin(2π/9). Select n with 3^(-n) comparable to t and at most t, so n≥m+2 eventually. Counting all their depth-n descendants and retaining the middle third of each next gap gives total length ≥c₁ t^Q. Damping there is bounded below by a fixed positive constant and the input density by a positive constant. Thus the positive phase contributes at least c₂t^Q.

The negative part is restricted to the gaps removed through level m and the two outer collars. There are O(k^d) such intervals; |sin(kx)|≤kr on their negative part. The integrated contribution in each is at most a fixed multiple of
$$
 k\int_0^{r_0} r\exp(-c k^2 V(r)/2)\,dr=O(k t^2).
$$
Use the trivial bound for r≤t and the assumed power scaling for r≥t, followed by r=t z. Compact collar pieces separated from E contribute exponentially little. Consequently the ratio of the negative upper bound to the positive lower bound is O(k^(d+1)t^(2−Q))=O((kt)^(1+d)), tending to zero independently of C/c. Also t≥c₃ k^(−2/q) by the same scaling bound with fixed s, so the lower bound is at least a power of 1/k. The compact source density from U must have exponential Fourier decay by the strip argument of round 186, a contradiction.

The argument admits distinct bounds on gaps. For example suppose v≥c r^P on the entire punctured collar and every depth-n central gap has a subset of relative length at least c₄ where v≤C 3^(−np). If 0<p≤P and 4/P>1+d+2Q/p, exactness fails. The latter implies p<2. Use n with 3^(−n) comparable to k^(−2/p) for the positive subcylinders (they then lie strictly deeper than level m); the negative cost is O(k^(1+d−4/P)). The subset may change arbitrarily with its gap. This extension needs no upper bound on the rest of each gap.
