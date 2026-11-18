# A spatial derivative at an exposed endpoint of a regular terminal cluster

Write p_b for the centered Gaussian density of variance b, g=p₁. Suppose ∫p_{v(x)+s}(y−x)g(x)dx=p_{H+s}(y) for all s>0 with H>0, v≥0. Let E be a compact d-regular subset of the zero set of v with 0<d<1, Q=1−d, and write a=max E, r(x)=dist(x,E). If c r(x)^p≤v(x)≤C r(x)^p throughout a punctured collar of E, then all exponents

    2(1−d)/(2−d) ≤ p < 2

are impossible, with arbitrary coefficients between the two comparison constants. There is no restriction on other zeros outside the collar or the behavior at infinity. The derivative argument below proves this for p<4/3, including its lower endpoint. For p≥4/3 (indeed all 1<p<2), the previously demonstrated one-wing fractional inverse-moment bound of round 183 applies directly at a, since on its right r(a+u)=u and v(a+u)/u²→∞. That bound gives ∫ v(a+u)^(−s)du<∞ for all 1/2<s<1, contradicted by choosing ps≥1.

**Uniform budget for a first spatial derivative.** For each s>0 differentiate the smoothed identity at a. The derivative kernel for a source at x is (x−a)p_{v+s}(a−x)/(v+s). All x<a contribute nonpositively. On x>a+δ its absolute value is uniformly O(δ⁻²), by maximizing u b^(−3/2)exp(−u²/(2b)) over b for u≥δ. On a<x<a+δ its positive integral is at most C∫₀^δ u v(a+u)^(−3/2)du, finite for p<4/3. The exact derivative p′_{H+s}(a) stays bounded as s↓0. Fatou therefore gives the finite left-wing budget

    ∫_{x<a, v(x)>0} (a−x)v(x)^(−3/2)exp(−(a−x)²/(2v(x)))g(x)dx < ∞.          (1)

This conclusion needs only the right-wing integral or a uniform positive derivative bound on that side.

**Local tube shells on the other side.** Fix η>0 small depending only on the regularity constants so μ(E∩[a−h/2,a−2ηh])≥c₁h^d for all small h. This follows from the mass lower bound at a with radius h/2, the mass upper bound with radius 2ηh, and E⊂(−∞,a]. For t=o(h) a maximal t-separated family of these points has ≥c(h/t)^d members by upper regularity. The disjoint radius-t/3 intervals around them lie in [a−h,a−ηh] and in {r≤t} and have total length ≥c₂h^d t^Q. Conversely the portion of [a−h,a−ηh] in {r≤θt} has length at most C h^d(θt)^Q, from a maximal θt-separated family on E∩[a−2h,a] and its disjoint θt/3 balls: each has mass ≥c(θt)^d and their union has mass ≤C'h^d. Choose a fixed 0<θ<1 so that Cθ^Q≤c₂/2. This proves

    |{a−h≤x≤a−ηh, θt≤r(x)≤t}| ≥ c₃ h^d t^Q.                               (2)

In (2) set t=h^(2/p), so t=o(h) for p<2. At all such points v is between positive multiples of h² and a−x is between ηh and h, making the integrand in (1) at least c₄h^(−2), since g is bounded below near a. Take h_n=h₀γⁿ with 0<γ<η so the annuli are disjoint. Formula (1) forces convergence of Σ h_n^(d−2+2Q/p), which fails if d−2+2Q/p≤0. This proves the claimed threshold including equality.

The same proof applies to a collar comparability profile cV(r)≤v≤CV(r). It suffices that the right-wing integral in (1)'s proof be finite and that some t_n=o(h_n) satisfy V(r) comparable to h_n² uniformly for θt_n≤r≤t_n. Then Σ h_n^(d−2)t_n^Q must converge. At the critical exponent p=2Q/(2−d), the profile V(r)=r^p(log(1/r))^β requires β>p/Q=2/(2−d): choose t_n=h_n^(2/p)(log(1/h_n))^(−β/p). More generally the same necessary series decides secondary logarithms when the first logarithmic exponent is critical.
