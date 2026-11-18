# Uniform eventual positivity on an entire spatial wing

Let f be a bounded nonnegative continuous exact graph. Suppose f>0 throughout [a,∞), where a≥0, and write l=liminf_{x→∞}f(x), h=limsup_{x→∞}f(x). Assume

    l>1,        h−l<(l−1)².

For every fixed integer n≥1 there is q₀<∞ such that (−∂q)ⁿ R(q,y)>0 for all q≥q₀ and **all** y≥a. In particular this holds with an arbitrary amount of oscillation on compact intervals and any positive-side limit exceeding one. The reflected version only needs the corresponding left-wing hypotheses. No information about contact at the origin or about the other tail is needed. Integer differences D_{j,n} are consequently strictly positive throughout this wing for integers j≥q₀, by integrating the derivative over the n-dimensional unit cube. Combining with the round 122 far cone, any wrong sign of order n on this wing must lie inside a bounded rectangle in (q,y).

**Positive integrand region.** Use the exact right-wing formula of round 122 with u=x−y≥0. Let z=√(2q), K=sinh(zu)/z, and m=inf_{[a,∞)} f>0 by compactness and the positive lower tail. Choose 1<r<l and ε>0 so that r<l−ε and (r−1)²>h−l+2ε. Consider the signed derivative of exp(−qf)K divided by exp(−qf)K (for u>0). For u/z≤δ this ratio is at least fⁿ[2−(1+δ/m)ⁿ] by the product bound 0≤∂q^j K/K≤(u/z)^j of round 122. Fix δ>0 small enough that this gives a positive uniform lower bound. On δ≤u/z≤r, the ratio is

    (f−u/z)ⁿ+O(z⁻²).

To verify the approximation, differentiate exp(zu)/z and exp(−zu)/z separately using ∂q=z⁻¹∂z with u held fixed. Successive normalized derivatives of the increasing exponential equal (u/z)^j+O(z⁻²) on the stated compact ratio interval: the exact multiplier is a polynomial in u/z and z⁻², with the leading term (u/z)^j. The decreasing exponential has the factor exp(−2zu) relative to the increasing one and all normalized polynomials there are bounded. Leibniz gives the displayed formula. As x≥a+δz on this interval, f(x)≥l−ε for large z. Consequently the signed derivative is bounded below by c exp(−qf)K on the whole region u≤rz, with c>0 uniform in y. Beyond that region its absolute value is bounded by exp(−qf)K(H+u/z)^n where H=sup f.

**Uniform comparison when a≤y≤z.** The packet z+1≤x≤z+2 lies in the positive region for large z. There u≥1, K≥c exp(zu)/z and f≤h+ε. Integration against g(x) gives a lower bound

    c z⁻¹ exp[−zy−q(h+ε−1)].

On u>rz we have x≥rz+y and f≥l−ε. Complete the Gaussian square after multiplying by exp(zu) and use u/z≤x/z. The bad contribution in absolute value is bounded by

    C(1+z)^n z⁻¹ exp[−zy−q(l−ε−1+(r−1)²)].

Indeed the remaining Gaussian tail starts at x−z≥(r−1)z and its fixed polynomial moment obeys this bound by integration by parts. The ratio tends to zero uniformly.

**Uniform comparison when y≥z.** Integrate instead on 1/(1+y)≤u≤2/(1+y). There K≥u, g(y+u)≥c g(y), and f≤h+ε. This gives the positive lower bound c g(y)exp[−q(h+ε)]/(1+y)². On u>rz use exp(zu)g(y+u)=g(y)exp[−(y−z)u−u²/2]. The absolute bad contribution is at most

    C(1+z)^n g(y) exp[−q(l−ε+r²)−rz(y−z)].

Here (1+y)² exp[−rz(y−z)]≤C(1+z)² for z≥1, so the bad/positive ratio tends to zero uniformly. This proves the claim including arbitrarily large spatial cuts. The strict inequality in the oscillation hypothesis is used only in the first packet comparison; no assertion of necessity is made.
