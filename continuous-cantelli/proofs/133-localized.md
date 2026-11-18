# Local mass and the large-parameter signed resolvent on fixed cuts

Let X be standard normal with density g, let f≥0 be bounded and continuous, and assume the exact tilt ν_z(dx)=exp(zx−z²f(x)/2)γ(dx) has total mass one for every real z. Fix 0<a≤b and assume f>0 on (0,2b+η] for some η>0. No spatial condition is imposed beyond that interval, on the other wing, or on the rate of approach to the origin. Let K(q,y)=exp(−zy)/z for y>0, z=√(2q), and let R be the signed heat resolvent of round 132. For fixed 0<r<a and integer n≥0,

    (−∂q)ⁿ R(q,y) / (−∂q)ⁿ K(q,y) = 1−ν_z((0,r))+o(1),

uniformly in a≤y≤b as z tends to infinity. Reflection covers negative cuts. In particular the normalized negative part tends to zero. If a positive fraction of the tilted mass stays outside (0,r) uniformly for large z, the fixed-order sign is strictly positive throughout the compact interval. The statement applies even when the local mass oscillates.

**Proof.** Denote the monic free derivative polynomials of round 132 by P_i. Divide the marked Leibniz sum by the free derivative. The integrand with respect to γ is

    exp[z(y−|y−x|)−qf] Σ_{j=0}ⁿ binom(n,j) (z²f)^j P_{n−j}(z|y−x|)/P_n(zy).

For z≥1 the polynomial sum is bounded by C(1+|x|)^n(1+zf)^n, uniformly on the stated cuts. On x≤0 the phase before subtracting qf is zx, and on x≥2b+η it is at most −zη. Absorb (1+zf)^n into exp(nzf) for large z. These two regions contribute O(1/z) and O(exp(−zη)), respectively. On [δ,2b+η] the positive minimum of f dominates every linear term in z, for every fixed 0<δ<r.

On (0,δ) the exponential is exactly the density of ν_z relative to γ. The part with zf≥4δ contributes o(1), even including its polynomial weights, since zx−(q−nz)f≤−zδ/2 for large z. Elsewhere zf≤4δ and y−x≥a−δ. The j=0 term is 1+O(δ+1/z) by monicity and y≥a (take δ<a/2), and all j≥1 together are O(δ) for 0<δ<min(1,a/2). The constants depend on a,b,n only. The ν_z mass is at most one. Finally ν_z([δ,r]) tends to zero. Taking the upper limit as z tends to infinity and then letting δ decrease to zero proves the assertion.

**Explicit critical profiles.** Write ℓ=log(1/x). Suppose for all sufficiently small positive x,

    f(x) ≥ x²/[2(ℓ+β log ℓ+c+o(1))].

For β<1/2 the normalized derivative tends to one. For β=1/2 its lower limit is at least 1−exp(c)/√2, uniformly on the compact interval. In particular c<log(2)/2 implies eventual strict positivity at every fixed order. To see this use round 28's one-sided saddle estimate on the upper majorant of the local mass, replacing c by c+ε and then letting ε decrease to zero. If equality holds in the contact formula up to o(1) in the denominator, two-sided comparison gives convergence to 1−exp(c)/√2 in the critical case. The contact asymptotic for mass alone was obtained in rounds 28–29. Here it determines the leading term of every fixed signed resolvent derivative at macroscopic cuts. At full saturation the leading term vanishes and this estimate does not decide the sign.

The accompanying integration uses u=zx and centers at u=2B, with B=ℓ+β log ℓ+c and log z=ℓ+log(2B). It checks the saddle normalization for two examples. It is a finite-window numerical check of the mass asymptotic. The proof of the resolvent statement uses the exact polynomial identity and localization above.
