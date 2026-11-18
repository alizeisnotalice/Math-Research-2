# First time derivative of the signed resolvent in a uniform strip and far cone

Let 0≤f≤H satisfy ∫exp(λx−λ²f(x)/2)dγ=1 for every real λ. Let y≥0, u=x−y, z=√(2q) and k_q(u)=sinh(zu)/z with k₀(u)=u. The right-wing resolvent is

    R(q,y)=2∫₀^∞ e^(−q f(y+u)) k_q(u) g(y+u)du.

All derivatives in q below are justified by bounded f and Gaussian domination. Reflection gives the same results on the left. Write μ_y for the probability measure on u>0 proportional to u g(y+u), A for its normalizing integral, m=Eμ u, m₂=Eμ u², C=y+m, and d=log(m/A).

## Uniform strict entropy margin at zero

The KL divergence relative to the Gaussian input law of x=y+u is Eμ log(u/A)≤d. The entropy variational inequality for the exact tilt yields λC−λ²Eμ f/2≤d. Optimization gives Eμ f≥C²/(2d). Integration by parts gives m₂=2−ym. Writing V=Varμ(u), differentiation yields m'=-V, C'=1−V=mC−1, and d'=2C'/m. The moment recurrence gives m₃=3m−y m₂. By Cauchy–Schwarz, m₂²≤m m₃. Substituting m₂=2−ym and cancelling gives 3m²+2ym≥4. Hence mC=m²+ym≥(4+ym)/3.

Consequently mC=2−V≥4/3, C'>0 and

    (3C²/4−d)'=(C'/m)(3mC/2−2) ≥0.
    3C²/4−d ≥ δ := 3π/8−log π >0.

For a purely rational verification of the last strict sign, π>3.14 and π<22/7; the exponential at 3·3.14/8 exceeds 22/7 already on taking the first five nonnegative Taylor terms (orders 0 through 4). Hence Eμ f≥(2/3)(1+δ/d). Since d increases and m₂ decreases (m₂'=−Cov(u²,u)<0), for every y≥0

    Eμ f ≥ κ m₂/3,    κ=min{1+δ/d(1), 2/m₂(1)}>1.                 (1)

## All spatial cuts in a positive-q strip

For each q≥0 the ratio k'_q(u)/(u³/3) is increasing in u by its nonnegative power series. The probability measure proportional to u³g(y+u) is an exponentially decreasing tilt of the same measure at y=0. The monotone covariance identity bounds the ratio of the two integrals by its value at y=0. Direct integration gives

    ∫ k_q(u)g(u)du /∫u g(u)du = ∫₀¹exp(q(1−s²))ds

where the left integrals are on u>0; this follows also by Gaussian completion of the square. Differentiating and using ∫₀¹(1−s²)ds=2/3 gives

    ∫ k'_q(u)g(y+u)du ≤ exp(q)∫ (u³/3)g(y+u)du.

The negative term in R' has magnitude at least 2exp(−qH)∫f u g(y+u)du since k_q≥u. Thus R'(q,y)<0 for every y≥0 as long as

    0≤q<log κ/(H+1).                                                (2)

## A cone reaching q proportional to squared distance

Set c_H=min(1/4,1/(32H)); normalization gives H≥1. There is Y_H finite such that

    y≥Y_H,   0≤√(2q)≤c_H y    implies    R'(q,y)<0.                (3)

The cutoff Y_H can depend only on H. Choose z₀=√(log κ/(H+1)); (2) covers z≤z₀. For z₀≤z≤c_H y define D_q=k'_q/k_q. From coth(v)−1/v≤1 we get D_q≤u/z. Also k'_q≤u³exp(zu)/3 by differentiating k_q=∫₀ᵘcosh(zv)dv and using sinh(zv)≤zv exp(zv).

The exact transform at λ=2y gives γ{x>y,f(x)≤1/2}≤exp(−y²). For 1/y≤u≤2/y the total Gaussian mass is ≥c exp(−y²/2)/y for y≥1. Removing the exceptional low-f set thus leaves at least half this mass for all sufficiently large y. On the remaining portion of this interval D_q≤1/4 for y≥8/z₀, so e^(−qf)k_q(f−D_q) integrated there supplies a positive reserve

    N ≥ c y⁻² exp(−y²/2−Hz²/2).

To bound contributions with the opposite sign, when 0<u≤z/4 they only occur on f<1/2 because D_q≤1/4. Their magnitude is at most C z³exp(z²/4−y²) by the exceptional-set bound. When u>z/4, discard f and use

    ∫_{u>z/4} k'_q(u)g(y+u)du
       ≤ C exp(−y²/2)∫_{u>z/4}u³exp(−(y−z)u−u²/2)du
       ≤ C exp(−y²/2−yz/16)

since z≤y/4. Both error bounds divided by N tend to zero uniformly in z₀≤z≤c_H y. For the second one, Hz²/2≤yz/64; for the first, (H/2+1/4)z²≤y²/8 by the choice of c_H. Thus N dominates, proving (3). No monotonicity of the graph or tail-gap hypothesis enters.

## Every fixed derivative order in a far cone

For each integer n≥1 put a_n=(2^{1/n}−1)/4 and c_{H,n}=a_n/(4H). There is Y_{H,n}<∞ such that

    y≥Y_{H,n},   0≤√(2q)≤c_{H,n} y    imply    (−1)^n ∂_qⁿ R(q,y)>0.      (4)

The same region works for all orders 1 through n after taking the maximum of their cutoffs. This gives the higher-order signs at q=0 for large cuts without any global higher-order entropy bound.

To prove (4), the convergent product k_q=u∏_{l≥1}(1+2qu²/(π²l²)) implies 0≤k_q^{(j)}/k_q≤D_q^j for j≥1. In fact differentiation restricts the j indices in the expansion of D_q^j to distinct ones, and the claim follows by taking limits of polynomial products on compact q-sets. Also D_q≤D₀=u²/3 and D_q≤u/z when z>0. Put v(z)=max(√(3a_n),a_n z). For every u≤v(z) one of the two upper bounds yields D_q≤a_n, including at z=0. Thus on f≥1/2 and u≤v(z), Leibniz and the triangle inequality give

    (−1)^n ∂_qⁿ(e^(−qf)k_q) ≥ e^(−qf)k_q f^n [2−(1+2a_n)^n].

The bracket is strictly positive since 1+2a_n=(1+2^{1/n})/2<2^{1/n}. The same positive packet 1/y≤u≤2/y, f≥1/2 supplies N_n≥c_n y⁻²exp(−y²/2−Hz²/2) for all sufficiently large y.

Throughout the half-line the absolute derivative is at most k_q(H+D_q)^n≤u exp(zu)(H+u²/3)^n. In the inner region u≤v(z), the remaining set f<1/2 has contribution at most C_{H,n}(1+z)^{2n+1} exp(−y²+zv(z)), by the same exact tilt at 2y. In the outer region u>v(z), since z≤y/4, integration of the polynomial against the Gaussian density gives C_{H,n}exp(−y²/2−y v(z)/4). Divide by N_n: for the outer term, Hz²/2≤y a_n z/8≤y v(z)/8 and v(z)≥√(3a_n)>0. For the inner term, zv(z)≤√(3a_n)z+a_n z² and z≤a_n y/(4H), with a_n≤1/4 and H≥1. Thus its relative logarithm is bounded above by −y²/4 for sufficiently large y, uniformly in the stated cone after absorbing polynomial factors. Both relative errors vanish, proving (4).
