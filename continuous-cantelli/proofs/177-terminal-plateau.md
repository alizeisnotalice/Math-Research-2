# A terminal plateau cannot open along a monotone flank or a flank with bounded rescaled Hölder seminorm

Let a continuous bounded graph f be exact in the normalization X~N(0,1), E exp(zX−z²f(X)/2)=1. Let H=sup f, let a be an attained maximum and set v(u)=H−f(a+u). Suppose v(u)=0 on [−r₀,0]. Write w(u)=v(u) on the right and assume w has positive values in every (0,r). Each of the following hypotheses contradicts exactness.

1. The function w is nondecreasing on (0,r₀). Continuity at 0 is the only modulus assumption here.
2. There are 0<β≤1 and K<∞ such that for all sufficiently small r the functions w(rt)/r² have β-Hölder seminorm at most K on 1/2≤t≤3.

Reflection gives the other orientation. The second hypothesis holds when w is locally absolutely continuous off zero and |w'(u)|≤Cu almost everywhere. It holds in particular if the full graph is locally C¹ with Lipschitz derivative, since its derivative at the plateau edge vanishes. Other maxima, including intervals and accumulation points away from a, are allowed.

## The phase at the edge

R021, proofs/021-contact-phase.md lines 3–8, gives a common volume fraction θ on all fixed rescaled intervals for {v(u)≤u²}, and vanishing relative mass on each fixed annulus for {η≤v(u)/u²≤1}, for every 0<η<1. The plateau forces θ=1. Consequently for every fixed 0<c<C and η>0, the length of {cr<u<Cr : w(u)≥ηu²} is o(r).

In case 1 a value w(r)≥ηr² would imply w(u)/u²≥η/4 throughout [r,2r]. In case 2 a value w(r)≥ηr² would imply w(rt)/r²≥η/2 on [1,1+d], where d>0 depends only on K,β,η. Either implication contradicts the preceding length estimate along a sequence r→0. Hence q(u)=w(u)/u²→0 uniformly as the upper endpoint of (0,r) goes to zero.

## A power gain for the area of a general flat flank

This step only requires a plateau on the left and q→0 on the right. Let p_s denote the Gaussian density of variance s, g(u)=p₁(a+u), K_s(z)=∫₀ˢp_t(z)dt and K₀=0. For a fixed small dyadic cutoff R₀ define

    F(y)=∫₀ᴿ⁰ g(u) K_w(u)(y−u) du,

where K_w(u) in this expression means the kernel K_{w(u)}. The full potential has distributional second derivative 2[p_H(a+y)−p₁(a+y)], by the heat equation and the exact spatial output law. It is therefore real analytic. All inputs |u|≥R₀ contribute an analytic function for |y|<R₀/4, by the separated complex Gaussian bound. The inputs −R₀<u≤0 give zero. Thus F is itself analytic near zero. In particular its odd remainder O(y)=F(y)−F(−y)−2yF'(0) is O(y³).

The prescribed nonzero log-power contact in R063 was used in proofs/063-opposite-side.md lines 129–134 to prescribe the leading term of O(y). The present cubic bound follows from the plateau and analyticity. The estimates at lines 136–144, especially (9), are estimates on the one active side and require no log-power comparison flank. Here are their details for this application.

Write ρ_j=2^(−j), R₀=2ρ_{j₀}, m_j=ρ_j^(−3)∫_[ρ_j,2ρ_j] w, e_j=sup_[ρ_j,2ρ_j] exp(−b/q), with value zero at q=0 and a sufficiently small absolute b>0. The numbers m_j are bounded and e_j→0. Let r=ρ_j and r/2≤y≤3r. Take C₀=16 and j sufficiently large that C₀r<R₀. In O(y), retain the positive integral of g(u)K_{w(u)}(y−u) on 0<u≤C₀r.

On this same interval the subtracted negative-observation and linear terms cost in absolute value at most

    C∫₀^(C₀r) [w(u)/u + r w(u)/u²] exp(−b/q(u))du.

Indeed K_s(−y−u)≤C√s exp(−u²/(2s)), and |K'_s(−u)|=2 Φ̄(u/√s). Decreasing b absorbs the powers of q. The bounded kernel derivative also justifies F'(0) under the integral. On C₀r<u<R₀ the third-order Taylor remainder of the odd difference costs at most

    C r³ ∫_(C₀r)^R₀ w(u)u^(−4) exp(−b/q(u))du.

Off the kernel cusp K_s'''=2p_s', and |u−y|≥u/2, which gives the stated bound after absorbing powers of q. On a shell of radius ρ≤r the two inner terms divided by r² are bounded by C[(ρ/r)²+ρ/r] e_n m_n. On shells of radius ρ≥C₀r the corresponding outer errors divided by r² are at most C(r/ρ)e_nm_n. The finitely many intermediate shell offsets obey the same geometric bound after increasing C. Therefore the retained positive integral differs from O(y) in absolute value by at most

    C r² Σ_{n≥j₀} 2^(−|j−n|) e_n m_n + C r³.                 (1)

The last term also covers any fixed separated analytic contributions. For u∈[r,2r] and 0<t≤w(u), a Gaussian N(u,t) lies in [r/2,3r] with probability at least 1/2 for all sufficiently small r, since sup_[r,2r]w=o(r²). Integrate the retained positive integral in y on that observation interval, restrict u to this shell and use inf g>0 there. Its lower bound is c r³m_j. Integrating (1) and using O(y)=O(r³) gives

    m_j ≤ C 2^(−j)+C Σ_{n≥j₀} 2^(−|j−n|) e_n m_n.             (2)

Fix 0<α<1. With W_j=min(2^(αj),2^(αN)), the ratio W_j/W_n is at most 2^(α|j−n|) uniformly in N. Choose J large enough that C sup_{n≥J}e_n times Σ_{l∈ℤ}2^(−(1−α)|l|) is <1/2. The finitely many n<J contribute a bounded quantity to W_j times the right side uniformly in N, since they give O_J(2^(−j)) for j≥J. Also sup_j W_j2^(−j) is bounded. The finitely many j excluded when deriving (2) are handled directly. Taking the supremum, absorbing the term with coefficient <1/2 and sending N→∞ yields

    ∫_[r,2r] w ≤ C_α r^(3+α)                                  (3)

first for dyadic r and then for every small r by covering with two dyadic shells. No smoothness or pointwise positive lower bound on w enters (1)–(3).

## From shell area to a pointwise upper bound

In case 1 monotonicity gives r w(r)≤∫_[r,2r]w, hence w(r)≤C_αr^(2+α). In case 2 put M=w(r)/r², which tends to zero. When M>0 there is a right interval from t=1 of length at least c M^(1/β) on which w(rt)/r²≥M/2. This gives ∫_[r,2r]w≥c r³ M^(1+1/β). Thus w(r)≤C r^(2+αβ/(1+β)). For K=0 the estimate follows immediately by integrating the constant rescaled profile. In either case there is a δ>0 with

    w(r) ≤ C r^(2+δ).                                          (4)

In case 2 the variance extended by zero through the plateau is locally β-Hölder. For comparable positive x,y, take r=x and use |w(x)−w(y)|≤C x^(2−β)|x−y|^β. For noncomparable points use w(u)=O(u²) and |x−y| comparable with the larger positive coordinate. The same bound applies across zero. Thus the local Hölder theorem in R031, proofs/031-first-hitting.md lines 9–16, forces sup_[0,r]w≥c r²/log(1/r) for all small r. This contradicts (4).

## First hitting with a monotone flank of arbitrary modulus

For case 1 only the analytic propagation from R031 is used. Its local disk bounds and chain, proofs/031-first-hitting.md lines 3–7 and 11–12, use the positivity of the potential, the plateau, and a uniform upper bound for w on the core. Hölder continuity is used in the separate hitting-point lower bound at line 9. We give a different lower bound.

Fix a small attained level h₀>0. Put h_j=2^(−j)h₀ and let r_j be its first hitting radius for j≥0. These radii decrease strictly to zero by continuity and nontrivial emergence. At j≥1 let R=r_{j−1}, Δ=R−r_j and x=R−Δ/2. The core potential S_R(y)=∫₀ᴿg(u)K_{w(u)}(y−u)du is twice the convention of R031; the fixed factor only changes constants. On its core w≤2h_j; on the segment −R/32≤y≤−R/64 the positive inputs have distance at least R/64, so S_R(y)≤C exp(−cR²/h_j). By R031's analytic interpolation chain,

    S_R(x) ≤ C exp[−c(Δ/R)^p R²/h_j]                           (5)

with a fixed p>2. This estimate is used even when its right side is only a constant. Monotonicity gives w(u)≥h_j on [r_j,R]. Integrating over |u−x|≤min(Δ/4,√h_j/4) and heat times h_j/2≤t≤h_j gives

    S_R(x) ≥ c √h_j min(Δ,√h_j).                               (6)

If Δ≥√h_j this lower bound has negative logarithm at most C(j+1), and (5) yields Δ≤C R^(1−2/p)[(j+1)h_j]^(1/p). If Δ<√h_j, the inequality already follows from Δ≤min(R,√h_j)≤R^(1−2/p)h_j^(1/p). Sum over all indices l≥j, using r_{l−1}≤R and ∑_{l≥j}[(l+1)h_l]^(1/p)≤C[(j+1)h_j]^(1/p). Since the gaps sum to R, this proves R²≤C(j+1)h_j. In particular log(1/R)≥c j for large j. On the other hand (4) at R gives 2h_j=w(R)≤C R^(2+δ). Combining yields 1≤C(j+1)R^δ, whose right side tends to zero by the exponential bound on R. This proves case 1 without Hölder regularity or an inverse-speed bound.

The proof uses only the local analytic final output and the R021 phase. It imposes no hypotheses on late temporal signs. It treats terminal intervals with the specified emerging flank and leaves isolated maxima and general continuous nonmonotone flanks open.
