# Two frequencies at the middle-thirds maximum set

Let E be the middle-thirds Cantor set in [0,1], d=log(2)/log(3), D=1−d. Assume the normalized bounded continuous Cantelli graph is exact and that, in a punctured collar of E, v=H−f satisfies A₋L^(−β)≤v≤A₊L^(−β), L=log(1/dist(x,E)), with β>1 and 0<A₋≤A₊. E consists of highest points. Write α=2/(β+1),

    J(q,A)=min_{u>0} [qu+A/(2u^β)].

The comparison depends only on A₊/A₋. Fixed affine copies obey the same conclusions by translation, rescaling the frequencies, and allowing arbitrarily small strict changes to the asymptotic coefficient bounds to absorb an additive constant inside L.

## A weighted test and its rates

For a fixed sufficiently small collar U as in R186, let χ_U(k)=∫_U exp(ikx−k²v(x)/2)g(x)dx. The spatial strip argument in proofs/186-nonradial-cantor.md, “Localization”, gives |χ_U(k)|≤C exp(−ck). Set k=π3^m, K=k^α, q=k²/2. Fix 0<ε<1/2 and τ≥0. Test

    T=Im χ_U(k)+ε exp(τK) Re χ_U(2k).

Then |T|≤C exp(−c'k), since K=o(k). With s=sin(kx), the integrand apart from the positive source density is

    exp(−qv)[s+ε exp(τK−3qv)(1−2s²)].                                 (1)

For all τ the positive contribution has lower exponential rate at least exp[−(J(D,A₊)+o(1))K]. Indeed each depth-m cylinder starts at an even integer multiple of 3^(−m). In its normalized [0,1] coordinate choose the same descendant [2/27,1/9], where sin(πx)≥c₀ and cos(2πx)≥c₀. In each much deeper descendant take a central fixed fraction of the next gap. At depth n their combined length is comparable to (2/3)^n and L=n log 3+O(1). Both terms in (1) are positive there. For arbitrary u>0 choose n log 3=uK+O(1) and apply the upper variance bound to the first term; minimizing gives the claimed rate.

For the negative contribution define

    N_c(τ)=min_{u>0} [2u+max(A₋/(2u^β),(u+τ)/3)],
    u₀(τ)=(3A₋/(2τ))^(1/β),
    N_d(τ)=inf_{u≥u₀(τ)} [D u+2A₋/u^β−τ]    (τ>0),     N_d(0)=+∞.      (2)

Its absolute mass is bounded by exp{−[min(N_c(τ),N_d(τ))−o(1)]K}. To see this, separate two cases.

* If s<0, x lies in a gap created by depth m or in an exterior collar. There are O(2^m) such intervals, and a nearest point of E is an endpoint with sin(kx)=0. Thus |s|≤k exp(−L). The region k exp(−L)≥1/4 contributes O(exp[−c k²/(log k)^β]) even after exp(τK). Else 1−2s²≥1/2 and a negative value of (1) requires 3qv≥L+τK−log k−C. Its negative part is bounded by exp(−qv)k exp(−L). Changing variables to L in each gap adds exp(−L)dL up to a fixed factor. The first bound in (2) now follows from qv≥A₋k²/(2L^β) after scaling L=uK; polynomial factors in k cost exp(o(K)).
* If s≥0 and τ=0, (1) is nonnegative, since 1−2s²≥−s on [0,1]. If τ>0 and (1) is negative, then 1−2s²<0 and s>1/√2, so 3qv≤τK+C. This forces L/K≥u₀(τ)−o(1). Its negative part is at most ε exp(τK−4qv). For t small, |{x:dist(x,E)<t}|≤C t^D by the depth-n interval cover at scale t. Divide L into unit intervals and use the lower variance bound in each. The resulting exponential rate is at least N_d(τ). In both Laplace estimates the range of u can first be restricted to a large fixed compact upper range since the volume factor decays exponentially in L; at u↓0 the term u^(−β) diverges. These bounds require no equality of the variances at different addresses.

Thus the sufficient exclusion is

    J(D,A₊)<min(N_c(τ),N_d(τ))     for some τ≥0.                         (3)

## Explicit improvement with a fixed multiplier

For τ=0, the old optimizer for 2u+A₋/(2u^β) is u_*=(βA₋/4)^(1/(β+1)). The two expressions inside the maximum cross at u_b=(3A₋/2)^(1/(β+1)). When β>6, u_*>u_b and the minimum in (2) occurs at the crossing, giving N_c(0)=(7/3)u_b>J(2,A₋). When β≤6 this test gives N_c(0)=J(2,A₋). Since J(q,A) is proportional to (q^β A)^(1/(β+1)), for β>6 the explicit sufficient condition is

    A₊/A₋ < (2/D)^β (6/β)[7β/(6(β+1))]^(β+1).                         (4)

The factor beyond (2/D)^β is strictly greater than one by the strict inequality for the minimum. It is about 1.0404 at β=8 and 1.3103 at β=12.

## An improvement below and at β=6

Put t=τ/u_* and use u_* as the unit of u. Then A₋K/(2(u/u_*)^β u_*^β), divided by K u_*, has coefficient 2/β. In these units J(2,A₋)/u_*=2+2/β. The old minimizer is at u/u_*=1. Strict improvement in N_c is equivalent to t>6/β−1 (with t≥0): at that point the new second entry (1+t)/3 exceeds 2/β. Compactness of bounded sublevels and uniqueness of the old minimizer make the increase strict. In the same units

    N_d(τ)/u_* ≥ D(6/(βt))^(1/β)−t     (t>0).

For every 119999/20000≤β≤6 choose t=1/100000. Then t>6/β−1, and (6/(βt))^(1/β)≥100000^(1/6)>27/4. Since 2^11<3^7 we have D>4/11, so N_d(τ)/u_*>27/11−1/100000>12/5>2+2/β. Both bounds in (3) strictly exceed the R186 rate. In particular there is a nonempty open interval of coefficient ratios strictly beyond (2/D)^β for each of these β. One fully specified new threshold is

    R(β,τ)=(2/D)^β [min(N_c(τ),N_d(τ))/J(2,A₋)]^(β+1).                (5)

Here τ=t u_*; the right side is independent of A₋ by homogeneity. No optimality for this interval of β or the two-frequency family is asserted.

## Localization on the address tree

Let B(x)=v(x)[log(1/dist(x,E))]^β in the ambient collar, with no global ratio assumed. For a ternary address a∈E let G_n(a) be the middle gap of its depth-n cylinder. Let B_n⁻(a),B_n⁺(a) be the infimum and supremum of B throughout G_n(a), for sufficiently large n that those gaps lie in the collar. These extended-real functions are constant on depth-n cylinders, whose intersections with E are clopen. Let R_* be any of the strict ratio thresholds above, or the R186 threshold when appropriate.

Suppose the set of addresses with

    0<liminf_n B_n⁻(a)≤limsup_n B_n⁺(a)<∞,
    limsup_n B_n⁺(a) < R_* liminf_n B_n⁻(a)

is not of first category in E. For rational 0<p<q with q/p<R_* and an integer N, the set ∩_{n≥N}{B_n⁻≥p, B_n⁺≤q} is closed. The union of these sets contains every indicated address. Some one of the closed sets has nonempty interior, hence contains a full cylinder. Descend sufficiently far inside it and choose a smaller cylinder strictly separated from both of its endpoints (append a 0 followed by a 2 in standard ternary digit notation, after any intermediate descent). Both gaps bordering this smaller cylinder lie in the original cylinder and belong to depths ≥N. In a sufficiently small fixed neighborhood of the smaller Cantor copy the distance to that copy equals distance to E, and all nonzero-variance locations lie in controlled gaps. The local result therefore applies with coefficients p,q. For an exact graph the indicated addresses necessarily form a meagre set; addresswise convergence of the two bounds to the same finite positive limit at a nonmeagre set of addresses suffices regardless of the global amplitude ratio. The argument requires control through each indicated gap, including near its endpoints.

## Exact onset of strict improvement for the rates in (2)

For 1<β<6 let t₀=6/β−1 in the normalized units and V₀=(6/(6−β))^(1/β). The coarse rate improves precisely when t>t₀. The deep rate is nonincreasing in t, since its lower integration cutoff decreases and the additive term is −t. At t₀ its value is at most D V₀+t₀/3 by evaluating at the cutoff (8/(β V₀^β)=4t₀/3). This is larger than 2+2/β only if D V₀>7/3. Conversely if this strict inequality holds, the unconstrained minimizer (8/D)^(1/(1+β)) is smaller than V₀: indeed V₀>7/(3D)>√(8/D), using D<1/2. The deep rate at t₀ then equals its endpoint value. It remains strictly above the old rate for t just above t₀ by continuity, which proves sufficiency. Equality gives no strict gain.

The function V₀ is strictly increasing in β, by the strictly increasing power series for −log(1−z)/z, 0<z<1. It starts below 7/(3D) at β=1 and diverges towards 6. Thus the rates (2) certify strict improvement precisely for β>β_c, where β_c is the unique solution in (1,6) of

    (6−β_c)[7/(3D)]^β_c = 6.

Bisection using mpmath with 50 decimal places gives approximately 5.99990602333901570618. This boundary concerns the estimates (2). It does not establish optimality among other estimates even for the same tests.
