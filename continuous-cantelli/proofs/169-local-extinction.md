# Local extinction time from a one-sided imaginary transform

Let γ(dx)=g(x)dx be standard Gaussian measure and let a bounded continuous function 0≤f≤H satisfy ∫exp(zx−z²f(x)/2)dγ=1 for all real z and hence for complex z. Write p_t for the centered Gaussian density of variance t, P_t for its semigroup, and

    Q_t=P_tδ₀−∫_{f(u)≤t}P_{t−f(u)}δ_u γ(du).

Here P₀ is the identity. The Fourier transform gives Q_t=0 for t≥H; the quadratic Taylor coefficient gives ∫f dγ=1. Let S be the closed support of the signed spacetime measure Q_t(dx)dt on ℝ×[0,H], and define E(x)=sup{t:(x,t)∈S}, with value 0 for an empty fiber. Thus for every s>E(x) the remainder is zero on a neighborhood of x throughout times above s, as a spacetime distribution. Compactness of the time interval is used in this characterization.

## Formula and continuity

For positive integers n put k_n=√(2n),

    M_n(x)=∫_{u>x} exp(nf(u)+ik_n u)γ(du),
    A(x)=limsup_n n⁻¹ log(max(1,|M_n(x)|)).

Then

    E(x)=max(f(x),A(x)),
    E(x)≤min(sup_{u≤x}f(u),sup_{u≥x}f(u)),
    osc_I E ≤ osc_I f for every interval I.                         (1)

In particular E is continuous and inherits any modulus of continuity of f. On every connected component of {E>f} it is constant. On a bounded component (a,b) this constant c satisfies f(a)=f(b)=c and f<c throughout (a,b). On a right-unbounded component (a,∞), c=f(a)=limsup_{u→∞} f(u). The reflected conclusion holds on the left. If the component is the entire line, both tail limsups equal c. Every one-sided record high (a point with f(x)=sup_{u≤x}f or f(x)=sup_{u≥x}f) satisfies E(x)=f(x). If f has locally bounded variation, E has locally bounded variation and dE=1_{E=f}df, by locality of the derivative of continuous BV functions on coincidence sets (equivalently, apply the local coarea formula to E−f on its zero set).

For the proof define F_n dx=∫₀ᴴ e^{nt}Q_t dt. The heat equation in distributions gives F_n″+k_n²F_n=2e^{nf}g−2δ₀. Explicitly, on either side of zero,

    F_n(x)=(2/k_n) Im{e^(−ik_n x)[M_n(x)−1_{x<0}]},
    F_n′(x)=−2 Re{e^(−ik_n x)[M_n(x)−1_{x<0}]}.                    (2)

One can verify (2) by direct differentiation of the outgoing sine integral: both it and F_n are integrable, their difference solves the homogeneous Helmholtz equation on the whole line, and no nonzero solution is integrable. Integrability of the outgoing expression at −∞ uses ∫exp(nf+ik_nu)dγ=1. At zero F_n is continuous with derivative jump −2. Either one-sided derivative suffices for its exponential rate. The right and left transforms add to 1, which gives the second bound in (1) for A, and trivially also for f. This enclosure is consistent with the slice-wise convex-hull theorem of round 58.

The part F_n from times t≤s has uniform absolute value ≤C_H e^{ns} and Lipschitz constant ≤2e^{ns}. Indeed ∫₀ᴴp_v(d)dv≤√(2H/π) and ∫₀ᴴ|p_v′(d)|dv≤1 uniformly in d (with the weak spatial derivative after integration). These bounds apply to each positive source and the initial kernel. Zero-age P₀ terms at a single time per input have zero spacetime mass by Fubini. Thus if s>E(x), F_n and either local one-sided derivative at x are bounded by these early estimates, yielding A(x)≤E(x). Also E(x)≥f(x): the distributional source (∂_t−½∂²_x)Q is −γ(du)δ_{f(u)}(dt) on positive times, extending Q by zero for t>H. Continuity of f puts every positive-height graph point in the support of this sink, hence in S.

Conversely, M_n(y)−M_n(x) has magnitude ≤γ([x,y])exp(n sup_[x,y] f) for x<y. Suppose max(f(x),A(x))<s<E(x). There is a neighborhood J of x and s'<s with f≤s' on J and A(x)<s'. The increment bound and (2) imply |F_n(y)|≤C exp(ns') for all sufficiently large n and y∈J. For every nonnegative smooth spatial test ψ supported there, the scalar signed compact-time measure ⟨Q_t,ψ⟩dt therefore has no support at times >s'. We used the signed moment fact that the limsup of the nth root of |∫vⁿdν| is the upper support endpoint for a nonzero signed measure on a compact subinterval of (0,∞), applied to v=e^t. For completeness, a strictly smaller growth radius would extend the Laurent series of its Cauchy transform analytically and with real values across an interval below the upper endpoint; integrating its imaginary part at positive imaginary height against tests and using the Poisson approximation to the identity forces the measure to vanish on that interval, a contradiction. Tests in J then imply S has no points of J×(s',H], contradicting E(x)>s. This proves the formula.

For x,y in an interval I the same increment estimate gives E(y)≤max(E(x),sup_[x,y] f) and its reverse. If E(y)>E(x), subtract E(x)≥f(x) to get E(y)−E(x)≤osc_[x,y] f. If E(x)>f(x), the bound in both directions implies E is constant near x. The topology statements follow, using the two-sided enclosure at an infinite endpoint. On an interval with M_n amplitude exp(nc+o(n)) along a subsequence and supremum of f at most m<c, the increment estimate also transfers a single spatial oscillation phase across the entire interval with absolute complex-amplitude error at most exp(nm); (2) turns it into the corresponding sine profile.

## A quantitative lower bound valid at every cut

If f is nonconstant, set

    B(x)=∫_{f(u)<1}∫₀^{1−f(u)}p_v(x−u)dv γ(du)>0.

Then

    B(x) ≤ ∫₁^{E(x)}p_t(x)dt,     E(x)>1 for every finite x.          (3)

To see this, the total occupation identity ∫₀ᴴQ_tdt=∫₀¹p_tdt follows by integrating the heat equation and using the common terminal spatial marginal γ. Equivalently, the second derivatives of both sides are 2(γ−δ₀), and both sides are integrable. Subtracting the explicit early occupation gives ∫₁ᴴQ_tdt=B(x)dx. Positivity of B uses nonconstancy and ∫f=1. The density B is continuous by the uniform bound on the integrated kernel and dominated convergence. For s>E(x), the late spacetime integral on a neighborhood of x is zero above s, so B(y)≤∫₁^s p_t(y)dt there whenever s>1, first as a measure inequality and then pointwise by continuity. If E(x)≤1, let s↓1 to get a contradiction. Letting s↓E(x) now proves (3). In particular the unique T(x)∈(1,H] solving ∫₁^{T(x)}p_t(x)dt=B(x) is a lower bound for both one-sided suprema of f. The bound covers points where f(x)=1 or lies below 1, and compact sets have a uniform strict lower bound for E−1.

For example write m=inf f<1 and h_±(x)=sup_{±u≥x}f(u)−1. Then both h_±(x)>0 at every finite x and

    liminf_{x→∞} x^(−2) log h_±(x) ≥ −m/[2(1−m)].                  (4)

If h_± does not tend to zero this is immediate by monotonicity and boundedness. Otherwise choose 0<2η<1−m and a fixed bounded interval on which f≤m+η. Integrating B over this spatial interval and 1−m−2η≤v≤1−m−η gives log B(±x)≥−x²/[2(1−m−2η)]−O_η(x+1). By (1),(3), B(±x)≤h_±(x)(2π)^(−1/2)exp(−x²/[2(1+h_±(x))]). Divide by x² and let η↓0. In particular when inf f=0 a Gaussian-rate upper decay of the suffix surplus on either entire wing is impossible.

There is a stronger local-in-space lower bound when an appreciable fraction of the real tilt is supplied at sublinear inputs. Put ν_x(du)=exp(xu−x²f(u)/2)γ(du) and c_x=ν_x(0≤u≤R(x)) for any positive R(x)=o(x), x→+∞. Then

    liminf_{x→∞} [x²(E(x)−1)−2log(1+c_x)] ≥ 0.                    (5)

By the score concentration of round 116, W=u−xf(u) tends to zero in ν_x-measure on this event, and so does f=(u−W)/x. At the upper age in B the exact square completion is

    p_{1−f}(x−u)/g(x)=(1−f)^(−1/2) exp(xu−x²f/2) exp[−W²/(2(1−f))].

For fixed K, integrate ages 1−f−v/x², 0≤v≤K. On the selected typical sources (u/x→0, f→0, W→0), the age-kernel ratio to p_{1−f}(x−u) is exp(−v/2)+o(1) uniformly on this v-range. Positivity and a diagonal selection of typical sets give x²B(x)/g(x)≥2(1−e^(−K/2))(c_x−o(1)); let K→∞. Along sequences on which a=x²(E(x)−1) stays bounded, x²∫₁^{E(x)}p_t(x)dt/g(x)=2(exp(a/2)−1)+o(1), uniformly on bounded a. Other sequences do not lower the liminf. This proves (5); the opposite wing follows by reflection.

## Each separated spatial interval has its own future transform

Let a<b satisfy E(a),E(b)≤h, allowing either end to be infinite and omitting its condition. For I=(a,b), s≥h and every real z,

    ∫_I exp(zy)Q_s(dy)=∫_{I∩{f>s}} exp(zx−z²(f(x)−s)/2)γ(dx).     (6)

Indeed for s>h the two finite endpoints have open collars with local extinction below s by continuity. There the spacetime distribution Q is zero after s, and f<s. Use a spatial cutoff equal to the indicator of I away from the collars, with transitions within them. Testing the heat equation from s to H against this cutoff times the backward heat solution exp(zy−z²(t−s)/2) leaves zero boundary terms. On each unbounded end use compact cutoffs tending to one: Gaussian spatial tails of Q_t in total variation integrated against exp(|z||y|) are uniformly small for t≤H, by Gaussian tails of the initial spatial law and the bounded heat ages. Thus (6) holds for s>h, including the strict future convention. To pass down to h, the positive past-source measure ∫_{f≤s}P_{s−f}δ_x dγ converges weakly from the right to its value at h by bounded continuous testing and dominated convergence (the endpoints f=h are included). At any s>0 the zero-age source is 1_{f=s}γ and all positive-age sources are Gaussian smoothings, so Q_s and the positive past-source measure have no spatial atoms. This allows the interval indicator, with uniform exponential integrability handling unbounded ends. The right side converges by dominated convergence. For a nonconstant graph h>1 at any finite endpoint by (3); for a whole-line interval (6) is the original future equality, with the known δ₀ at s=0.

In particular each connected component I of {E>h} obeys (6) on its own. The masses, first moments and second moments there agree with the future tests 1, x, x²−(f−s), respectively. The round-58 signed convex-hull argument applied to (6) gives supp(Q_s restricted to I) contained in the closed convex hull of the essential support of 1_{I∩{f>s}}γ. This statement applies to bounded components and rays. If I has no future mass at that level, the restricted slice vanishes by Fourier uniqueness. Conversely any discrepancy in (6) calculated from the graph at time s proves that at least one finite endpoint has E>s. These are graph tests and structural properties; they do not assert positivity of the local signed laws.
