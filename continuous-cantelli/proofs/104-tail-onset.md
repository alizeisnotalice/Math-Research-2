# Exact onset of the signed heat trace and a global exclusion for oscillating wings

Let g be the standard Gaussian density, γ its law, p_t the centered heat kernel of variance t. Let f be bounded, continuous and nonnegative and satisfy

    ∫ exp(zx−z²f(x)/2) dγ(x) = 1                 (z real).

Write H=sup f and Q_t=P_tδ₀−∫_{f(x)≤t}P_{t−f(x)}δ_x dγ(x). This signed measure vanishes for t≥H by the exact complex transform, as in R083. All supports in time below are supports of finite signed measures. For y>0 use the finite time trace

    dμ_y(t)=[p_t(y)−∫_{f(x)<t}p_{t−f(x)}(y−x) dγ(x)] dt,    0<t<H.

The positive and subtracted terms are integrable since ∫₀ᴴ p_s(v) ds≤C√H uniformly in v. Deterministic zero-age terms at time t do not affect this spacetime trace, by Fubini over the graph. Formulas for negative y follow by reflection.

## 1. Exact starting time

Define the right tail numbers

    ℓ = liminf_{x→+∞} f(x),
    S_K(z) = ∫_{x>K} exp(zx−z²f(x)/2) dγ(x),
    c = lim_{K→∞} liminf_{z→∞} [−(2/z²)log S_K(z)].

Then, for every y>0,

    inf supp μ_y = min(c, inf_{x>y} f(x)),                          (1)
    max(0,ℓ−1) ≤ c ≤ min(ℓ,H−1).                                   (2)

The K limit is increasing. Exactness gives S_K≤1. The lower bound ℓ−1 follows from the unrestricted Gaussian moment generating function and f≥ℓ−ε eventually. For every K, a fixed compact interval above K where f≤ℓ+ε gives the upper bound ℓ+ε on the corresponding exponential rate. The interval [z,z+1] contributes at least const·exp((1−H)z²/2) for z>K.

To prove (1), put q=z²/2. By R083 §3,

    L_y(q)=∫e^(−qt) dμ_y(t)=(2/z)∫_{x>y}e^(−qf(x))sinh(z(x−y)) dγ>0.

On each bounded interval (y,K), the exponential rate −lim (log integral)/q is inf_{y<x<K}f(x), by restriction to a compact subinterval and continuity, also when the infimum is approached at y. For x>K>y the hyperbolic sine is uniformly comparable to exp(z(x−y)) for large z. Rates of a sum of two positive terms combine by the minimum when one rate is a limit and the other a lower limit. Sending K→∞ proves that −limsup_q (log L_y(q))/q is the right side of (1). For a finite compactly supported signed measure this number is its first support time whenever its Laplace transform is positive on q>0. Indeed an exponential upper bound gives the lower support bound by R083 §3 (or by the Cauchy transform of its image under t↦exp(−t)). The opposite inequality follows directly from the total variation bound.

In particular all slices vanish on the open exterior y>max(0,sup{x:f(x)≤t}) for 0<t<c on the right; the left uses its own c. To see that this holds for individual slices, fix an exterior point and t<c. There is a neighborhood of y and a time margin for which the running infimum stays above t: on a bounded segment use compactness, on the tail use ℓ≥c. Fubini and (1) give zero on the open spacetime neighborhood. The heat kernels there have positive spatial separation from nearby-time sources and are continuous, including at age zero by extending the off-diagonal heat kernel by zero. This gives the assertion pointwise.

## 2. An exact formula using sublevel tail masses

For s>ℓ let ν_s=1_{f≤s}γ on the positive half axis and put

    β(s)=liminf_{r→∞} (−2/r²)log ν_s((r,∞)),     1/∞=0.

Then β(s)≥1, and

    c = inf_{s>ℓ} [s−1/β(s)].                                        (3)

To verify the tail calculation independently, any nonzero finite positive measure ν with unbounded right support dominated by γ has, for every fixed K>0,

    limsup_{z→∞}(2/z²)log∫_{x>K}e^(zx) dν = 1/β,
    β=liminf_{r→∞}(−2/r²)log ν((r,∞)).

For 0<b<β, integrate the tail bound exp(−br²/2) by parts to obtain an upper exponent 1/b. If β is finite take r_j with ν((r_j,∞))≥exp(−(β+ε)r_j²/2) and z_j=(β+ε)r_j, yielding lower exponent 1/(β+ε). For β infinite any fixed positive tail interval gives the zero lower exponent. Continuity and s>ℓ ensure that ν_s has unbounded right support.

Restricting S_K to {f≤s} now gives c≤s−1/β(s). Conversely, for δ>0 choose K so large that f≥ℓ−δ on x>K. Partition its range into finitely many slabs of width at most 2δ, with first upper endpoint ℓ+δ, all upper endpoints s_i>ℓ, last endpoint ≥H. Every slab contributes at most exp(−(s_i−2δ)z²/2)∫_{x>K} e^(zx)dν_{s_i}. The rate of their finite sum is at least inf_{s>ℓ}(s−1/β(s))−2δ. Let δ↓0. This proves (3), including plateaus at H and infima attained only at s↓ℓ.

## 3. A necessary oscillation on at least one tail

Assume that at each nonzero point x one of its two flanks obeys the R087 local condition

    liminf_{r↓0} log(1/r) sup_{0≤v≤r}|f(x±v)−f(x)|=0.              (4)

Locally Hölder functions off the origin satisfy (4), regardless of changes of their local constants toward infinity. For a nonconstant exact graph,

    min(c₊,c₋) < H−1.                                               (5)

Equivalently, on at least one wing there is s>ℓ with s−1/β(s)<H−1. This gives a quantitative necessary width of deep tail sublevels. In particular both tails cannot attain the maximal rate H−1, even if their pointwise values oscillate indefinitely.

Here are the full changes required in the earlier-band results. Suppose c₊=c₋=T=H−1>0 (Ef=1 and continuity give H>1 in the nonconstant case). By (2) the strict sublevels below T are bounded and the closed sublevels {f≤t}, t<T, compact. The exact positive-z transform forces at least one zero: if the compact part has no zero it decays exponentially in z², the favored tail S_K decays exponentially since c₊>0 and K can be chosen large, and the unfavored negative half contributes O(1/z). A nonzero zero is excluded by the R087 one-flank restriction of the exact real identity as in (4). Thus the origin is the unique zero.

Formula (1) supplies the entire exterior support hypothesis in R087. All remaining steps of its separated-component birth bound use only this exterior support, compactness, positivity of the Gaussian source measure, the semigroup identity and (4). The endpoints have no jumps and sweep all intermediate positions. Thus every sublevel below T is an interval. In the origin upper-envelope proof of R083 §2 the sole use of the far-tail gap is an o(1) bound for S_K(z) in the favored direction; it follows here from c₊>0 and on the reflected wing from c₋>0. All local monotonicity requirements in that proof concern sufficiently small values of f and hold by the interval property. Now the no-plateau, steepness, boundary-trace and first-exit uniqueness arguments in R083 §§3,5 and R084 apply on every compact time subinterval of (0,T). They use the exterior support and origin estimate just verified. Hence Q_t is positive and supported on the early interval [−A(t),B(t)] for 0<t<T.

Both limiting radii A,B at T are positive finite. Indeed if a whole wing had f≤T then its fixed-tail tilted mass would have exponent at most T−1 (restrict to [z,z+1] on the reflected or favored wing). This contradicts c=T. Thus each wing contains a point where f>T, which bounds its respective early radius. The R084 endpoint test (x+A(t))(x−B(t)) and dominated convergence now imply

    ∫_{outside [−A,B]}(f−1−T) dγ ≥ A [g(B)−B γ((B,∞))]
                                       + B [g(A)−A γ((A,∞))] > 0.

The integrand is ≤0 since 1+T=H. This proves (5). If connected early sublevels and a unique zero are assumed directly, the same proof dispenses with (4) entirely.

## 4. Example of a previously uncovered tail geometry

On either wing suppose f tends uniformly to its global supremum H outside disjoint intervals centered at the integers n≫1 of lengths ≤C exp(−a n²), with f≥d on those intervals. Here a>0 and H−1+1/(1+2a)≤d<H. No interior geometry is prescribed except (4). For d≤s<H, the γ-mass past r of these intervals is at most exp(−(1+2a)r²/2+O(r+1)), by the Gaussian density, their lengths, and summation from n≥r−1. Thus β(s)≥1+2a. For s≥H the complete tail has β(s)=1. Formula (3) gives c=H−1. If the hypothesis holds on both wings, (5) rules the nonconstant graph out. For instance H=1.6, a=1, d=0.95 satisfy the inequalities and allow the two tail lower limits 0.95, despite the all-wing effective rate being 0.6. An identical bound using a level-dependent β covers arbitrary placements and depths via (3).

This is an early-time support and rigidity result. The sign of an arbitrary individual slice after the two effective tail rates remains open.
