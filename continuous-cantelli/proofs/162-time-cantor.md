# A time Cantor set of recoveries for one deterministic Gaussian graph

Let γ have standard Gaussian density g, write p_t for the centered Gaussian density of variance t and P_t for convolution with it, with P₀ the identity. Let E be the closed middle-thirds Cantor set affinely placed in [3/4,7/8]. There is a bounded nonnegative Borel function f on ℝ, satisfying

    ∫ exp(zx−z²f(x)/2) γ(dx) = 1                          (z real),

whose signed heat remainder

    Q_t = p_t dx − ∫_{f(x)≤t} P_{t−f(x)}δ_x γ(dx)

has the following properties.

* At every t∈E, its density q_t satisfies q_t≥p_t/2 on the whole line. Every complementary gap of E in [3/4,7/8] contains an open interval where Q_t has a strictly negative smooth spatial test. Consequently negative times accumulate at every point of E.
* The early time distribution on (0,1) has no atoms. On this time interval t↦Q_t is continuous in total variation, and each q_t has a continuous spatial representative. The construction gives a spatially smooth q_t at every t∈E and a jointly smooth density near each selected observation of a negative peak.
* There is a single compact spatial interval J=[2,2.1] and a constant c>0 such that the peak in each time gap has value at most −c somewhere in J. At every t∈E the densities are discontinuous in time in the norm L∞(J), and the points exhibiting the negative peaks can all be chosen to approach one spatial accumulation point as the gaps shrink.

These statements also hold with independently prescribed disjoint open pulse windows in [3/4,7/8], with the global positive bound outside the windows for 0<t<1 and a negative observation in each window. The accumulation assertion applies at limit points of distinct windows. The example is a measurable graph. Its positive-measure set of discontinuities is described below.

## Space and time budgets

Enumerate the gaps (a_j,b_j), j≥1, in any order. Choose s_j<v_j strictly inside that gap, for example at fractions 1/3 and 2/3 of its length, and put d_j=v_j−s_j. The early mark on the spatial patch K_j will lie in [s_j,s_j+δ_j] with 0<δ_j<d_j/4, and it will be tested for negativity at u_j=s_j+2δ_j. For 3/4≤h≤7/8 and y∈J, direct completion of the square gives

    sup_x p_{t−h}(x−y)/p_t(x) = √(t/(t−h)) exp(y²/(2h)),    t>h.

Thus for v_j≤t≤1 and s_j≤h≤s_j+δ_j this ratio is at most 27/√d_j. Also p_s(y)≤(1−D)g(y) on [3/4,7/8]×J, for a fixed D>0: the maximum ratio there is √(8/7)exp(−2/7)<1.

Take a fixed compact fat Cantor set K₀ in J, obtained by deleting the middle relative fraction 2^(−n−2) at stage n=1,2,… from each remaining interval. It has positive length and every cylinder has positive length of K₀. Follow the leftmost chain of cylinders to its endpoint x_*. At successively larger depths choose the right sibling of the chain as the jth spatial interval; its intersection K_j with K₀ is compact and has positive Gaussian mass m_j. The siblings are disjoint, with x_* as their only accumulation point. They can be made arbitrarily short at each choice. Arrange that

    (27/√d_j)m_j ≤ 2^(−j−2),       m_j(d_j/2)^(−(j+1)/2) ≤ 2^(−j).

The second constraint is used only for spatial derivatives at recovery times. The set K={x_*}∪⋃K_j is closed and nowhere dense. It inherits the following uniform lower bound on the length of its complement in each interval I, where Δ=diam K₀:

    |I\K| ≥ c₀ min(|I|, |I|²/Δ).

Indeed if the central half of I meets K₀, take a nested cylinder about such a point with length first below |I|/4, for |I|≤4Δ. Its length is comparable to |I| since the contraction ratios are bounded above and below. Stage-n lengths are comparable to Δ2^(−n) by positivity of the infinite retained product, so the next central gap has length at least c|I|²/Δ and lies inside I. If the central half misses K₀, or if |I|>4Δ, the linear bound follows directly.

## Atomless pulses and their sign

Let η_j=1_{K_j}γ. At a Lebesgue density point of each K_j, positive unit-integral smooth approximate identities give a nonnegative smooth test ψ_j supported in the interior of J, arbitrarily close to that point, with

    ∫ψ_j p_{s_j} − ∫ψ_j dη_j < −D min_J g/2.

Choose the supports also within distance 1/j of K_j. On the cylinder [l_j,r_j] of K_j set h(x)=s_j+δ_j(x−l_j)/(r_j−l_j) for x∈K_j. Its image law under η_j is atomless and has bounded density, with a bound depending on j. Shrink δ_j as needed. At u_j the contribution of this pulse to the test is ∫P_{u_j−h(x)}ψ_j(x)dη_j; the heat ages lie in [δ_j,2δ_j], and the heat averages tend uniformly to ψ_j as δ_j↓0. The source p_{u_j} tends to p_{s_j}. The signed test at u_j is consequently less than −D min_J g/3, since every previously activated pulse is subtracted. Set h(x_*) to any time in [3/4,7/8]. This singleton has zero mass.

For t<1, use only these early graph marks in q_t. At a time outside all gaps, every pulse that has contributed already satisfies t≥v_j. The kernel-ratio bound and the summable budgets yield q_t≥(1−Σ 2^(−j−2))p_t≥p_t/2. The same reasoning applies on (v_j,b_j) after the jth pulse. These estimates include all accumulation times and all gap endpoints in E. Exclusion of t=1 here allows for immediate stops in the restart used below.

## Density regularity

Let η=Ση_j. As h_*η has no atoms, Q_t has an L¹ density for every 0<t<1. For each source x with h(x)≠t the past-kernel term varies continuously in L¹ at t, by strong continuity of the positive-age Gaussian kernel or because the mark remains in the future. Its L¹ norm is at most one. Integration against η yields total-variation continuity.

At t∈E, each active pulse has t−h≥d_j/2. For any spatial derivative of order k its contribution has uniform norm at most C_k m_j(d_j/2)^(−(k+1)/2), by Gaussian scaling. This is summable: the terms with j≥k obey the second budget since d_j/2<1. Termwise differentiation thus gives a C∞ density in the spatial variable at each such t. This argument does not assert joint smoothness across a time in E.

In the interior of one gap, pulses in all gaps to its left have already ended at or before its left endpoint, and those on the right are still in the future. On any compact subinterval of this gap the older kernels therefore have a common strictly positive age. All their fixed space and time derivatives are uniformly bounded after integration by constants times Σm_j. The current pulse also has positive age separated from zero in a neighborhood of u_j, since u_j=s_j+2δ_j. This proves joint smoothness there and persistence of the negative test on an open interval. During the ramp the same pulse gives a continuous density: its contribution with 0<t−h<ε is bounded in uniform norm by C_j∫₀^ε a^(−1/2)da using the bounded density of h_*η_j; use a continuous cutoff in age for the older, jointly continuous part, then let ε↓0. Times strictly outside [3/4,7/8] have either no active sources or a positive common minimum age on compact time neighborhoods. This proves the continuity claims at every early time.

Since ∫ψ_j=1 and supp ψ_j⊂J, some point there has q_{u_j}≤−D min_J g/3. It can be taken inside supp ψ_j, hence it tends to x_* along any sequence of distinct gaps. For every t∈E there are distinct gaps shrinking to t. At their observation times q_{u_j} has the uniform negative value, whereas q_t≥p_t/2 on J. Continuity of each spatial density makes this a lower bound for the essential supremum norm of the difference as well.

## Completing the bounded deterministic graph

At time 1 just before any restart stops the remaining source and target measures are

    μ=γ−∫_K P_{1−h(y)}δ_y γ(dy),      λ=1_{Kᶜ}γ.

The pulse bound gives μ≥γ/2. They have the same mass and first moment, and μ precedes λ in convex order because λ−μ is a sum of mean-preserving Gaussian spreads minus their starting masses. The potential difference is the strictly positive function

    U(x)=∫ |x−v|(λ−μ)(dv)=∫_K∫₀^{1−h(y)} p_s(x−y)ds γ(dy).

Normalize the two laws and use their regular Root embedding, as in round 154. The barrier r is lower semicontinuous and its stop is uniformly integrable. It is finite and continuous locally on Kᶜ, where the target density is locally bounded below. Its terminal law gives no mass to K, so the stopping time equals r at the terminal position almost surely. These are the general-start Root existence and local continuity results used in round 154. [Cox–Obłój–Touzi](https://arxiv.org/abs/1505.03169) · [Bayraktar–Bernhardt](https://arxiv.org/abs/2010.14695).

For clarity, the bound on r only needs the vacancy estimate and the uniform upper heat age 1−h≤1/4. Gaussian decay gives U(x)≤C exp(−dist(x,J)²/(2·(1/4))), with harmless changes to the constant, by integrating the s^(−1/2) prefactor. This is o(g(x)/(1+x²)) on both tails. In unnormalized units U is the total expected occupation density by Tanaka's formula. For |x| large set l=1/(1+|x|). Each of the intervals [x−l,x] and [x,x+l] contains a point with r≤l²: otherwise the initial density μ≥γ/2 dominates a constant times g(x) times the first Dirichlet sine on that interval. Its killed evolution up to l² would give U≥c g(x)l² on the middle third, contradicting the tail bound there. Choose two such points on opposite far tails to trap a fixed central compact region by time 1. All surviving components of {r>1} now have uniformly bounded length and a uniform ratio of maximum to minimum g within each: the central ones lie in a fixed compact interval, and any far point x still active lies between the two comparison points within l of it.

A bounded component I of {r>t} with width D contains precisely λ(I) surviving mass at time t. No path could stop there before t; from t on, paths cannot cross the two barrier endpoints, and terminal endpoint masses are zero. Within I every subinterval of length D/2 has λ-mass at least c min(1,D/Δ)λ(I) by the density ratio and vacancy estimate. To stay active for h more time within I has probability at most C exp(−c'h/D²), uniformly over starting points; this follows by iterating the static interval survival bound over periods D². Consequently after C' D²(1+log_+(Δ/D)) each surviving subcomponent has width at most D/2. Iterate; this cost is increasing in D and the sum over dyadic widths is finite uniformly for bounded starting width. No open surviving component remains after that common finite time. This proves bounded r.

Set f=h on K and f=1+r on Kᶜ. The two branches are carried by mutually singular parts of γ, giving a bounded deterministic Borel graph. Optional sampling of the bounded-time exponential Brownian martingale for the residual transport gives ∫exp(zx−z²r(x)/2)λ(dx)=∫exp(zx)μ(dx). Multiplication by exp(−z²/2) and addition of the early term yield the stated exact transform by direct Gaussian convolution. The assigned values on K are ≤7/8, while its dense open complement has assigned values ≥1 and f is continuous there. Thus its discontinuity set is exactly K. No assertion about time atoms in the later Root branch is needed.
