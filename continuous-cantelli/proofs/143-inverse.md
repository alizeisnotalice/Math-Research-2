# Resolving individual remote wells at the natural Gaussian width

Let the graph agree with the canonical core up to b at height T>0 (round 137). On any Borel source group A⊂{x≥b+d₀}, d₀>0, put d=x−b, v=f−T, m=d/v. Points with v=0 have zero Gaussian mass by round 138. Define the positive center measure ν by pushing √(2π/v) exp(d²/(2v))γ restricted to A forward under x↦m. Keep d as an auxiliary mark before this pushforward. The contribution at parameter t>0 is

    H_A(t)=∫_A exp(td−vt²/2)dγ=∫ p_{m/d}(t−m) dν(m,d).

The full edge output is asymptotic to J/t² and bounds this contribution. The center measure is finite: integrating in t>0 captures at least half of every Gaussian's mass, and the full edge output is integrable there. For R large, every interval of length √R in [R/2,4R] has center mass ≤C R^(−3/2). Indeed for those centers the Gaussian standard deviation is at most 2√(R/d₀); a fixed multiple enlargement of the interval captures at least half their mass, and throughout the enlargement the output is O(R⁻²).

It follows uniformly for intervals I⊂[R,2R] that

    |ν(I)−∫_I H_A(t)dt| ≤ C R^(−3/2).

Inside [R/2,4R], group centers by distance to the interval endpoints in steps of √R. Their probability of crossing an endpoint decreases as C exp(−c k²) in the step index k, and the preceding local mass bound sums these errors. Outside [R/2,4R], the probability of entering I is O(exp(−cR)) since m/d≤m/d₀. This compares sums over arbitrarily many wells at any mesoscopic width much larger than √R.

**Sharp inverse at a single location.** Suppose in addition d stays in a compact subset of (0,∞) on A and d→D>0 uniformly as m→∞ on A (essential uniformity suffices). Let μ_R be the measure R^(3/2)ν on the coordinate s=(m−R)/√R. For any j≥0,

    t² H_A(t)→j    if and only if    μ_R → j ds vaguely as R→∞.

The unit-interval bounds give local vague compactness and translation-bounded subsequential limits. At t=R+√R w the rescaled output R²H_A(t) is the convolution against μ_R of the variable-variance Gaussian with variance m/(Rd) on the scaled coordinate. On a compact center range this tends uniformly to p_{1/D}(w−s), locally uniformly in w. The truncated integral errors tend uniformly to zero: on R/2≤m≤4R use the unit-interval mass bound and a summable Gaussian majorant; outside, the kernel √R p_{m/d}(t−m) is exponentially small uniformly even after multiplication by R^(3/2). For the latter bound m is bounded below because d≥d₀ and f is bounded; also d has a finite upper bound. Consequently every joint subsequential limit of the local output is p_{1/D}*μ. If the output converges to j the convolution equals the constant j. The Gaussian is injective on translation-bounded measures: Fourier transform the difference from j ds as a tempered distribution and divide by the nonvanishing Gaussian multiplier on each compactly supported smooth frequency test. Thus μ=j ds. The converse follows by the same convolution limit at w=0.

If j>0, sampling x from A under its normalized contribution exp(td−vt²/2)γ/H_A(t) yields (m−t)/√t converging to N(0,1/D), by multiplying the center limit by the Gaussian transfer kernel. In particular tv→D in probability. No asymptotic expansion of an individual valley profile is assumed. This result conditions on an observed limiting fraction of edge flux; the existence of that fraction is not automatic.
