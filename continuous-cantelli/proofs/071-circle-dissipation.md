# Low-variance tilts and local phase cancellation on the circle

Use Haar probability dx/(2π) on the circle of length 2π and subtract min v from a continuous stationary graph. Thus v≥0, its zero set S is compact and nonempty, and for every nonzero integer j one has ∫ exp(ijx−j²v/2)dx=0. Write Z_t=∫exp(−t²v/2)dx with consistently normalized dx and let μ_t be the corresponding probability tilt. Every weak limit of μ_t for t→∞ is supported on S: for any ε>0 the set v<ε/2 has positive measure and dominates the weight of v≥ε exponentially. In particular M(t)=∫v dμ_t→0.

## Integrable local Fourier errors

For every pair of fixed integers m≥1, k, the sequence

    b_{m,k}(n)=∫ exp[i(mn+k)x] exp[−(m²−1)n²v(x)/2] dμ_n(x)

satisfies ∑_{n≥1}|b_{m,k}(n)|²/n<∞. Indeed |b|≤|∫ exp[i(mn+k)x]dμ_{mn}(x)| since Z_{mn}≤Z_n. The same unweighted harmonic integrates to zero under μ_{mn+k} when mn+k>0. The density derivative ∂_t μ_t=−t(v−M(t)) μ_t has L¹ norm at most √I(t), I(t)=t²Var_{μ_t}(v). The exact dissipation law M'(t)=−I(t)/t gives ∫_T^∞ I(t)dt/t=M(T). By Cauchy–Schwarz between mn and mn+k,

    |b_{m,k}(n)|² ≤ |k| ∫_[min(mn,mn+k),max(mn,mn+k)] I(t)dt.

For n≥N large enough that mN>2|k|+2 and n≥1, each window has 1/n≤2m/t there; any t belongs to at most |k|/m+2 windows. Thus the weighted tail sum is bounded by 2m|k|(|k|/m+2)M(mN−|k|). The case k=0 is exactly zero. The finitely many low frequencies are harmless.

There is a single exceptional set E⊂N, with ∑_{n∈E}1/n<∞, such that b_{m,k}(n)→0 off E simultaneously for every fixed m,k. For each j take N_j increasing so that beyond N_j the union of exceptions |b_{m,k}|>1/j for the first j pairs in an enumeration costs at most 2^(−j) in harmonic mass. Take the union of those exceptional tails. Equivalently a growing finite collection of nearby harmonics is controlled at each nonexceptional n.

## Conditional harmonics at the minimum set

Lift μ_n to the joint probability of (x,θ,s)=(x,nx mod 2π,n²v(x)/2), compactifying s in [0,+∞]. Every joint weak limit Π along n→∞ outside E has its first marginal supported on S, and

    E_Π[ exp(imθ) exp(−(m²−1)s) | x ] = 0      for every m≥1.      (1)

At m=1 take the second factor identically one, including at s=∞. For m≥2 take its value at infinity to be zero. Indeed testing against exp(ikx) gives exactly b_{m,k}(n), and these are continuous tests on the compact product. Fourier uniqueness on the first circle and disintegration give (1). In particular the conditional first harmonic vanishes at almost every macroscopic location. If the mark s is finite and is a measurable function of x under Π, all conditional harmonics of θ vanish, so its conditional law is Haar.

For arbitrary S, take 0<c<π and the centered angular arc A(α,c)={d(θ,α)≤c}. For Π-almost every x, every α and c the conditional law gives mass at most 1/(1+cos c) to A(α,c), since the real part of the centered phase is ≥cos c on the arc and ≥−1 elsewhere. This is nontrivial for c<π/2 and sharp for the abstract first-moment condition. Moreover ∫d(θ,α)² dΠ(θ|x)≥π²/4, using t²≥π²/4−π cos t on [0,π] and the same conditional moment. The scalar inequality follows by differentiating; sin(t)/t decreases on (0,π), so the minimum of t²+π cos t is at π/2.

For a fixed finite set A={a₁,…,a_J} containing S, this gives a scale restriction independent of the shape of the wells:

    liminf_{n→∞, n∉E} n² ∫ d(x,A)² dμ_n ≥ π²/4.                  (2)

In fact pass along any subsequence attaining finite liminf, extract the joint limit and limits of all na_j modulo 2π, and use disjoint small fixed neighborhoods of the a_j, whose boundaries have zero limiting mass. On each such neighborhood of a_j the nearest center is a_j and n d(x,a_j)≥d(nx,na_j). Sum the limits of the smaller bounded functions over the neighborhoods and apply (1). If a limiting atom belongs only to a_j, its conditional second moment enters weighted by that atom. Likewise the upper limit of μ_n({n d(x,A)≤c}) is at most 1/(1+cos c), c<π/2. If S is infinite, the same argument applies locally at each isolated point a of S: along a subsequence with limiting mass p at a, within a fixed small neighborhood containing no other point of S, the limiting upper mass of d(x,a)≤c/n is at most p/(1+cos c). These bounds need no simultaneous Diophantine returns and no quantitative lower bound on the volume of sublevel sets.

## Quantitative localization by bandwidth

With m=1, over a dyadic frequency block N≤n≤2N, 1≤H≤N/4, the previous dissipation bound yields

    ∑_{N≤n≤2N} (1/n) ∑_{|k|≤H}|∫e^{i(n+k)x}dμ_n|²
          ≤ C H³ M(N−H).                                          (3)

Let F_H be the Fejér kernel of integral one (with Fourier support |k|≤H), and let λ_n=e^{inx}μ_n be the modulated complex measure. Parseval yields ∥F_H*λ_n∥²₂ bounded by the inner sum in (3); hence the same holds for the squared L¹ norm. One can choose H=H(N)→∞ sufficiently slowly so that H³M(N−H)→0 on all dyadic blocks, implying small spatially localized modulation on most integers of each block. A fixed H also has a finite total weighted error across all n by the tail estimate above. Thus localization can be taken on genuinely shrinking spatial windows 1/H(N), with its guaranteed rate governed by the tilt energy. This controls absolute localized mass; normalization by arbitrarily small window masses would require an additional lower bound.

The estimates concern necessary spectral and spatial constraints. They do not classify arbitrary continuous stationary profiles or exclude an arbitrary compact null minimum set.

## Centered polynomial test and the full-band bound

For integer H≤n/8 let

    A_H(n) = sup_{deg h≤H, ||h||∞≤1} |∫ e^{inx}h(x)dμ_n(x)|,

where the supremum ranges over complex trigonometric polynomials. Put λ=n²/2, τ=nH+H²/2, L(λ)=log Z_{√(2λ)} and
D=exp([L(λ+τ)+L(λ−τ)]/2−L(λ))≥1. Then

    A_H(n) ≤ min[1,D−1+sqrt(D²−1)].                                (4)

Write B=∂_x²/2+i n ∂_x on the space of degree-H trigonometric polynomials. Bernstein's trigonometric inequality, iterated once, bounds its sup-norm operator norm by τ. On exp(ikx), B has eigenvalue −k²/2−nk. Hence stationarity shows for every real number a that

    ∫ e^{inx} [exp((v(x)−a)B)h](x)dμ_n(x)=0.

Indeed on each harmonic this is the exact (n+k)-th tilted identity times the fixed scalar exp[a(k²/2+nk)]/Z_n; k=−n is excluded. The exponential is defined on a finite-dimensional space for any sign of v−a. Its difference from the identity is bounded in sup norm by exp(τ|v−a|)−1 for each x. Choose a=(log E_mu exp(τv)−log E_mu exp(−τv))/(2τ) so E_mu exp(±τ(v−a))=D. With X=τ(v−a), Cauchy–Schwarz gives E|sinh X|≤sqrt(E(cosh X−1) E(cosh X+1))=sqrt(D²−1). This proves (4), since exp(|X|)−1=cosh X−1+|sinh X|. For H=0 use the original exact identity.

In particular A_H(n)²≤C log D, by splitting log D≤1 and log D>1. Since L''(λ)=Var_{μ_{√(2λ)}}(v), its centered difference satisfies

    log D = (1/2)∫_{−τ}^τ (τ−|s|)Var_{μ_{√(2λ+2s)}}(v)ds
          ≤ C H ∫_{n−2H}^{n+2H} I(t)dt.

Use t=√(2λ+2s) and τ=O(nH); t stays between n−2H and n+2H. Therefore, with a universal constant and 1≤H≤N/8,

    ∑_{N≤n≤2N} A_H(n)²/n ≤ C H² ∫_{N−2H}^{2N+2H} I(t)dt/t.        (5)

This strengthens (3) for L¹ localization. In fact convolution against the positive normalized Fejér kernel maps every bounded test function of sup norm one to a degree-H polynomial of sup norm at most one, so ∥F_H*(exp(in·)μ_n)∥₁≤A_H(n) by duality. A degree-H bounded test may vary with n, the minimum locations, or the data. For fixed H the total harmonic sum is finite. One may choose a slowly growing H(N) such that the right side tends to zero; even exceptional frequencies of finite harmonic sum suffice if the widths are increased sufficiently slowly by diagonal selection over all n.

If Z_t≥c t^(−P) for large t, then for every n and H=o(n/log n) one has A_H(n)→0. In particular this lower bound on Z holds when v has a positive-power one-sided upper bound at just one minimum. To see the estimate, monotonicity of M gives M(t)≤[log Z_{t/2}−log Z_t]/(3t²/8)=O(log t/t²). In the proof of (4) choose a=0 instead. Then A_H(n)≤Z_{√(n²−2τ)}/Z_n−1 ≤ exp[τ M(√(n²−2τ))]−1=o(1). More generally, if Z_{t/2}/Z_t is bounded for all large t, H=o(n) suffices. Thus these hypotheses give spatially localized cancellation for every frequency with resolution approaching the single wavelength. No regularity at other minima is used.

## A geometric exclusion for arbitrarily many minimum clusters

Suppose Z_n≥c n^(−P) and n d(X,S)→0 in μ_n-probability. The latter holds if v≥c₁d(·,S)^α, α<2: the tilt mass where d≥η/n is at most c^(−1)n^P exp(−c₁η^α n^(2−α)/2) for every fixed η>0. If for infinitely many n one can partition S into finitely many nonempty clusters (the number may grow without bound) of covering radii ≤b/n with fixed b<π/2, and pairwise distances ≥G_n where nG_n/log n→∞, stationarity is impossible.

In fact choose H=o(n/log n) with HG_n→∞. On the G_n/3 neighborhood of each cluster set h₀=exp(−ina_j) where a_j is a covering center, and assign zero elsewhere. These neighborhoods are disjoint. The normalized Fejér convolution h=F_H*h₀ has sup norm ≤1 and degree H. On each G_n/6 neighborhood it differs from exp(−ina_j) by at most 2∫_{|t|>G_n/6}F_H(t)dt/(2π)=O(1/(HG_n)), independently of the cluster count (the standard pointwise Fejér estimate is C/(H t²)). Since n d(X,S)→0 in probability, the phase e^{inX} differs from e^{ina_j} by an angular deviation ≤b+o_P(1) in the nearest cluster. Thus Re∫ e^{inx}h dμ_n≥cos b−o(1), contradicting (4) with the polynomial partition-function lower bound. For example the perfect uncountable set of all ∑_{j≥1} ε_j 2^(−2^j), ε_j∈{0,1}, has clusters of radius O(ℓ_j²) separated by at least ℓ_j/2 after j digits, where ℓ_j=2^(−2^j), so n_j near ℓ_j^(−3/2) works. The phases of different clusters can be arbitrary.

There is also a fixed separation ratio version. Suppose c₁d(·,S)^α≤v≤C₁d(·,S)^α globally for 0<α<2. Then Z_{n/2}/Z_n≤D=(4C₁/c₁)^(1/α) uniformly for every closed S. Indeed |S_{qr}|≤q|S_r| for q≥1, where S_r is the distance-r neighborhood: the complementary gaps of lengths l_j contribute min(l_j,2r) each to |S_r|, in addition to |S|. Use the exponential layer-cake formula ∫exp(−t² c d(x,S)^α/2)dx=∫₀^∞ exp(−s)|S_{(2s/(ct²))^(1/α)}|ds. Also Z_n has the required polynomial lower bound at any minimum. Fix b<π/2. There is a finite B depending only on b,α,C₁/c₁ such that stationarity fails whenever for arbitrarily large n the cluster radii are ≤b/n and the pairwise gaps are ≥B/n. To see this, use H=floor(δ n), δ≤1/8. For t between √(n²−2τ) and n one has M(t)≤8 log D/(3t²), so A_H(n)≤exp(Cδ log D)−1 by (4)'s uncentered proof. First choose δ small to make that bound <cos(b)/4, then choose B large to make the Fejér error C/(δB)<cos(b)/4. Taking the same n-scale concentration argument with fixed arbitrarily small η gives the contradiction. In particular this excludes two-branch Cantor minimum sets with a fixed sufficiently small contraction ratio: their level-j cells have radii ≤Cρ^j and gaps ≥c(1−2ρ)ρ^(j−1), allowing integer n satisfying both fixed inequalities at every sufficiently large j. The Cantor set may be embedded in a proper arc to avoid the identified endpoints.

## The natural measure test for equally scaled Cantor minima

Let S be the two-branch, equal-ratio Cantor set on an arc of length l≤2π, keeping the two end intervals of relative length 0<ρ<1/2 at every stage. Endpoints of the full circle may be identified. Let ν be the natural probability giving each level-j cell mass 2^(−j). If v=c d(x,S)^α with 0<α<2, c>0, then, without assuming stationarity,

    ∫e^{inx} dμ_n − ∫e^{inx}dν → 0.                                (6)

In particular stationarity requires the natural measure to have Fourier coefficients tending to zero. To prove this set σ=n^(−2/α), d=log(2)/log(1/ρ). There is a constant c₀>0 with Z_n≥c₀σ^(1−d): take a cell generation at which lρ^j is between ρσ and σ. Each of its cells lies in S_σ, their interiors are disjoint, their count is 2^j comparable to σ^(−d), and the integrand is at least exp(−c/2) there. Instead choose a generation j=j(n) with σ≪lρ^j≪1/n; this exists since α<2 and ρ is fixed. Within each cell the distance to S is exactly the distance to its own scaled copy of the Cantor set, since its two endpoints belong to S. Thus every cell has identical tilt mass. The complement consists of at most 2^j gaps, each with integral ≤2∫₀^∞exp(−c(u/σ)^α/2)du=O(σ). Normalized discarded mass is O(2^j σ^d)=o(1). The phase e^{inx} varies by o(1) across each chosen cell. Integrating against μ_n and ν now gives (6).

For the standard middle-thirds Cantor set scaled onto the full circle, its natural measure is the law of 2π∑_{j≥1}2B_j/3^j for independent uniform zero-one digits B_j. For n=3^h, terms with j≤h contribute unit phase and those with j>h reproduce the same characteristic function as at n=1. The modulus equals ∏_{j≥1}|cos(2π/3^j)|>0, since no factor vanishes and the sum of (1−|cos(2π/3^j)|) converges. This gives a continuous profile with a perfect uncountable minimum set excluded directly by (6). The claim concerns this class of graph shapes; no classification of all variance profiles with Cantor minima is asserted.

The same Rajchman necessity holds if v/d(·,S)^α has a positive continuous limit on S, in the precise sense v(x)=[a(x)+o(1)]d(x,S)^α uniformly as d(x,S)→0, with a continuous positive function a on the circle. Take the same intermediate cells. Their common constant-amplitude mass is W_n(z)=∫_cell exp(−n² z d(x,S)^α/2)dx, identical across all cells for any z>0. Fix 0<m≤a≤M, and include 1 in this fixed parameter range if necessary. For z within an enlarged compact positive range,

    0≤−∂_z log W_n(z)≤ (2/z) log[W_n(z/2)/W_n(z)]
                          ≤ (2 log 2)/(α z).

The first inequality bound uses convexity of log W_n. For the second apply the layer-cake proof of neighborhood doubling inside the cell; its endpoints belong to S and all internal gaps have the same truncated-length formula. Thus ω_n(z)=W_n(z)/W_n(1) is uniformly bounded above and below and uniformly Lipschitz on the fixed compact range, as is 1/ω_n.

Choose a representative x_C∈S in each intermediate cell C. Uniform continuity of a and the uniform relative remainder show that the actual cell mass is [1+o(1)]W_n(a(x_C)), uniformly in C: bracket the variance there between (a(x_C)±ε_n)d^α, ε_n→0, and use the log derivative bound. The complement again contributes o(1) to the normalized tilt mass, since a global two-sided comparability c₁d^α≤v≤C₁d^α follows by continuity and positivity off S. Let h_n(x)=1/ω_n(a(x)); this is a uniformly bounded equicontinuous family. For a stationary v, ∫e^{inx}h_n(x)dμ_n→0 along all n: polynomial partition mass makes every fixed nearby Fourier error tend to zero, and Arzelà–Ascoli with uniform trigonometric approximation applies to h_n. On the other hand, summing the cells using nρ^j→0 gives

    ∫e^{inx}h_n dμ_n
      = [2^j W_n(1)/Z_n] ( ∫e^{inx}dν + o(1) ),

where all measures are taken with consistent normalization; the factor is uniformly bounded above and below, by the cell mass comparisons and the negligible complement. Hence ν̂(n)→0 as asserted. In particular positivity and continuity of the coefficient trace suffice at the whole middle-thirds minimum set, with no derivative bounds on the relative remainder.
