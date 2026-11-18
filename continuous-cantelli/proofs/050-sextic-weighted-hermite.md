# Sextic convergent jets with arbitrary real flat remainder

Suppose v≥0 is continuous and bounded, the Gaussian-source output is Gaussian, and near x₀ one has v(x₀+u)=V(u)+O_J(|u|^J) for every fixed J, where V(u)=cu⁶H(u)², c>0, H real analytic at zero and H(0)=1. This is impossible. No differentiability of the actual remainder is needed.

## Coefficient identity

Write g(u)=p₁(x₀+u). Fix a small core |u|<r. Its output has an analytic density A(w) near w=0, w=y−x₀, since the complement exterior positive Gaussian kernels are holomorphic and dominated by C/|u| on |w|<r/4, while zero-variance exterior output cannot reach this interval. The real asymptotic coefficients satisfy
[w^j] A(w)=∫_{|u|<r}g(u) [w^j]p_{v(x₀+u)}(w−u)du
+ [w^j] ∑_{k=0}^{⌊j/4⌋} (2^k k!)^−1 ∂_w^{2k}[g(w)V(w)^k].       (1)
The first integral uses pointwise off-diagonal kernel coefficients for u≠0. To prove it for the analytic reference alone fix a finite degree j. On |u|>8|w|, Cauchy Taylor remainder on the input disk of radius |u|/4 is bounded by C_j |w|^{j+1}|u|^−(j+4)exp(−c'/u⁴), integrable. On |u|≤8|w| and |u−w|≥|w|/2, the unexpanded kernel contributes smaller than all powers by v≤C w⁶ and the Gaussian distance estimate. On |u−w|<|w|/2 put u=w+|w|³t. Then the integrand in dt is g(u)p_{V(u)/w⁶}(t). Taylor expand in u about w to any fixed order L≥j: on this local segment V(u)/w⁶ is bounded above and below, and its qth u-derivative is O_q(|w|^{-q}). Gaussian derivatives in the variance give a fixed polynomial times an integrable Gaussian. The integrated remainder is O_L(|w|^{2L+2}). The Taylor terms may be extended to all t at superpolynomial error; their even Gaussian moments are (2k−1)!!(V(u)/w⁶)^k, odd moments zero. This gives exactly the finite diagonal sum shown; terms k with 4k>j start too late. All fixed off-diagonal Taylor coefficients over the omitted inner range are negligible by the Cauchy estimates. To replace the reference by the actual graph in the first integral, apply the same exterior estimates and bound the inner difference of the two unexpanded Gaussian kernels by C|v−V||u|^−9 through the heat equation, uniformly on the real axis. Flatness makes the integral smaller than every fixed power. Comparing asymptotic coefficients along the full-measure equality set proves (1).

## A Gaussian multiplier and positivity of every off-diagonal kernel

Take N→∞ along multiples of four, a=(cN)^−1/4, W=w/a, z=u/a. Fix
t₀=3/8, z₀=t₀^{1/4}, λ=(9t₀−2)/(4t₀^{3/2}).
Inspect P_N=[W^N] e^{λNW²/2} A(aW). On |W|=1/√λ the physical argument tends to 0, giving |P_N|≤B exp(N f₀), f₀=(1+logλ)/2.

Write S=Nv(x₀+az)/a²=z⁶h(az), where h(u)=v(x₀+u)/(c u⁶)→1. An individual off-diagonal weighted coefficient is a positive prefactor times
[W^N] exp[N z W/S + (N/2)(λ−1/S) W²].
For λS≥1 this is nonnegative at every real z since N is even. Otherwise its sign is the sign of He_N(√N z/√[S(1−λS)]). All roots of He_N lie in (−2√N,2√N), by the tridiagonal recurrence with off-diagonal entries √j and the maximum absolute row sum bound. We have
S(1−λS)/z² = h(az)^{1/3} S^{2/3}(1−λS)
≤ h(az)^{1/3} (3/5)(2/(5λ))^{2/3} < 1/4
on the core after shrinking r. Indeed λ>(2/5)(12/5)^{3/2}, as follows by squaring the positive constants. Thus all off-diagonal contributions to P_N are nonnegative.

More strongly their integral is at least exp[(f₀+δ)N] for a fixed δ>0 and large N. For real z in a closed positive interval close to z₀ put p=z/S, d=1/S−λ and R=2/(p+√(p²−4d)). These quantities and d are strictly positive there; the discriminant has a uniform positive gap. The small positive saddle R solves pR−dR²=1 and 1−dR²>0. The uniform one-variable saddle estimate gives
[W^N]exp[N(pW−dW²/2)]
= exp[N(pR−dR²/2−log R)] /√[2π N(1−dR²)] ·(1+o(1)).
Indeed on its circular Cauchy contour the real exponent relative to its maximum at angle 0 is −(1−cos θ)[pR−dR²(1+cos θ)], with a uniform quadratic gap away from zero; rescaling θ by √N gives the positive Gaussian factor. The actual exponent including e^{−Nz²/(2S)} is
F=−(R−z)²/(2S)+λR²/2−log R.
At z=z₀ and S=z₀⁶, R=2z₀/3 is the smaller root (1−dR²>0 from t₀>1/3), and F−f₀=(1/2)log(27/11)−4/9>0. Positivity of this last elementary number follows e.g. from log x=2∑_{l≥0} ((x−1)/(x+1))^{2l+1}/(2l+1), x=27/11, retaining the first three positive terms. By continuity F≥f₀+2δ on a fixed interval of z around z₀ for large N, uniformly for the actual S since h(az)→1 there. Kernel prefactors and the source Jacobian du=a dz cost only polynomial factors. This proves the stated lower estimate.

## The diagonal for a pure sextic jet

For H=1 the following short argument completes the proof without a second saddle. Write b∼N(0,c), so E b^{2k}=(2k)!c^k/(2^k k!). After scaling and weighting (1), group diagonal coefficients by k. For each k the term from the constant Taylor coefficient g₀ is positive, proportional to
D_{N,k}=binom(6k,2k)(2k−1)!! N^−k · (λN/2)^{(N−4k)/2}/[((N−4k)/2)!].
A higher even source coefficient g_{2j} divided by that same positive factor has modulus bounded by |g_{2j}|(9a²/(4λ))^j: use binom(6k+2j,2k)/binom(6k,2k)≤(3/2)^{2j} and the falling-factorial bound for the Gaussian multiplier's coefficient ratios. Odd source coefficients make no contribution. Since g is analytic, the sum for j≥1 is O(a²), uniformly over k. Each diagonal group is positive for large N, and hence the off-diagonal lower bound contradicts |P_N|≤B e^{N f₀}.

## Any convergent sextic jet

For general H the kth diagonal group, up to the positive factor (2k−1)!! N^{−k}, is the two-variable coefficient
[W^N T^{2k}] e^{λNW²/2} (W+T)^{6k} g(a(W+T)) H(a(W+T))^{2k}.
This follows directly by replacing ∂^{2k} by (2k)![T^{2k}] at W+T in (1). The support restricts 0≤k≤N/4. If H,g are replaced by 1, the full kth group is D_{N,k}. For any k its absolute value is ≤C D_{N,k} exp(CNa) by the previous binomial-ratio argument applied to the absolute Taylor majorants g^+(q)H^+(q)^{2k} at q=3a/(2√λ). Uniform Stirling gives N^−1 log D_{N,k}=f(k/N)+O(log(N+2)/N), where
f(ν)=ν(6log6−2log2−4log4)+ν(log(2ν)−1)+(1/2−2ν)[1+log(λ/(1−4ν))].
The continuous endpoint values are f(0)=f₀ and f(1/4)=¼(6log6−2log2−4log4+log(1/2)−1)<f₀ (direct substitution). Therefore the total absolute contribution for k/N in sufficiently small fixed neighborhoods of 0 and 1/4 is ≤exp[(f₀+δ/2)N] for large N.

For k/N=ν in the remaining compact subinterval of (0,1/4), every full kth group is positive for all sufficiently large N, uniformly in ν. We give the details of this real-saddle claim. At a=0 the Cauchy torus with positive radii R=√[(1−4ν)/λ], U=R/2 has precisely two maximal phase points: (W,T)=(R,U) and (−R,−U). Equality in |exp(λNW²/2)| forces the W angle to be 0 or π and equality in |W+T|≤R+U forces T to have the same angle. In logarithmic radii x,y the phase is
Φ₀(x,y)=λ e^{2x}/2+6ν log(e^x+e^y)−x−2ν y.
At (log R,log U) its gradient vanishes. Its Hessian is positive definite, the sum of diag(2λR²,0) and 6ν RU/(R+U)² times [[1,−1],[−1,1]], uniformly on the compact ν-range.

For each sign σ=±1 add 2ν log H(σa(e^x+e^y)) to Φ₀; the real branch of the logarithm exists for small a. The perturbation and all needed fixed derivatives are O(a) on fixed neighborhoods. By the implicit function theorem the two real log-radius saddles exist O(a) from the original, preserving positive-definite Hessians. In two disjoint small angular patches of the torus about (0,0) and (π,π), deform both log radii to their respective stationary real values, using a smooth cutoff equal to one near the saddle and zero outside the patch. This torus deformation is legal: the original bivariate Cauchy integrand is holomorphic for nonzero W,T on the swept compact region (use its original product, so W+T=0 creates no pole), and the Cauchy two-form is closed. The transition zones and the region off the patches are exponentially smaller with a uniform real-action gap from the unperturbed saddles: the normalized unperturbed modulus has a strict maximum only at the two designated points, compactness gives a gap, and the radii and the H factor change its logarithmic rate by only O(a). At possible zeros of W+T the same upper bound follows directly by continuity of normalized modulus.

On the constant-radius part of each patch, expand the holomorphic logarithm of the phase in the two real angles. The first derivatives vanish by construction, the real part has negative Hessian equal to minus the positive log-radius Hessian, and the third derivatives are uniformly bounded. Two-dimensional Gaussian domination after rescaling angles by √N gives its contribution
exp[N Φ_σ(x_σ,y_σ)] {g(σa(e^{x_σ}+e^{y_σ}))/(2πN√det Φ_σ''(x_σ,y_σ))} (1+o(1)),
uniformly in ν. The sign phases at the negative point are +1 since N,2k,6k are even. Both displayed leading numbers are real positive. Each remainder is relative to its own saddle and the rest of the contour has a fixed extra gap. Thus their sum is positive, as is the original real coefficient.

Combining the nonnegative middle diagonal, the endpoint error and the strictly larger off-diagonal lower bound contradicts the analytic upper bound. This covers general convergent sextic jets with completely arbitrary real all-order-flat value error, and does not require angular analyticity of that error.
