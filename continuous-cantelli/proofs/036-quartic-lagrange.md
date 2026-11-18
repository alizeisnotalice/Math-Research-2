# Independent quartic large-jet audit by Lagrange inversion

Write g(u)=p₁(x₀+u), b∼N(0,d), K_b(w) the conditional density of u+b u² when u has density g. Let N be even, C_N=4^{-N}binom(2N,N) and B_N=(4√d)^N C_N E|Z|^N. The real asymptotic coefficient A_N of the unrestricted quartic density satisfies
A_N/B_N = 2g(0) + g''(0)/(2d(2N−1)) + O_{x₀,d}(N^{-2}).

Indeed let t(w) be the solution at zero of t+b t²=w. Then K_b(w)=t'(w)[g(t(w))+g(−1/b−t(w))] near zero for b≠0. By the change of contour w=t(1+bt) in Cauchy's coefficient formula,
[w^N]K_b = [t^N] (g(t)+g(−1/b−t))(1+bt)^{−N−1}.
Divide by binom(2N,N)(−b)^N. The resulting finite sum is
Σ_{j=0}^N r_{N,j} (−b)^−j [g^{(j)}(0)+(−1)^j g^{(j)}(−1/b)]/j!,
with r_{N,j}=∏_{l=0}^{j−1}(N−l)/(2N−l) ≤2^−j. For |b|≥2 the source Taylor coefficients are bounded by H=exp(1/2)/√(2π) on all real unit disks, so this is uniformly 2g(0)+O(1/|b|). More precisely after averaging b and −b, its expansion through order b^−2 equals 2g(0)+r_{N,2}g''(0)b^−2+O(b^−4) uniformly in N≥4: the constant and linear-j parts of the expansion of g(−1/b) and −[g'(0)−g'(−1/b)]/(2b) cancel at order b^−2. The j=2 term supplies the asserted coefficient; b^−1 and b^−3 drop out by symmetry. Uniformity of the fourth-power remainder follows by bounded first four derivatives on |u|≤1 and the summed 2^-j bounds for j≥4.

The conditional expansion may be integrated to any finite order even though K_b has a moving square-root singularity. Fix cutoff A=2. On |b|≤A the analytic branch near w=0 has both root imaginary parts bounded uniformly (the far root is −1/b−t), so this contribution is analytic, with geometric coefficient bounds. On A<|b|≤|w|^−1/2 take Taylor series in q=4bw; on |q|≤3/4 the two roots and inverse discriminant stay uniformly bounded. Taylor's remainder in |q|≤1/2 is O_N(|bw|^{N+1}), integrable against the Gaussian law of b. On |b|>|w|^−1/2 the actual conditional density is at most 2‖g‖∞ 1_{1+4bw>0}(1+4bw)^−1/2. Extract exp[−1/(4d|w|)] and integrate the remaining Gaussian-weighted factor, bounded O(1+1/|w|) by splitting the discriminant at 1/2. The finite Taylor terms in this tail are also negligible to all algebraic orders.

Normalized signed averaging and integration use E|b|^{N−2}/E|b|^N=1/[d(N−1)] and E|b|^{N−4}/E|b|^N=1/[d²(N−1)(N−3)]. The contribution to these moments of |b|≤2 and its entire analytic contribution are geometrically bounded and negligible after normalization by B_N. Since r_{N,2}=(N−1)/(2(2N−1)), this proves the formula. Truncating input u at any fixed r>0 also changes the transform by a density analytic near zero and therefore preserves it. The leading factor 2g(0) comes from both inverse branches; using only the local inverse would halve the coefficient.
