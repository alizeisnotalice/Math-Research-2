# Noncancellation for infinitely many arbitrarily placed separated near-maxima

Let v=M−f≥0 be bounded and continuous and suppose the exact Gaussian-output identity gives, for all real t,
∫ exp(itx−t²v(x)/2)γ(dx)=exp(−Mt²/2), M>0.
The following tail theorem allows the infimum of v to be zero without attainment.

**General criterion.** Let measurable disjoint input sets I_j be contained in [x_j−W,x_j+W], W<∞, with centers separated by |x_j−x_k|≥d>0. Let v≥d₀>0 outside their union. Assume inf_j ess inf_{I_j}v=0, essential with respect to γ. Define A_j(t)=∫_{I_j}e^{−t²v/2}dγ, Z(t)=(∑ A_j(t)²)^{1/2}, p_j=A_j²/Z², and the probability µ_{j,t} on I_j proportional to e^{−t²v/2}dγ, omitting zero-mass sets. If
t²∑ p_j(t) ∫_{I_j}(x−x_j)² dµ_{j,t}(x)→0,                       (1)
then exactness is impossible. More generally it suffices that the unscaled second moment in (1) tend to zero and the ℓ² norm of the complex normalized vector in (2) below have positive lower limit.

**Distribution of pairs and stable frequency intervals.** Give (j,x,y) probability proportional to 1_{I_j}(x)1_{I_j}(y)e^{−t²Y}γ(dx)γ(dy), where Y=(v(x)+v(y))/2. Denote this law Q_t. Its normalizer is Z². Conditional on j, x and y are independent with common law µ_{j,t}; the index law is p_j. Boundedness of v permits differentiation, giving
d E_{Q_t}Y/dt=−2t Var_{Q_t}(Y), hence ∫₁∞ t Var_{Q_t}(Y)dt<∞.
For every fixed length L>0 there exist disjoint intervals [n,n+L], n→∞, on which ∫ t²Var_{Q_t}(Y)dt→0. Otherwise on all intervals in the partition [jL,(j+1)L] the integrals would eventually be bounded below, contradicting integrability with weight 1/t.

Let
b_j(t)=(A_j/Z) E_{µ_{j,t}} exp[it(x−x_j)].                          (2)
With ̅v=E_Q Y,
b'_j=(A_j/Z)E_{µ_j}[e^{it(x−x_j)}(i(x−x_j)−t(v−̅v))].
Jensen and variance decomposition bound ∑|b'_j|² by 2∑p_j E_{µ_j}(x−x_j)²+4t²Var_Q(Y). Indeed ∑p_j E_{µ_j}(v−̅v)² is at most twice the pair variance because conditional pair fluctuations have half the single-input variance. Thus ∫_{n}^{n+L}‖b'‖₂ dt→0 on the chosen intervals by Cauchy–Schwarz. Also ‖b(t)−(A_j/Z)_j‖₂²≤t²∑p_j E_{µ_j}(x−x_j)², so ‖b(n)‖₂→1 under (1). On any finite t interval the vector is continuously differentiable in ℓ² by boundedness of v, x−x_j, and termwise domination by the square-summable vector of prior masses times a constant (Z is bounded below there).

**Elementary separated-frequency bounds.** For finite coefficients a_j and any interval J of a sufficiently large fixed length L=L(d), there are c,C>0 depending on L,d alone such that
c‖a‖₂ ≤ ‖∑a_j e^{itx_j}‖_{L²(J)} ≤ C‖a‖₂.
To verify this, choose a half-length T=L/2 and integrate on the centered interval against the triangle (1−|t|/T)_+. Its transform is T sinc²(Tu/2), so its diagonal contribution is T∑|a_j|² and the off-diagonal row sum is ≤4π²/(3T d²) by ordering the separated points to either side of each fixed point. Translation changes only phases. If T is large the weighted lower bound gives the claimed lower bound. For the upper bound, twice the wider triangle of half-length 2T dominates the interval indicator; apply the same row estimate. Completion extends the bounds to ℓ² coefficients.

On a stable interval J=[n,n+L] first freeze b=b(n). The lower inequality gives an L² norm bounded away from zero. For the remaining error write b_j(t)−b_j(n)=∫_n^t b'_j(s)ds and use Minkowski and the upper inequality on [s,n+L] (bounded above by the full J). The error is ≤C∫_J‖b'(s)‖₂ds=o(1). Start with finite truncations, then pass to the limit using the same bounds on ℓ² tails and dominated convergence for ‖b'‖₂. For each t the actual series for the transform is also absolutely convergent since ∑|b_j(t)|≤1/Z(t), locally uniformly in t. Thus the normalized core transform ∑ b_j(t)e^{itx_j} has L² norm bounded below on J.

But Z(t)≥c_η e^{−ηt²/2} for every η>0, using a single positive-mass subset of one I_j with v≤η. The omitted input contributes at most exp(−d₀t²/2) in absolute value. Exactness bounds the normalized core transform above by [exp(−Mt²/2)+exp(−d₀t²/2)]/Z(t)→0 uniformly on these intervals. This contradicts the lower bound. No common phases of the x_j, spacing growth, or uniform lower mass per near-maximum was used.

**Escaping power-shaped wells.** Let I_j=(x_j−w_j,x_j+w_j), 0<w_j≤W, with the above fixed center separation and disjointness. Suppose v(x_j)=ε_j>0 with inf ε_j=0, v≥d₀>0 off all wells, and throughout each well
c₁(|u|/w_j)^α ≤ v(x_j+u)−ε_j ≤ c₂(|u|/w_j)^α,   0<α<2,
with constants independent of j. This configuration is excluded, regardless of how slowly the ε_j decrease or how rapidly the widths w_j shrink.

Fix A>0 and first take |x_j|≤At. Put l=w_j t^{−2/α}. The integral defining A_j is at least const·φ(x_j)exp(−t²ε_j/2)l by |u|≤l, since At l→0 uniformly. On l≤|u|≤w_j, the Gaussian source tilt |x_j u| is ≤(c₁/4)t²(|u|/w_j)^α for all large t: the maximum ratio is bounded by A W t^{1−2/α} if α≥1 and A W/t if α≤1. Thus an integrable exp[−c(|u|/l)^α] bounds all conditional scaled tails after division by this mass, and E_{µ_j}u²≤C l². It follows that t² E u²→0 uniformly for these centers. For |x_j|>At, ∑ A_j≤P(|X|>At−W) and Z(t)=exp[−o(t²)] by the lower test above and Z≤1, so the pair-index mass of these wells is exponentially small. Since |u|≤W, (1) holds. A finite collection of exceptional wells where v has strictly positive minimum may be absorbed into the complement by decreasing d₀. The general continuous case, including arbitrarily crowded peaks and very flat contacts, remains open.
