# Convex Bayes risks of the real-tilt experiment

Let `γ=N(0,1)`. Assume that `f:ℝ→[0,M]` is measurable and satisfies
\[
\int \exp(tx-t^2f(x)/2)\,d\gamma(x)=1 \qquad (t\in\mathbb R).
\]
Write `P_t(dx)=exp(tx−t²f(x)/2)γ(dx)` and `q(t)=E_{P_t}f`. Assume `f` is nonconstant modulo `γ`. The experiment with `f=1` is denoted by `G_t=N(t,1)`.

**Tilt information.** Boundedness of `f` and Gaussian tails permit two differentiations of normalization on every compact set of parameters. Consequently `Eγ f=1` and
\[
E_{P_t}(X-tf)=0,\qquad E_{P_t}(X-tf)^2=q(t).
\]
In particular `q(t)` is the Fisher information in one observation from `P_t`. The first identity also gives
\[
D(P_t\Vert N(tq(t),1))=\frac{t^2}{2}q(t)[1-q(t)].
\]
Indeed subtract the logarithmic density of `N(tq(t),1)` relative to `γ` from the logarithmic density of `P_t`, and use `E_{P_t}X=tq(t)`. Every `P_t` is equivalent to `γ`, so `q(t)>0`. If the displayed entropy vanished for `t≠0`, equality of the two densities would make `f` affine almost everywhere. Boundedness would make it constant. Thus `q(0)=1` and `0<q(t)<1` for `t≠0`, as also shown in round 85.

## Gaussian priors

Fix `r>0` and `m∈ℝ`. Take `Θ∼N(m,1/r)` and observe `X|Θ=t∼P_t`. Let `Q` be the marginal law of `X` and `q̄=E_Q f`. Tonelli gives `q̄=E q(Θ)<1`. Completing the square in the joint density gives
\[
\Theta\mid X=x \sim
N\!\left(\frac{rm+x}{r+f(x)},\,\frac1{r+f(x)}\right).
\tag{1}
\]
It also gives
\[
\frac{dQ}{d\gamma}(x)=
\sqrt{\frac r{r+f(x)}}\exp\!\left(-\frac{rm^2}{2}+
          \frac{(rm+x)^2}{2(r+f(x))}\right).
\]
Thus `Q` is equivalent to `γ`. The all-real identity normalizes the joint law.

**Loss comparison.** Let `ℓ:ℝ→ℝ` be finite convex, coercive, and of at most polynomial growth, and use the translation-invariant loss `ℓ(d−Θ)` for decision `d`. The minimum Bayes risk is strictly higher for `P` than for `G` under every Gaussian prior above. This includes asymmetric absolute loss and even powers. It also holds for an even finite nonconstant convex loss whenever its Gaussian expectation is finite at scale `1/√r`, with no further growth requirement.

For the first claim write `L(s)=inf_a Eℓ(a+sZ)` for a standard normal `Z`. The function `(a,s)↦Eℓ(a+sZ)` is jointly convex and continuous by polynomial growth. It is coercive in `a` because Jensen gives `Eℓ(a+sZ)≥ℓ(a)`. On compact scale sets its minimum is therefore attained in a fixed compact set, by comparing to `a=0`. Convexity is preserved under partial minimization, so `L` is convex, and symmetry of `Z` makes it even. It is strictly increasing for `s≥0`. For `s_2>s_1`, choose a minimizer at `s_2`, represent `s_2Z` as `s_1Z_1+√(s_2²−s_1²)Z_2` with independent normals, and apply strict Jensen conditionally on `Z_1`. Strictness holds since `ℓ` is convex and nonaffine and the added Gaussian has full support. For an even loss covered by the second claim the optimal shift is zero, `L(s)=Eℓ(sZ)`, and the same convexity and strict increase follow using integrability up to `1/√r`.

Now `J(v)=L((r+v)^(−1/2))` is convex and strictly decreasing for `v∈[0,M]`, by convexity of `(r+v)^(−1/2)` and monotonicity and convexity of `L`. The conditional minimizer for a polynomial-growth loss is measurable. Convolution with a full-support normal makes the shift objective strictly convex for each `s>0`. Equality at two distinct shifts would make `ℓ` affine on every translate of the interval between them, by continuity and full support. The unique minimizer is continuous in `s>0` by continuity of the objective and uniform compactness of minimizers on compact scale sets. Thus (1) and Jensen give the minimum Bayes risk for `P` as
\[
E_Q J(f) \ \ge\ J(\bar q)\ >\ J(1),
\tag{2}
\]
where `J(1)` is the minimum Bayes risk for `G`.

**Squared loss and information.** For squared loss, `MMSE_P=E_Q (r+f)^(−1)`. Expanding the reciprocal exactly around `q̄` gives
\[
\mathrm{MMSE}_P-\frac1{r+1}=
 \frac{1-\bar q}{(r+\bar q)(r+1)}
 +
 E_Q\frac{(f-\bar q)^2}{(r+f)(r+\bar q)^2}.
\tag{3}
\]
The mutual information of the joint experiment is explicit:
\[
I_P(\Theta;X)=\frac12 E_Q\log(1+f/r)
    <\frac12\log(1+1/r)=I_G(\Theta;X).
\tag{4}
\]
To see this, integrate the relative entropy of the Gaussian posterior in (1) to the Gaussian prior. The variance and squared mean-shift terms cancel by the law of total variance and `E_Q E[Θ|X]=m`. For a quantitative inequality, the second derivative of `log(r+v)` is at most `−(r+M)^(−2)` on `[0,M]`. Taylor's inequality at `q̄` and concavity between `q̄` and 1 give
\[
I_G-I_P\ \ge\
  \frac{1-\bar q}{2(r+1)}
  +\frac{\operatorname{Var}_Q(f)}{4(r+M)^2}.
\tag{5}
\]

## General regular prior and squared loss

Suppose the prior on `Θ` has a positive locally absolutely continuous probability density `π` with finite Fisher information `I_π=∫(π'/π)^2π`. Its minimum Bayes risk satisfies
\[
\inf_d E(d(X)-\Theta)^2\ >
  \frac{1}{I_\pi+E_\pi q(\Theta)}\ >
  \frac{1}{I_\pi+1}.
\tag{6}
\]
This specializes the van Trees bound. The joint score in the prior coordinate is `S=π'(Θ)/π(Θ)+X−Θ f(X)`. Its second moment is `I_π+E_π q(Θ)` by the first two tilt identities. For an estimator of finite mean squared error, multiply `d(X)−Θ` by a smooth cutoff `χ(Θ/R)` equal to 1 around zero. Integration by parts in `Θ` against the joint density gives the expected cutoff minus the expectation of `(d(X)−Θ)χ'(Θ/R)/R`. The second term tends to zero by Cauchy–Schwarz and the bounded derivative of the cutoff. The score pairing converges by its own Cauchy–Schwarz bound. Hence `E[(d(X)−Θ)S]=1`, and Cauchy–Schwarz gives the first bound with a weak inequality.

Equality would give `d(x)−t=c[π'(t)/π(t)+x−t f(x)]` almost everywhere in the joint law, with `c=(I_π+E_π q(Θ))^(−1)>0`. The joint law is equivalent to the product of the prior and `γ`. By Fubini, subtract the equality at two good values `x_1,x_2`. A constant then equals `−ct[f(x_1)−f(x_2)]` for almost every `t`, so `f(x_1)=f(x_2)`, a contradiction. The strict inequality holds for the minimum risk. If an estimator has finite mean squared error, its conditional squared error is integrable almost everywhere, so the conditional laws have finite second moments and the conditional mean attains the minimum. If every estimator has infinite risk, (6) follows immediately. Strictness of the second inequality follows as the prior gives zero mass to `t=0`.

For a general prior, the rightmost expression in (6) is also a lower bound for the Gaussian shift risk. The exact comparison in (2) uses a Gaussian prior, for which this Gaussian lower bound is attained. Equation (3) gives the strict excess above the stronger van Trees bound for the Gaussian prior.

## Provenance and scope

The pointwise information law `q(t)<1` is already in [round 85](085-tilt-sign.md), which also proves that the two experiments are incomparable under fixed Markov transformations. The risk and information comparisons here follow from the exact Gaussian posterior and that pointwise information law. The regular-prior estimate (6) uses the classical Bayesian information bound. See [Gill and Levit, Applications of the van Trees inequality](https://dspace.library.uu.nl/bitstream/handle/1874/1190/van_trees.pdf?isAllowed=y&sequence=1), §1. All arguments apply to bounded measurable graphs up to null sets, including the normalized counterexamples of [Kleptsyn and Kurtzmann, Theorem 1.1 and Proposition 2.1](https://arxiv.org/pdf/1202.2250). The paragraph after Theorem 1.1 formulates the continuous problem.

Local symbolic verification used SymPy 1.14.0 to check the posterior square completion and the reciprocal identity in (3). The risk and information comparisons follow from the analytic proof above.
