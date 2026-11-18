# Convex Bayes risks of the real-tilt experiment

Let `γ=N(0,1)`. Assume that `f:ℝ→[0,M]` is measurable and satisfies
[
int exp(tx-t^2f(x)/2),dgamma(x)=1 qquad (tinmathbb R).
]
Write `P_t(dx)=exp(tx−t²f(x)/2)γ(dx)` and `q(t)=E_{P_t}f`. Assume `f` is nonconstant modulo `γ`. The experiment with `f=1` is denoted by `G_t=N(t,1)`.

**Tilt information.** Boundedness of `f` and Gaussian tails permit two differentiations of normalization on every compact set of parameters. Consequently `Eγ f=1` and
[
E_{P_t}(X-tf)=0,qquad E_{P_t}(X-tf)^2=q(t).
]
In particular `q(t)` is the Fisher information in one observation from `P_t`. The first identity also gives
[
D(P_tVert N(tq(t),1))=
rac{t^2}{2}q(t)[1-q(t)].
]
Indeed subtract the logarithmic density of `N(tq(t),1)` relative to `γ` from the logarithmic density of `P_t`, and use `E_{P_t}X=tq(t)`. Every `P_t` is equivalent to `γ`, so `q(t)>0`. If the displayed entropy vanished for `t≠0`, equality of the two densities would make `f` affine almost everywhere. Boundedness would make it constant. Thus `q(0)=1` and `0<q(t)<1` for `t≠0`, as also shown in round 85.

## Gaussian priors

Fix `r>0` and `m∈ℝ`. Take `Θ∼N(m,1/r)` and observe `X|Θ=t∼P_t`. Let `Q` be the marginal law of `X` and `q̄=E_Q f`. Tonelli gives `q̄=E q(Θ)<1`. Completing the square in the joint density gives
[
Thetamid X=x sim
N!left(
rac{rm+x}{r+f(x)},,
rac1{r+f(x)}
ight).
	ag{1}
]
It also gives
[

rac{dQ}{dgamma}(x)=
sqrt{
rac r{r+f(x)}}exp!left(-
rac{rm^2}{2}+

rac{(rm+x)^2}{2(r+f(x))}
ight).
]
Thus `Q` is equivalent to `γ`. Formula (1) uses the all-real identity to obtain the correctly normalized joint law.

**Loss comparison.** Let `ℓ:ℝ→ℝ` be finite convex, coercive, and of at most polynomial growth, and use the translation-invariant loss `ℓ(d−Θ)` for decision `d`. The minimum Bayes risk is strictly higher for `P` than for `G` under every Gaussian prior above. This includes the usual asymmetric absolute loss as well as even powers. It also holds for an even finite nonconstant convex loss whenever its Gaussian expectation is finite at scale `1/√r`, with no further growth requirement.

For the first claim write `L(s)=inf_a Eℓ(a+sZ)` for a standard normal `Z`. The function `(a,s)↦Eℓ(a+sZ)` is jointly convex and continuous by polynomial growth. It is coercive in `a` because Jensen gives `Eℓ(a+sZ)≥ℓ(a)`. On compact scale sets its minimum is therefore attained in a fixed compact set, by comparing to `a=0`. One may compute the conditional infimum over rational shifts and obtain a measurable estimator arbitrarily close to it by taking the first such shift in an enumeration. Thus the minimum risk value is the expectation of the conditional infimum. Convexity is preserved under partial minimization, so `L` is convex, and symmetry of `Z` makes it even. It is strictly increasing for `s≥0`: for `s_2>s_1` choose a minimizer at `s_2`, represent `s_2Z` as `s_1Z_1+√(s_2²−s_1²)Z_2` with independent normals, and apply strict Jensen conditionally on `Z_1`. Strictness holds since `ℓ` is convex and nonaffine and the added Gaussian has full support. For an even loss covered by the second claim the optimal shift is zero, `L(s)=Eℓ(sZ)`, and the same convexity and strict increase follow using integrability up to `1/√r`.

Now `J(v)=L((r+v)^(−1/2))` is convex and strictly decreasing for `v∈[0,M]`, by convexity of `(r+v)^(−1/2)` and monotonicity and convexity of `L`. Equation (1) and Jensen imply that the Bayes risk for `P` is
[
E_Q J(f)  ge J(ar q) > J(1),
	ag{2}
]
where `J(1)` is exactly the Bayes risk for `G`. For the polynomial-growth loss the pointwise minimizer is measurable as well: strict Jensen using the full-support normal makes the objective strictly convex in the shift for every `s>0`. Indeed if equality held at two distinct shifts, `ℓ` would be affine on every translate of the interval between them, by continuity and full support. The unique minimizer is continuous in `s>0` by continuity and uniform compactness of minimizers on compact scale sets.

**Squared loss and information.** For squared loss, `MMSE_P=E_Q (r+f)^(−1)`. Expanding the reciprocal exactly around `q̄` gives the identity
[
mathrm{MMSE}_P-
rac1{r+1}=

rac{1-ar q}{(r+ar q)(r+1)}
 +
 E_Q
rac{(f-ar q)^2}{(r+f)(r+ar q)^2}.
	ag{3}
]
The mutual information of the joint experiment is also explicit:
[
I_P(Theta;X)=
rac12 E_Qlog(1+f/r)
    <
rac12log(1+1/r)=I_G(Theta;X).
	ag{4}
]
To see this, integrate the relative entropy of the Gaussian posterior in (1) to the Gaussian prior. The variance and squared mean-shift terms cancel after taking expectation by the law of total variance and `E_Q E[Theta|X]=m`. This establishes (4) even if one does not use differential entropy notation. More quantitatively, the second derivative of `log(r+v)` is at most `−(r+M)^(−2)` on `[0,M]`. Applying one Taylor inequality at `q̄` and concavity between `q̄` and 1 gives
[
I_G-I_P ge
rac{1-ar q}{2(r+1)}
  +
rac{operatorname{Var}_Q(f)}{4(r+M)^2}.
	ag{5}
]

## General regular prior and squared loss

Suppose the prior on `Θ` has a positive locally absolutely continuous probability density `π` with finite Fisher information `I_pi=int(pi'/pi)^2pi`. Every estimator `d(X)` satisfies
[
E(d(X)-Theta)^2 ge
rac{1}{I_pi+E_pi q(Theta)} >

rac{1}{I_pi+1}.
	ag{6}
]
This is the van Trees bound specialized to the exact tilt experiment. To justify it on the whole real line, the joint score in the prior coordinate is `S=pi'(Theta)/pi(Theta)+X-Theta f(X)`. Its second moment is `I_pi+E_pi q(Theta)` by the first two tilt identities. For an estimator of finite mean squared error, multiply `d(X)−Θ` by a smooth cutoff `χ(Θ/R)` equal to 1 around zero. Integration by parts in `Θ` against the joint density gives the expected cutoff minus the expectation of `(d(X)−Θ)χ'(Θ/R)/R`. The second term tends to zero by Cauchy–Schwarz and the bounded derivative of the cutoff. The score pairing converges by its own Cauchy–Schwarz bound. The resulting identity is `E[(d(X)−Θ)S]=1`, so a final application of Cauchy–Schwarz proves the first inequality. Strictness of the second follows as the prior gives zero mass to `t=0`. For a general prior, the rightmost expression is itself a lower bound for the Gaussian shift risk. The exact comparison in (2) uses the Gaussian prior. Equation (3) also gives the strict excess above the refined van Trees bound for that prior.

## Provenance and scope

The pointwise information law `q(t)<1` is already in [round 85](085-tilt-sign.md), which also proves that the two experiments are incomparable under fixed Markov transformations. The Gaussian-prior risk and information comparisons use the exact posterior plus that information law. The regular-prior estimate (6) is an application of the classical Bayesian information bound. See [Gill and Levit, Applications of the van Trees inequality](https://dspace.library.uu.nl/bitstream/handle/1874/1190/van_trees.pdf?isAllowed=y&sequence=1), §1. All arguments apply to bounded measurable graphs up to null sets, including the normalized counterexamples obtained in [Kleptsyn and Kurtzmann, Theorem 1.1 and Proposition 2.1](https://arxiv.org/pdf/1202.2250). Continuity of the coefficient is the unresolved condition formulated in their paragraph after Theorem 1.1. These risk comparisons impose no additional regularity.

Local symbolic check used SymPy 1.14.0 to reduce the difference between the exponents in the posterior square completion and the remainder after subtracting the three terms of the reciprocal Taylor identity to identically zero. The risk and information comparisons are analytic.
