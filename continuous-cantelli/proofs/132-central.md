# Central positivity at large Laplace parameter

Write γ for the standard Gaussian law with density g. Let (X,F) be any joint law with X standard normal and F nonnegative. The signed heat remainder subtracts the heat flow started at the marked times from the free heat kernel. Its Laplace transform at q>0 has density

    R(q,y)=K(q,y)−E[e^(−qF) K(q,y−X)],    K(q,x)=exp(−z|x|)/z,   z=√(2q).

This follows by Laplace transformation of each positive heat flow separately. For n≥0 define P₀=1 and Pₙ₊₁(a)=(2n+1+a)Pₙ(a)−aPₙ′(a). Differentiation gives

    (−∂q)ⁿ K(q,x)=exp(−a)Pₙ(a)/z^(2n+1),    a=z|x|.

Each Pₙ is monic of degree n with strictly positive coefficients, since each coefficient of degree i≤n acquires the positive multiplier 2n+1−i in addition to its shifted predecessor. Let

    Cₙ=2g(0) Σ_{j=0}ⁿ binom(n,j) Mⱼ ∫₀^∞ exp(−u)Pₙ₋ⱼ(u)du,
    M₀=1,    Mⱼ=(2j/e)^j for j≥1.

Leibniz's rule for the marked term produces factors (z²F)^j exp(−z²F/2)≤Mⱼ after multiplication by z^(2n+1). This pointwise estimate permits arbitrary dependence of F on X. Integrate the remaining positive kernels in X using g≤g(0) and substitute u=z|y−x| on each half-line. This bounds the normalized marked term by Cₙ/z at every y. All differentiations under the expectation follow from the same bounds on a compact neighborhood of q>0, or from uniform boundedness there of the exponential-polynomial factors. Consequently

    z^(2n+1)(−∂q)ⁿ R(q,y) ≥ exp(−z|y|)Pₙ(z|y|)−Cₙ/z.

This global lower estimate gives a strictly positive central region at every fixed derivative order. More explicitly choose Bₙ with exp(Bₙ)>Cₙ. For all sufficiently large z the region includes z|y|≤log z+n log log z−Bₙ. Indeed on a bounded a-interval the free profile has a positive minimum. For a≥max(1,n), use Pₙ(a)≥aⁿ and monotonic decrease of exp(−a)aⁿ. At the proposed boundary z exp(−a)aⁿ tends to exp(Bₙ). This conclusion also holds for n=0. The exact exponential identities, continuity and an upper bound on F are unnecessary for this estimate. The strip shrinks in physical space as q grows, and leaves the intermediate spatial cuts open.
