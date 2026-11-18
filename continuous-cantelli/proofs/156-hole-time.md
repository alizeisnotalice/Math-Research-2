# Quantitative restrictions on filling the holes of the exceptional set

Consider the completion used in round 154 from μ to λ=1_{Kᶜ}γ at time 1. If a Brownian embedding of these unnormalized marginals has duration bounded by H, then on every interval I of length D where μ≥c_I dx,
[
 λ(I) ≥ (2/π)c_I D\exp[-π² H/(2D²)].
]
Indeed start the killed Brownian flow inside I from the subdensity c_I sin(π(x−inf I)/D). Every path staying in I up to H must stop there. Its total probability is exactly the displayed quantity by the first Dirichlet eigenfunction. Extend the Brownian sample after stopping for this comparison. Thus even a shrinking sequence of holes cannot approximate a whole early plateau on a fixed interval while retaining uniform duration and a fixed positive starting density. If the early sink fully occupies I and hence λ(I)=0, no such bounded completion exists.

For the explicit construction μ≥γ/2, giving directly
[
 |I∖K|/D ≥ [\min_I g/(π\max_I g)]\exp[-π²H/(2D²)] .
]
For fixed I the minimum required H therefore tends to infinity at least logarithmically in the inverse vacant fraction. This obstruction concerns any Brownian completion from the specified positive restart law. Earlier signed evolution has already been included in that law.

**A sufficient condition on a compact part.** Suppose every interval I of length D on the chosen compact region and every subinterval J⊂I of length D/2 have λ(J)≥θ(D)λ(I) with 0<θ(D)≤1. Once a surviving Root component has width ≤D₀ and stays in this compact, its lifetime is bounded above by
[
 C\sum_{n≥0}D_n²[1+\log(1/θ_n)],\quad D_n=2^{-n}D₀,qquad
 θ_n=\inf_{D_n/2≤D≤D_n}θ(D),
]
provided this sum is finite. Indeed at stage n each component already of width ≤D_n/2 needs no further time. For every larger one the static exit inequality of round 154 shows that after the n-th displayed time increment its surviving subcomponents all have width ≤D_n/2. A point remaining active after infinitely many halvings would have an open active neighborhood, impossible. The same far-tail initial localization as in round 154 bounds the initial bounded components outside this compact, where the full Gaussian target is present.

For example on compacts it is enough to have the uniform vacancy fraction |I∖K|/D≥c exp(−C D^(−p)) at sufficiently small D, for any fixed 0<p<2, even if the target density vanishes on a positive-measure set. Indeed density comparability and the condition at D/2 give the required θ. The series still converges under exp(−C/[D²(log(1/D))^(1+δ)]) for δ>0. In the other direction, holes whose vacancy fraction has unbounded D²log(D/|I∖K|) along shrinking intervals in a fixed compact cannot have bounded completion time. The criterion at the remaining borderline is not identified here.
