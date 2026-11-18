# Connected early sublevels and a global excess budget

Write f=1+Ea−a, so the exact equation is E exp(zX−z²f(X)/2)=1 and Ef=1. Suppose f is bounded continuous nonnegative, its unique zero is the origin, T>0, and both far tails have liminf f≥1+T. Assume that for each 0<t<T the sublevel set {f≤t} is an interval containing 0. No monotonicity or comparability is assumed on the part above level T. Then for 0<t<T the signed heat remainder Q_t is the subprobability law of Brownian motion from zero still inside the expanding early interval. If −A and B are the limits of its two endpoints as t↑T (A,B>0), a necessary condition is

    ∫_{x<−A or x>B} (f(x)−1−T) γ(dx)
           ≥ A [φ(B)−B barΦ(B)] + B [φ(A)−A barΦ(A)] >0.           (1)

In particular excesses above 1+T anywhere in the two remaining tails must pay the displayed gap. Negative deviations there count against the budget. If f≤1+T everywhere on the two tails, there can be no such graph.

**Reduction.** The capped profile min(f,T) is nondecreasing as one moves outward on either half-line: every strict sublevel below T is an interval containing zero. On either wing it obeys f(±r)~r²/[2log(1/r)] at the origin, by the restriction and upper-envelope contradiction in Section 2 of the preceding radial proof. For the favored tilt the opposite half-axis contributes at most O(1/|z|); small positive radii use only monotonicity below T. On each fixed compact interval away from the origin f has a positive minimum, and at infinity the tail lower bound strictly exceeds 1. Thus that proof applies on each wing separately.

For t<T and y>B(t)=max{f≤t}, the signed resolvent formula (4) of the radial proof uses only x>y. Its absolute value decays faster than exp(−qt) by some exponential margin, since on any bounded [y,R] one has min f>t, and on the far tail f≥1+t+δ for some δ>0. The same argument gives pointwise zero remainder throughout the open exterior in time and on the left. The boundary no-plateau and steepness proof uses only this vanishing, the strict positivity of the one-sided resolvent and monotonicity near the early edge. The ordinary Brownian first exit has the same past potential as the graph measure on the exterior. Integrability of the two traces at zero, the positivity argument for the symmetrized spacetime heat kernel, and uniqueness on each [0,T'] for T'<T are as in the preceding proof. Distinct left and right radii require no symmetry: the comparison of the two past potentials uses only that each past spatial mean lies between the two present boundary positions. The good-set diagonal estimates work within each wing, and the opposite wings are separated after excluding a neighborhood of time 0.

**Endpoint moment.** Both radii increase to finite positive limits since T is below both far-tail lower limits. The exact Gaussian-transform identity gives the mass, first and second moments of Q_t by testing the future sources with 1, x, x²−(f−t), respectively. Since (y+A(t))(y−B(t))≤0 on its support and Q_t≥0, these moments give

    ∫_{x<−A(t) or x>B(t)} [(x+A(t))(x−B(t))−(f(x)−t)]γ(dx)≤0.

Let t↑T. Dominated convergence applies since f is bounded and the endpoints converge. Integration by parts gives ∫_{x>B}(x+A)(x−B)γ(dx)=barΦ(B)+A[φ(B)−B barΦ(B)], and on the left gives the reflected expression with A and B exchanged. This proves (1); both bracketed quantities are positive Gaussian call values. No positivity of Q at times above T was used.

The same argument gives a family of global sign estimates: for every polynomial p nonnegative on [−A(t),B(t)],

    ∫_{f>t} [exp(−(f(x)−t)∂²/2)p](x) γ(dx) ≥0.

The exponential is a finite sum for polynomials, so no backward heat evolution has to be defined on arbitrary functions. The formula follows by comparing derivatives of the exact exponential transforms, and has the same conclusion at T by taking the limit with polynomials nonnegative throughout the limiting interval. It includes nonnegative localized polynomial tests of any degree. It does not imply positivity of Q on a general disconnected past set.
