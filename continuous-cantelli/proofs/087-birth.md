# An exponential mass cost for creating a separated component of a sublevel set

Let f=1+Ea−a be a continuous exact normalized graph and assume both tails have liminf f≥1+T for some T>0. The set {f=0} is nonempty: otherwise f has a strictly positive minimum on every compact set, and splitting the real transform between a large compact interval and tails with f≥1+T/2 would give a zero limit. All sublevels below T are compact. Put C_t=conv{x:f(x)≤t}=[L(t),R(t)]. The intervals are nested, continuous from the right (decreasing intersections of compact sets commute with extreme endpoints here), but at this stage their outward jumps are allowed. The set at t=0 need not yet be connected.

**Outer heat support.** The one-sided resolvent identity in the single-well proof gives supp Q_t ⊂ C_t for 0<t<T, provided C_t contains 0. More generally one may always use conv(C_t ∪ {0}) for this statement. Indeed for y>max(0,R(t)), choose a compact tail segment [y,M] on which min f>t. Beyond it the Gaussian integral costs exp(q) and the far-tail lower bound exceeds 1+t by a fixed margin. The resulting exponential time decay of the signed resolvent forces its time support at y to begin after t. Spatial separation from all sources of age zero gives pointwise vanishing throughout the open exterior as in the preceding proof. On the left reflect the argument. If f vanishes only at zero the two versions coincide.

**Birth bound.** Suppose R jumps outward at t₀∈(0,T), from b=R(t₀−) to x₀=R(t₀)>max(b,0), and put D=x₀−max(b,0). Then for each fixed 0<r<D and all sufficiently small s>0,

    γ{x: |x−x₀|≤r, f(x)≤t₀+s} ≤ C_r exp[−(D−r)²/(2s)].             (1)

There is a reflected statement on the left. In particular an outward jump requires, on each side of its new endpoint,

    liminf_{δ↓0} log(1/δ) max_{0≤v≤δ}|f(x₀±v)−f(x₀)| ≥ D²/2.       (2)

To prove (1) take the weak left limit Q_{t₀−}=p_{t₀}−∫_{f<t₀}P_{t₀−f(x)}δ_x dγ(x); dominated convergence against bounded continuous tests proves this limit, and total variation is at most two. Its support is to the left of max(b,0) by the outer support bound. Continuity and extremality give f(x₀)=t₀; all points in the fixed r-neighborhood have f≥t₀. Fix η>0, y=x₀+η, λ=(D+η)/(D−r)>1. For small s, y is outside the hull at t₀+λs by right continuity. The heat semigroup identity then gives, at this y,

    (P_{λs}Q_{t₀−})(y)=∫_{t₀≤f(x)≤t₀+λs} p_{t₀+λs−f(x)}(y−x) dγ(x).

Any age-zero source lies away from y and contributes zero there. The left side is at most 2(2πλs)^(−1/2)exp[−(D+η)²/(2λs)]. For the selected r-neighborhood sources of height ≤t₀+s, the ages belong to [(λ−1)s,λs] and |y−x|≤η+r. Their kernels are at least (2πλs)^(−1/2)exp[−(η+r)²/(2(λ−1)s)]. All other source terms are positive. The difference of the two exponential rates is exactly (D−r)²/2 by the choice of λ. In fact C_r=2 works, with the threshold for s depending on r and η.

For (2) fix a small r. The one-sided interval of radius δ<r has Gaussian mass ≥cδ and belongs to the sublevel with s equal to the maximum indicated in (2). This maximum cannot vanish since (1) holds as s↓0. Taking logarithms of (1), then δ↓0 and finally r↓0, proves (2). For a nonzero zero x₀ of f the same conclusion holds with D=|x₀| by the original exponential identity directly: on the small interval choose the favorable signed parameter of magnitude (|x₀|−r)/s, and use f≤s. A steep contact is forced even along a single flank.

**Consequence for arbitrary shapes.** Suppose at every nonzero point x at least one flank satisfies

    liminf_{δ↓0} log(1/δ) max_{0≤v≤δ}|f(x±v)−f(x)|=0.               (3)

The flank can depend on the point. This holds in particular for locally Hölder graphs on the punctured line; no modulus is prescribed at zero. Then zero is the unique minimum and (2) excludes every jump of either hull boundary before T. Both endpoints tend to zero at time zero by compactness. Their continuity now forces every early sublevel to equal its hull: each positive intermediate position is swept by the right endpoint at an earlier or equal time when it belongs to the closed sublevel, and likewise on the left. The global excess budget and all polynomial tests of the preceding connected-sublevel theorem therefore hold without any a priori geometric hypothesis on the sublevels.

In particular if f tends on both tails to H=sup f, any nonconstant exact graph would have H>1. Taking T=H−1 makes the integrand in the global excess budget nonpositive everywhere, contradicting its strictly positive lower bound. In the original variables: if a is continuous, tends at both ends to its global infimum, and obeys (3) off the origin, an exactly Gaussian output forces a constant. No restriction on the number or spacing of peaks and valleys is required.
