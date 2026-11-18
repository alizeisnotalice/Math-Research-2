# A global sign law, two-point information orders and the heat extinction time

Write γ=N(0,1), G_t=N(t,1), and suppose f is bounded, nonnegative, nonconstant and continuous with ∫ exp(tx−t²f(x)/2)dγ=1 for all real t. Put h=f−1 and P_t(dx)=exp(tx−t²f(x)/2)γ(dx). These are statements for an arbitrary exact graph with no extrema, monotonicity, spatial derivatives or tail assumptions. All Gaussian differentiation used below is valid by boundedness of f. The first four derivatives at 0 give E_γ f=1, E_γ[Xf]=0, and E_γ[X²f]=(1+E_γ f²)/2.

## Strict sign at every spatial center and every coupling strength

    sign[∫exp(−b h)dG_t−1] = sign[b(b−t²/2)]                 (1)

for all real t,b. At fixed t the left expression before taking its sign is a strictly convex function of b (h is nonzero on a set of positive Gaussian mass). It vanishes at b=0 and b=t²/2 by the exact graph identity. For t=0 these coincide, and its derivative vanishes there since E_γ h=0. This proves (1), including the strict signs in all regions off the zero set. For each fixed b>0 the Gaussian convolution of exp(−b h)−1 is positive on the centered interval |t|<sqrt(2b) and negative outside its closure. For b<0 it is positive everywhere.

Let B=sup h>0, U(t)=∫h dG_t. There is also a quantitative local-variance restriction

    U(t) ≥ [(exp(−bB)−1+bB)/(bB²)] ∫h² dG_t,   b=t²/2>0.    (2)

Indeed exp(−b y)−1+by=b² y² ∫₀¹(1−s)exp(−sby)ds; the integral is decreasing in y≤B. Substitute the exact expectation. Parseval in the Gaussian Hermite basis gives ∫h²dG_t=∑_{j≥0}|U^(j)(t)|²/j!, since U^(j)(t)=∫He_j(x−t)h(x)dG_t(x). Thus (2) controls all translated squared jets at once.

## Two-point entropy and the Rényi sign reversal

For q(t)=E_{P_t}f one has E_{P_t}X=tq(t) by differentiating the normalization. Direct likelihood-ratio integration gives

    D(P_t || P_s) = (t−s)² q(t)/2,
    D(G_t || P_t) = t² U(t)/2,
    D(P_t || G_{tq(t)}) = t² q(t)[1−q(t)]/2.

All laws here have positive Lebesgue densities. For t≠0 the last relative entropy is strictly positive: equality forces P_t to be a unit-variance translated Gaussian, so its logarithmic density ratio against γ is affine in x and the tilt formula makes f affine. Boundedness then forces f constant. The same applies to the middle relative entropy. Therefore 0<q(t)<1<E_{G_t} f for t≠0.

For t≠s, α>0 with α≠1, put u=αt+(1−α)s, B_α=αt²+(1−α)s², and define A_α(P,Q)=∫p^α q^(1−α). All orders are finite since f is bounded and the log-likelihood ratios have at most linear growth. Completing the square gives exactly

    A_α(P_t,P_s) / A_α(G_t,G_s) = ∫exp[−B_α h/2]dG_u.

Since B_α−u²=α(1−α)(t−s)², (1) determines its sign relative to 1 at every pair of tilts. In terms of the Rényi divergence D_α=(α−1)^(−1) log A_α,

- for 0<α<1 it is strictly smaller for (P_t,P_s) than for (G_t,G_s);
- for α>1 it is strictly smaller when B_α>0, strictly larger when B_α<0, and equal on the two rays B_α=0.

In particular D_α(P_0||P_s)>D_α(G_0||G_s) for s≠0 and every α>1, whereas D_α(P_s||P_0)<D_α(G_s||G_0). The two parameterized experiments cannot be transformed one into the other by fixed Markov kernels in either direction. This follows already at α=2, because conditional expectation contracts the squared L² norm of a likelihood ratio under a Markov transformation.

## Extinction of the negative heat phase

Write u(w,t)=∫h(x)p_w(t−x)dx for w>0. Formula (2) makes u(1,t)>0 at t≠0 and u(1,0)=0. For w>1 the convolution semigroup makes u(w,t)>0 at all t. For every 0<w<1 the function u(w,·) takes both signs: its average against the strictly positive p_(1−w) is u(1,0)=0, and if it vanished everywhere then so would u(1,·). Thus time one is the exact first time at which the heat average of h is globally nonnegative. The same conclusion at w=0 follows from E_γ h=0 and nonconstancy.

Writing V=Var_γ f>0, the fourth moment identity gives u_tt(1,0)=V/2, the third gives u_t(1,0)=0, and the heat equation gives u_w(1,0)=V/4. Hence

    u(w,t)=(V/4)(t²+w−1)+o(t²+|w−1|).

For w<1 sufficiently close to one there are exactly two zeros in a fixed sufficiently small neighborhood of 0, with positions ±sqrt(1−w)(1+o(1)); the function is negative between them. Indeed strict positivity of u_tt throughout a small spacetime neighborhood makes it strictly convex there; use the expansion at t=0 and t=±(1+ε)sqrt(1−w), then at t=±(1−ε)sqrt(1−w). This describes the local components at extinction and does not presume that faraway negative components are absent at earlier times.

These are global necessary conditions and consequences for every nonconstant exact graph. The sign and information constraints alone have not ruled out such a graph.
