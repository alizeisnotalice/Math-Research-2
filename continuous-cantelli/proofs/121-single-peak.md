# The common cubic moment excludes unequal monotone wings

Let γ be standard Gaussian measure with density g. Suppose the bounded continuous normalized graph f≥0 satisfies ∫exp(zx−z²f(x)/2)dγ(x)=1 for every real z.

**Theorem.** Suppose f is nonincreasing up to its minimum and nondecreasing afterwards (the minimum interval and mode location are initially arbitrary). Then f is constant. In particular no finite positive separation of its two tail limits permits a nonconstant single-well graph.

**Reduction for all one-peak and one-valley conventions.** If the original coefficient is the conditional variance a (or the nonnegative standard deviation), a single peak becomes a single valley of f=1+Ea−a. A single valley of the variance becomes a single peak of f. Passing from a nonnegative standard deviation to its square preserves all weak orderings. These definitions allow plateaus and endpoints at infinity.

For a bounded continuous normalized hill or valley the two one-sided limits exist. The R083 identities for a nonconstant exact f give Ef=1, E[Xf(X)]=0, and E_{N(z,1)}[f−1]>0 for every z≠0 (strict Jensen). Letting z tend to either infinity shows that both limits are at least one, by dominated convergence after translation. If f is a hill, f is everywhere at least the smaller of those two limits, so Ef=1 forces f≡1. A hill or valley with its extremum only at infinity is globally monotone. Indeed for a valley with infimum approached along a sequence tending to +∞, quasi-convexity gives f(y)≤max(f(x),f(x_n)) for x<y<x_n, and then f(y)≤f(x). Reflect the line or reverse the values for the other possibilities. A globally monotone nonconstant f has strictly signed Cov(X,f(X)) by the independent-copy identity and the full support of γ, contradicting E[Xf(X)]=0.

For the remaining finite valley, each wing limit strictly exceeds one. If a wing tends to one, it is bounded above by one all the way to the minimum. Since a nonconstant exact graph has min f<1, there is a fixed interval with f−1 strictly negative lying farther toward that wing than any point where f−1 is positive. Its contribution to ∫exp(zx)(f−1)dγ dominates all positive contributions as z tends to infinity toward that wing, contradicting strict Jensen; if the minimum interval has positive length choose the negative test interval on the indicated side of the opposite flank. This is R083 §6. With both far-tail limits strictly greater than one, a positive compact minimum implies that the exact transform tends to zero: on the far tails use f≥1+δ and the full Gaussian moment generating function, on each compact use the positive minimum. Hence the minimum is zero. If the zero set contained an interval, a smaller interval away from zero would make the favorably tilted integral unbounded. The zero is therefore unique. If its location c were nonzero, tilt away from c: its compact neighborhood contributes at most exp(−|zc|/2), the remaining compact pieces have a fixed positive minimum, and both far tails again decay. This contradicts exactness and places the unique zero at the origin. These arguments also record where R083's early-band hypotheses come from before using symmetry.

**Proof of the tail contradiction.** Assume nonconstant. R083 §6 gives finite tail limits H₋,H₊ strictly greater than one, a unique zero at the origin, and the Brownian first-exit early band for all 0<t<T=min(H₋,H₊)−1. Reflect to put H₋≤H₊, and let L=H₋=1+T. By R088 §1 the early interval is symmetric, with endpoints ±b(t), and f(r)=f(−r) for r<b(t). Its increasing limit b at T is finite and strictly positive (the two tail limits exceed T). Hence f(r)=f(−r) for 0<r<b. The R084 terminal budget states

    ∫_{|x|>b} (f(x)−L)dγ(x) ≥ 2b [g(b)−b barΦ(b)] >0.            (1)

Only t<T is used to obtain this by a limit. No sign of Q_{T+} or regularity at the splitting boundary is assumed.

On the other hand the cubic coefficient of the original exact transform is zero. Since the centered Gaussian has third moment zero, this says ∫x f(x)dγ(x)=0. Cancel the symmetric interior to obtain, with R having the conditional law of X given X>b, and F₊(r)=f(r), F₋(r)=f(−r),

    E[R F₊(R)] = E[R F₋(R)] ≤ L E[R].

The inequality uses only F₋≤L for r>b. By monotonicity of the high wing, Cov(R,F₊(R))≥0: twice this covariance equals E[(R−R')(F₊(R)−F₊(R'))] for independent copies. Therefore E F₊≤E[R F₊]/E[R]≤L. Also E F₋≤L. The left side of (1) is nonpositive, a contradiction. The constant normalized graph is f≡1 by the second moment. This proof includes the equal-tail case.

**General form.** The contradiction uses only a symmetric positive early band up to T with limiting radius 0<b<∞, a low side f(−r)≤1+T for r>b, and Cov(R,f(R))≥0 on the opposite tail. No global shape on the outer wings is required beyond these last two assumptions.

There is also a quantitative overshoot formula when the low side can rise above L. For 0<t<T write b=b(t), m=E[R]=g(b)/barΦ(b), Y for the survivor Q_t normalized by its mass 2barΦ(b), and C₊=Cov(R,F₊(R)). The cubic identity still gives E[R F₊]=E[R F₋]. The exact future second moment gives

    ½ E(F₊+F₋)=1+t+bm−EY²,

since E[R²]=1+bm. Eliminating the unweighted high-wing mean using the cubic identity gives the exact identity

    E[(1+R/m)(F₋−1−t)] = 2b(m−b) + 2(b²−EY²) + C₊/m.             (2)

Here b²−EY²≥0 by the support of the positive survivor. In particular, whenever the high-wing covariance is nonnegative, the remaining low wing must exceed 1+t in a weighted mean by at least 2b(m−b). For t↑T the strict endpoint margin in (1) remains positive. The same underlying two wings satisfy these constraints at every early time, and the odd moment enforces the same weighted mean in both wings at every cut.

**Size of a necessary descent on the high side.** Keep the finite terminal radius b and F₋≤L on r>b; drop monotonicity of F₊. In (2) let t↑T, using its inequality with b²−EY² omitted. This forces

    −Cov(R,F₊) ≥ 2mb(m−b),       m=g(b)/barΦ(b),  R∼X | X>b.       (3)

Put D_b=E[(R−m)_+]=(g(m)−m barΦ(m))/barΦ(b). Since E(R−m)=0, the magnitude of this negative covariance is at most D_b times the essential oscillation of F₊ on r>b. Its oscillation must consequently exceed or equal 2mb(m−b)/D_b. If F₊ has bounded variation, (3) gives the stronger bound

    ∫_{(b,∞)} [g(u)−m barΦ(u)]/barΦ(b) dF₊^−(u) ≥ 2mb(m−b),

where dF₊^− is the downward variation measure. Indeed integrate the conditional covariance kernel K(u)=E[(R−m)1_{R>u}] against dF₊. The kernel is zero at b and at infinity, increases to its maximum D_b at u=m, and is strictly positive in between. Thus the required weighted descent must occur on the active tail. For finite t<T the same inequality holds with 2mb(m−b) decreased by 2m(T−t), giving a family of tests of the same F₊ using the cut radii b(t).
