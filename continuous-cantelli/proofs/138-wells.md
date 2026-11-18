# A simultaneous budget for all exterior wells of a canonical core

Use the conclusion (1) of round 137 at a positive radius b and time T>0. Set v(x)=f(x)−T on |x|>b. The canonical surviving measure has a density η(y) bounded by C(b−y) near the right edge, by reflection at the static upper wall b; it is bounded everywhere. Reflect for the other edge. For −1<p<1, Tonelli applied to (1) gives the exact positive budget

    ∫_{|x|>b} ∫₀^∞ z^p exp[z(x−b)−v(x)z²/2] dz γ(dx)
       = Γ(p+1)∫_{−b}^b η(y)/(b−y)^(p+1) dy.                    (1)

The right side is finite, and bounded by C/(1−p) uniformly for 1/2≤p<1. The identity charges all exterior wells against a single finite quantity, even for infinitely many components. It implies v>0 almost everywhere on x>b. Moreover (1) of round 137 and the same density bound give the pointwise right side after multiplication by exp(−bz) at most C/z² for z≥1.

Fix d>0, A=d²/2, and let M(s)=γ{x≥b+d:0≤v(x)≤s}. The following two estimates have different quantifiers:

    M(s) ≤ C_d s² exp(−A/s)                                      (2)
    ∫_ε^{s₀} [M(s) exp(A/s)/s^(5/2)] ds/s ≤ C_d log(1/ε).       (3)

Here s₀>0 can be fixed sufficiently small and (3) holds as ε decreases to zero. In particular liminf as s decreases to zero of M(s) exp(A/s)/s^(5/2) is finite. No monotone flanks, isolated well or first-birth assumption enters.

For (2) use z=d/s in the pointwise estimate. For (3) integrate the inner z integral in (1) just on [d/v,d/v+1/√v] at x≥b+d, v>0 small. Completing the square gives a lower bound c_d v^(−p−1/2)exp(A/v), uniformly for 1/2≤p<1. Write B_p(s)=s^(−p−1/2)exp(A/s). Since −B_p'(s)≥A s^(−p−5/2)exp(A/s), the layer-cake identity for the positive decreasing B_p implies

    ∫₀^{s₀} M(s) exp(A/s) s^(−p−5/2) ds ≤ C_d/(1−p).

Take p=1−1/log(1/ε). For ε≤s≤s₀<1 we have s^(1−p)≥exp(−1), giving (3). For example no candidate can have M(s)≥c s^α exp(−A/s) for every small s with α<5/2.

**Critical cusp on a remote well.** If x₀=b+d and f(x₀)=T, a putative outer flank cannot satisfy, for all sufficiently small u>0,

    f(x₀+u)−T ≤ A/[L−β log L+C],     L=log(1/u),     β<5/2,     (4)

with any fixed C. Indeed take L_s=A/s+β log(A/s)+K for K sufficiently large. The denominator at L_s is ≥A/s and is increasing for all larger L. Hence an interval of length exp(−L_s) lies in the sublevel, giving M(s)≥c s^β exp(−A/s) by the positive Gaussian density there. This contradicts (3). Reflection treats left-hand outer flanks. The logarithmically averaged estimate allows exceptional height sequences; a pointwise s^(5/2) upper bound from round 114 required additional hypotheses.
