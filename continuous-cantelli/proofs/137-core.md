# One-sided domination identifies a whole canonical core

Let f be a bounded nonnegative continuous exact graph and let c(r) denote the inverse radius of the symmetric standard-Gaussian Brownian first-exit boundary used in round 136. Choose any b>0 in its range, T=c(b)>0. Assume only

    f(x)≤c(|x|) for |x|<b.

Then f=c throughout [−b,b], f≥T everywhere outside it, and for every real z, with q=z²/2,

    ∫_{|x|>b} exp(zx−q(f(x)−T)) γ(dx) = ∫_{−b}^b exp(zy) η_T(dy).       (1)

Here η_T is the true canonical surviving subprobability at time T. In particular all strict-past signed heat slices before T agree with canonical survival. These assertions impose no conditions on the number of additional wells, tail limits or a modulus at any point.

Indeed, stopping the Brownian exponential martingale at the minimum of T and the canonical exit time gives the sum of the canonical local integral and exp(−qT) times the right side of (1) equal to one. Its survivor term is at most exp(b|z|−qT). By exactness, the nonnegative excess of the local integral for f over the canonical one is bounded by that term. If f<c at an interior point, there is a closed interval of positive length with f≤T−ε and c−f≥ε' for positive ε,ε'. Its excess is at least a constant times exp(−b|z|−q(T−ε))(1−exp(−qε')), which contradicts the upper bound as z tends to infinity. Continuity gives the endpoint values. Subtraction now gives (1) exactly. Any point outside with f<T would likewise give a fixed compact interval whose left contribution grows like exp(εq−O(|z|)), contradicting the supported survivor. The stopped sources for earlier times therefore coincide pointwise with the canonical ones, proving the assertion on heat slices.

In particular if the first radius where either wing rises strictly above c has infimum b>0, all earlier heights are canonical, and every discrepancy below the canonical profile further out has height at least c(b). The result is valid even if f has infinitely many further local minima. Round 136 covered two-sided ordered germs under high-tail and isolated-zero assumptions; for this direction of ordering the full interval and all exterior mass are identified.
