# A convex-hull support restriction for signed intermediate heat mass

Let f≥0 be bounded and satisfy the complete identity ∫exp(zx−f(x)z²/2)γ(dx)=1 for complex z. For 0≤t≤sup f, define the finite signed measure
Q_t=P_tδ₀−∫_{f(x)≤t}P_{t−f(x)}δ_x γ(dx).
Then supp Q_t is contained in the closed convex hull of the essential support of the input measure γ|_{f>t}. In particular if a tail alone approaches an unattained supremum and the opposite tail has a strict gap, then for t sufficiently close to sup f this measure is identically zero on a half-line, without any assumption of positivity of Q_t.

For real s, the signed measure has Laplace variation bound
∫e^{−sy}|Q_t|(dy)≤ e^{t s²/2}+∫_{f≤t}e^{−sx+(t−f(x))s²/2}dγ(x) ≤2e^{(1+t)s²/2}.
Using the full identity gives the entire transform
∫e^{−zy}Q_t(dy)=∫_{f>t}e^{−zx−(f(x)−t)z²/2}γ(dx).
Suppose the future measure on the right is supported on [a,∞). Let F(z)=e^{az}∫e^{−zy}Q_t(dy). Along s≥0 one has 0≤F(s)≤1. Throughout the complex plane the variation bound gives |F(z)|≤C_a exp(C_t(Re z)²) by completing a square. We show this suffices to bound F throughout the right half-plane.

For 0<ε<π/4 let φ=π/2−ε, B_ε=C_t sin²ε/sin(2ε). On the sector 0≤arg z≤φ the function G=F(z)exp(i B_ε z²) is bounded by max(1,C_a) on the two rays: modulus of the extra factor on the second ray is exp(−B_ε |z|² sin(2φ)). Inside the sector it has growth at most exp(C'|z|²). Since π/φ>2, the sector maximum principle implies the same uniform bound in its interior. For a direct proof choose 2<p<π/φ, multiply G by exp[−δ exp(−ipφ/2)z^p] and use the maximum principle on large truncated sectors, then let δ↓0. The real part of exp(−ipφ/2)z^p is bounded below there by cos(pφ/2)|z|^p>0. Sending ε↓0 gives |F|≤max(1,C_a) in the first quadrant because B_ε→0; conjugation covers the fourth quadrant.

To convert this bound to support let σ be the translate of Q_t by −a and ψ smooth with compact support in (−∞,0). For s>0, Fourier inversion on the finite measure e^{−sy}σ(dy) gives
|∫ψ(y)σ(dy)|≤(sup_u |F(s+iu)|)/(2π) ∫| Fourier[e^{s·}ψ](u)|du.
Since ψ is supported in y≤−δ<0, its weighted test and its second derivative have L¹ norms O((1+s²)e^−δs). Two integrations by parts bound the displayed Fourier L¹ norm by the same order. Let s→∞, giving zero for every such test. This proves supp Q_t⊂[a,∞). Reflection gives the upper endpoint, and the empty-future case is immediate from Fourier uniqueness. The argument does not assert nonnegativity or erase bounded gaps inside the convex hull.
