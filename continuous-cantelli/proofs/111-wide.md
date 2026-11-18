# Wide windows, exponentially weighted mass and the deepest wells

Use round 101 notation. For b>0 let j_b(x)=J(1,x,b), R=√(2b), and define T_b(λ)=∫exp(λx) j_b(x)dx, allowing −∞. Its positive part has compact support. Then the limit β∞=lim_{w→∞}β(w) exists in (0,∞] and

    β∞≤b  if and only if  T_b(λ)≥0 for some real λ.

Indeed the sign of J(1+d,dλ,b) is the sign of
A_d(λ)=∫exp(λx+(R²−x²)/(2d)) j_b(x)dx.
For fixed λ this is strictly decreasing in d to T_b(λ), by the exact two signs of j_b and monotone convergence of each part. For all d≥1 the closed sets {λ:A_d(λ)≥0} lie in a single compact interval: the positive part on [−R,R] has uniformly bounded mass before exponential tilting, while each of two fixed intervals strictly beyond ±R has a uniformly positive absolute negative mass. If β∞≤b then all these sets are nonempty by round 101. Their decreasing intersection gives the asserted λ. Conversely T_b≥0 at some λ forces positive J somewhere for every finite d. This proves the equivalence even when the tail integrals diverge. In particular, whenever ∫exp(λx)|h(x)|dx is finite, T_b(λ)=exp(λ²/2)∫exp(λx)(exp(−bh(x))−1)dx/b by absolute Fubini. If h belongs to ordinary L¹, β∞ is finite: for large b the negative h on a fixed interval forces this last integral at λ=0 to become positive, while its negative part has size at most ∫|h|.

**Depth bound without tail hypotheses.** Put D=−inf h>0. For every ε∈(0,D), take a fixed finite interval I of positive length on which h≤−D+ε. Then P_w exp(−β(w)h)(0)≤1, and the contribution of I is at least c_I w^(−1/2) exp((D−ε)β(w)) for large w. Thus limsup β(w)/log w ≤1/(2D).

**Sharp scale with separated positive tails.** If h≥η>0 outside a compact interval, then β(w)→∞. For each fixed b>0 the outside value exp(−bη) is below one, and the compact excess contributes O_b(w^(−1/2)) uniformly in the center. At the critical b=β(w), define M(b)=∫_{h<0}exp(−bh(y))dy. Uniformly in the center the contribution from the outside tails is o(1), and that from the compact nonnegative region is O(w^(−1/2)). The maximum of the contribution from the compact negative region lies between (1−O(1/w))M(b)/√(2πw) (test the center zero) and M(b)/√(2πw). At threshold the total maximum is one, hence M(β(w))∼√(2πw). The elementary Laplace principle then gives β(w)∼log(w)/(2D).

For example, if the finitely many global minima x_i satisfy h(x_i+u)=−D+c_i|u|^{p_i}+o(|u|^{p_i}), c_i,p_i>0, put p=max p_i and C=Σ_{p_i=p}2Γ(1+1/p)c_i^(−1/p). Local rescaling using comparison to exponents c_i±ε and exponential suppression off the minima gives M(b)∼C exp(Db)b^(−1/p). Therefore

    Dβ(w)=(log w)/2+(log log w)/p+log(√(2π)/C)−log(2D)/p+o(1).

These are conditional necessary relations for exact graphs. The multiwell Laplace formula describes the threshold if these hypotheses hold; it does not assert existence of a nonconstant exact graph.
