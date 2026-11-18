# Sparse far-tail windows and the second-order threshold

Use the notation, assumptions and inequality (1) of [round 98](098-directed-enthalpy.md).

A full low-tail limit can also be replaced by sparse visible critical dips. Keep the early two-sided lower limits and the high-side gap. Suppose on a sequence λ_j→∞ there are measurable E_j⊂[λ_j−1,λ_j+1] of length exp(−o(λ_j²)) where T+h₀<f(−r)≤1+T+o(1) uniformly. Use just these points in the proof of (1), at λ=λ_j. Their transformed mass is exp(−o(λ_j²)) or greater, again forcing Δ>0 arbitrarily close to b.

**Second order at the slow-tail threshold.** Suppose 0<p≤1/2 and
f(−r)≤L+B r^{−p}+C r^{−2p}+o(r^{−2p}).
Complete the Gaussian square and restrict the far integral to a unit window about r=λ+(Bp/2)λ^{1−p}. Taylor expansion of r^{−p}, r^{−2p}, and the Gaussian exponent on this unit window gives
log I_far ≥ −Bλ^{2−p}/2 + [B²p²/8−C/2+o(1)]λ^{2−2p}.
Indeed an offset vλ^{1−p} contributes (Bpv/2−v²/2−C/2)λ^{2−2p}, maximized at v=Bp/2.
Take γ=(2−p)/p and k=B^{1+γ}γ^γ/[2(1+γ)^{1+γ}]. If
Δ≤exp(−k/h_*^γ+K₁/h_*^{γ−1}+o(h_*^{1−γ})),
then it is necessary that
K₁ ≥ [Bγ/(1+γ)]^{γ−1} (B²p²/8−C/2−b·1_{p=1/2}).
For details the minimizing h_* in qh_*+k/h_*^γ is h₀=[Bγ/(1+γ)]λ^{−p}. Upon setting h_*=h₀ z the leading phase has a unique strict minimum at z=1 and diverges at both extremes. The order λ^{2−2p} correction K₁/h_*^{γ−1} is lower order on compact z-intervals; on z<1/2 it is absorbed into an arbitrarily small fraction of k/h_*^γ when h_* is small, and for z>2 it is bounded by a constant times λ^{2−2p}, dominated by the leading fixed gap. Thus the optimizing z tends to one and the corrected upper coefficient is at most K₁[Bγ/(1+γ)]^{−(γ−1)}. The spatial factor contributes at most (b+o_h(1))λ; this is negligible for p<1/2 and at p=1/2 one sends h↓0. In particular for p<1/2, C<B²p²/4 rules out the critical upper bound Δ≤exp(−k/h_*^γ) itself.
