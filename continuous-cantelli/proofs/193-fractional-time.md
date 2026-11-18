# A fractional time difference at a regular maximum set

Let X,G be independent standard Gaussians, let g be the standard Gaussian density, and suppose X+sqrt(v(X))G has density p_H, where p_b is the centered Gaussian density of variance b, H>0 and v≥0 is bounded and continuous. Let μ be a probability on a compact set E contained in {v=0}, with c t^d ≤ μ(B(a,t)) ≤ C t^d for a in E, 0<t≤diam E, where 0<d<1 and diam E>0. Set Q=1−d and r=dist(x,E).

**Theorem.** For each fixed K>0 there exists 0<κ<1 depending only on K and the regularity data such that, if r≤K sqrt(v) throughout 0<r<δ, then
```
∫_{0<r<δ} v(x)^(−Q/2−κ) dx < ∞.
```
More generally this holds on the good portion r≤K sqrt(v) if the complementary portion of that collar has finite integral of r^(−Q−2κ). One may always shrink the collar in these statements; under the assumed lower bound the rest of the collar then costs a finite amount (on the good portion).

**Proof.** Put h_b(x)=∫p_b(x−a)dμ(a), j=Q/2. The ball upper bound extends to arbitrary centers and radii by changing the constant. Summing over annuli of width sqrt(b), and absorbing polynomial factors from differentiation in a wider Gaussian, gives for l=0,1
```
|∂_b^l h_b(x)| ≤ C_l b^(−j−l) exp(−r²/(8b))   (b>0).
```
If b is small and r≤K sqrt(b), the ball of radius sqrt(b) about a nearest point gives h_b(x)≥c_K b^(−j). Define F(b,x)=∫₀^∞(h_b−h_{b+t})(x)t^(−1−κ)dt. On 0<t<b use the age derivative bound. On t≥b retain the entire h_b term and use h_{b+t}≤C₀ t^(−j). Thus for the indicated centers and ages
```
F(b,x) ≥ b^(−j−κ) [c_K/κ − C₁/(1−κ) − C₀/(j+κ)].
```
Since j>0 this is strictly positive with a fixed margin at sufficiently small κ. Everywhere |F|≤C b^(−j−κ) by the same split. Also |F|≤C r^(−2j−2κ) for r>0: for t<r² bound all age derivatives by C r^(−2j−2) using the exponential estimate, and for t≥r² bound both kernels by C r^(−2j). All constants in the last two estimates may depend on the fixed κ.

Smooth the output by independent variance s>0. Fubini (ages ≥s, bounded derivative at short times and bounded kernel at long times) gives exactly
```
∫ g(x) F(v(x)+s,x) dx
 = ∫ dμ(a) ∫₀^∞ [p_{H+s}(a)−p_{H+s+t}(a)] t^(−1−κ)dt .
```
The right side is bounded uniformly for 0<s<1, using H>0. On a sufficiently small collar v+s is small for all sufficiently small s by continuity. On its good part the positive lower bound applies also after smoothing. On its bad part and outside the collar the negative part has a uniform integrable majorant by the distance estimate and the stated hypothesis. E has Lebesgue measure zero, by a disjoint-ball cover of cardinality O(t^(−d)). Fatou's lemma and the positive minimum of g on the bounded collar give the theorem.

**Consequence.** The tube volume satisfies c₁t^Q≤|{r<t}|≤C₁t^Q, by maximal separated subsets of E and the two mass bounds. Fix a>1 sufficiently large. Each shell t/a<r<t has length at least c₂t^Q. If on a fixed positive length fraction of infinitely many disjoint such shells v≤C t^p, the budget is impossible for p(j+κ)≥Q, provided the good-set hypothesis holds throughout the collar (or on those fractions with the integrable exceptional set above). In particular, when v≥c r² throughout, any upper bound v≤Cε r^(2−ε) for every small ε is impossible. This includes all two-sided r²(log(1/r))^β profiles for β≥0. The threshold for a strict power improvement can depend on E and K. No lattice, gap symmetry or uniform distribution of input across the collar is needed.
