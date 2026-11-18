# Complex Gaussian averaging on a homogeneous Cantor set of terminal maxima

Let E be the attractor of m≥2 equal-ratio similarities S_i on the line, each of contraction ratio ρ∈(0,1), with strongly separated first-level intervals inside the convex hull I=[l,r] of E. Assume mρ<1. The dimension is d=log(m)/log(1/ρ)∈(0,1). Orientations may be reversed. Let f be bounded, continuous and nonnegative, H=sup f, v=H−f and Z={v=0}. Suppose E⊂Z and the output from standard normal input X, with independent standard normal G, is X+√v(X)G~N(0,H). Write g=p₁.

Suppose that on a punctured neighborhood of E there is one positive function V with

    v(x)=V(dist(x,E)),     V(h)[log(1/h)]^β→A∈(0,∞),     β≥1.              (1)

For β>1 these assumptions are incompatible. For β=1 the same conclusion holds if dist(E,Z outside E)>√(2A(1−d)), with distance to the empty set infinite. No lower bound on v is imposed at spatial infinity. Formula (1) specifies exact radial dependence through V, while its logarithmic asymptotic coefficient can have an arbitrary uniform relative o(1) error as a function of distance. An address-dependent relative o(1) error in v is not covered by the main statement.

## Choice of test

For s>0, exactness implies ∫p_{v(x)+s}(z−x)g(x)dx=p_{H+s}(z) throughout the complex plane, by Gaussian smoothing and dominated holomorphic extension on compact sets. Let μ be the equal-weight self-similar probability on E, and let σ(da)=g(a)μ(da). Average the identity over z=a+iD against σ. For bounded D and s≤1 the target stays bounded uniformly, since H≥1 by the normalized second-moment identity. For any b>0 the corresponding model with the same real spatial measure in both arguments is nonnegative:

    M(b,D)=∬p_b(a+iD−x)σ(da)σ(dx)
          =(2π)^(−1)∫_ℝ exp(−bξ²/2+Dξ)|σ̂(ξ)|²dξ ≥0.                      (2)

Fourier inversion is absolutely integrable for fixed b. The sign of Dξ depends on convention and can be reflected since |σ̂|² is even.

Choose c strictly between 1−d and 1. At β=1 choose c also small enough that D=√(2cA)<dist(E,Z outside E). At each large integer j set L=j log(1/ρ) and take

    β=1:     D_L=√(2cA),                s_L=A L^(−2);
    β>1:     D_L=√(2cβA L^(1−β)),       s_L=(β−1)A L^(−β).                  (3)

In both cases D_L stays bounded. For some c' with c<c'<1 and all sufficiently small h (uniformly in L),

    exp(D_L²/[2(V(h)+s_L)]) ≤ h^(−c'),
    (V(h)+s_L)^(−3/2) ≤ C [log(1/h)]^(3β/2).                              (4)

For β=1 this follows directly from (1). For β>1 put U=log(1/h). The inequality β t^(β−1)≤1+(β−1)t^β for t=U/L>0 follows by differentiation or weighted arithmetic-geometric means. It gives D_L²/[2(A U^(−β)+s_L)]≤cU, with equality at U=L. In (1) replacing A by (1−η)A in the variance bound costs at most the factor 1/(1−η); choose η so c/(1−η)<c'. Decrease h₀<1 so (1) and (4) apply for 0<h<h₀.

## Comparison with positive model terms

The components of I outside the first-level intervals form finitely many basic open gaps G. Every x∈I outside E and the countable family of endpoints belongs to S_ω(G), for one word ω of length n≥0. Write x=S_ω(u) with u∈G and h_G(u)=dist(u,∂G). Since the two endpoints belong to E and no point of E is inside this gap, dist(x,E)=ρ^n h_G(u). The Jacobian is ρ^n. Lebesgue measure of E is zero because mρ<1.

For ρ^n h_G(u)<h₀, let b=s_L+V(ρ^n h_G(u)). It is independent of ω. Compare the average over all words of length n of
g(S_ω(u))∫p_b(a+iD_L−S_ω(u))σ(da) with M(b,D_L). Indeed if ω is uniform over these words and Y has law μ independently, S_ω(Y) has law μ by self-similarity; also |S_ω(u)−S_ω(Y)|≤|I|ρ^n. For x,a∈I, g and g' are bounded and

    |p_b(a+iD_L−x)|+|∂_x p_b(a+iD_L−x)|
       ≤ C(1+b^(−3/2)) exp(D_L²/(2b)).

The coupling therefore incurs error at most Cρ^n(1+b^(−3/2))exp(D_L²/(2b)). Using the same weight g in both measures of (2) is essential for positivity. Summing over the words multiplies the average by (mρ)^n after the change of variables. By (4), the sum of integrated absolute errors is bounded uniformly in L by a constant times

    Σ_G Σ_{n≥0} (mρ)^n ρ^n ∫_G
      (ρ^n h_G(u))^(−c') [1+n+|log h_G(u)|]^(3β/2)du.                       (5)

The u integral with h_G(u)^(−c') is finite even at the two endpoints since c'<1. The geometric factor in n is ρ^((2−d−c')n); its exponent is positive. This proves the uniform bound, including points arbitrarily close to endpoints of the coarse gaps.

The two parts of the punctured neighborhood just outside I have distances x−r or l−x from E. Their full absolute contributions are bounded uniformly by integrating h^(−c')[1+log(1/h)]^(3β/2) over (0,h₀). This also uses (4). For the remaining inputs dist(x,E)≥h₀, consider first |x−a|>D_L+η for fixed small η>0. The absolute Gaussian kernel is bounded uniformly in b because (x−a)²−D_L²≥η². For β>1 one can take D_L<h₀/2 eventually, so this bound covers the entire remaining set. For β=1 choose η so D_L+2η is less than the distance in the hypothesis. The compact annulus h₀≤dist(x,E)≤D_L+η contains no variance zeros and v has a positive minimum there. This completes the uniform bound for all remaining sources, including any distant maxima or small variances in the tails.

By (2), every retained model term is real and nonnegative. For fixed L the integrals and word sums are absolutely convergent (s_L>0). The bounded target and the uniformly bounded errors just established imply

    Σ_G Σ_{n≥0} (mρ)^n ∫_{u∈G:ρ^n h_G(u)<h₀}
          M(s_L+V(ρ^n h_G(u)),D_L)du ≤ C.                                   (6)

## One layer violates the bound

The measure σ is nonzero, finite, compactly supported and singular. There exist arbitrarily large ξ>0 with |σ̂(ξ)|≥ξ^(−1). Otherwise the Fourier transform is square integrable and Plancherel gives an L² density for σ. Since ∫|a|σ(da)<∞, σ̂ has bounded derivative. Thus along such a sequence |σ̂(k)|≥(2ξ)^(−1) throughout a frequency interval of length at least c₁/ξ next to ξ.

Set k_*(L)=D_L L^β/(βA), a positive constant times L^((β+1)/2). For each frequency ξ of this sequence, take j to be the nearest integer to k_*^(−1)(ξ)/log(1/ρ) and put L=j log(1/ρ). Then j→∞ and ξ/k_*(L)=1+O(1/L). Fix a closed subinterval of positive length strictly inside a basic gap. For u there, log(1/(ρ^j h_G(u)))=L+O(1) uniformly and

    b=s_L+V(ρ^j h_G(u))=β A L^(−β)(1+o(1)).

For k in the frequency interval just obtained,

    D_L k−b k²/2=cL+o(L),

uniformly in u and k. Indeed D_L k_*=2cL and βA L^(−β)k_*²/2=cL. Integrating (2) over this frequency interval gives M(b,D_L)≥C ξ^(−3)exp(cL−o(L)). The contribution to (6) at stage n=j over the chosen basic subinterval is at least

    C ξ^(−3) exp[(c−1+d)L−o(L)],

since (mρ)^j=exp[−(1−d)L]. This goes to infinity and contradicts (6).

The intermediate strict choices c<1 and c>1−d use 0<d<1. This proof gives an exclusion for a spatially extended cluster even when no individual maximum is isolated. It does not require Fourier decay or an arithmetic relation between contraction ratios and the frequency sequence; the same contraction ratio at every address of a given level was used to keep b fixed within each average.
