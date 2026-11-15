# R60: exact rectangle means for joint Cauchy kernels and the Jensen gap

2026-11-15. Independent bounded analytical/Python support for a joint radius/shape/time search. This note creates no new candidate and no new general C₂ upper bound. R55 files are preserved; no MATLAB or LP was run, and no main-workspace/commit/push operation was performed.

## 1. Entire rectangle integrals reduce to one dimension

Use U(x)=1/(|x₁|+|x₂|), coordinate-Cauchy density p_s, and V_t=U*p_t=2π∫_t^∞p_s ds. Let B=[l₁,h₁]×[l₂,h₂] with positive widths. For a source z, replace l_i,h_i by l_i−z_i,h_i−z_i below. Define

Δ_i(s)=atan(h_i/s)−atan(l_i/s)
      =atan2(s(h_i−l_i),s²+l_ih_i),  s>0.

The atan2 expression lies in (0,π), including intervals crossing zero, and avoids cancellation for narrow far intervals. Positivity and Tonelli give the exact formula

V_B(t):=∫_B V_t(x)dx=(2/π)∫_t^∞Δ₁(s)Δ₂(s)ds.              (1)

The complementary positive occupation integral is

Q_B(t):=∫_B(U−V_t)dx=(2/π)∫₀^tΔ₁(s)Δ₂(s)ds.              (2)

At s→0+, Δ_i tends to π if l_i<0<h_i, to π/2 if exactly one endpoint is zero, and to zero otherwise. Thus source points on cell corners/edges are legal and all rectangle integrals are finite. The t=0 value in (2) is zero. Neither representation discards the potential singularity.

For t>0 transform s=t/u and define

D_i(u)=atan2(tu(h_i−l_i),t²+l_ih_i u²)/u,  u>0,
D_i(0)=(h_i−l_i)/t.

Then

V_B(t)=(2t/π)∫₀¹D₁(u)D₂(u)du.                            (3)

The stable endpoint limit is explicit. More importantly,

D_i(u)=(1/t)∫_(l_i)^(h_i)[1+(xu/t)²]⁻¹dx

is nonnegative and decreasing on [0,1], for any position of the interval relative to zero. Formula (2), after s=tu, is often better for a small-t far cell because it avoids subtracting nearly equal U and V means.

The U rectangle mean is elementary. With 0 log0=0, set

H(x,y)=sgn(x)sgn(y)[(|x|+|y|)log(|x|+|y|)−|x|log|x|−|y|log|y|].

Then

U_B=H(h₁,h₂)−H(l₁,h₂)−H(h₁,l₂)+H(l₁,l₂),                (4)

because H is the integral from 0 to x and 0 to y of U. This remains valid through either axis. Therefore a normalized finite orbit σ_j has the exact basis mean

c_(B,j)=U_B−average_(z∈orbit_j)V_(B−z)(t_j).

Equivalently it is the compact atomic potential-difference mean plus the nonnegative average of Q_(B−z)(t_j). The latter separates the small-time contribution; cancellation in the compact log formula at very remote cells still needs precision or stable algebra.

## 2. Direct rigorous enclosures and a faster quadrature sandwich

For g(u)=D₁(u)D₂(u), monotonicity gives lower/right and upper/left endpoint Riemann sums on any partition. At n uniform subintervals their exact gap is

(2t/(πn))[g(0)−g(1)]≤2|B|/(πtn).                         (5)

Directed arctangent and π evaluation turns these algebraic inequalities into a rigorous computable enclosure; an ordinary floating Riemann sum does not. Adaptive dyadic subdivision can focus only where the upper/lower products differ.

A faster structural alternative follows by v=u². Put

F(v)=D₁(√v)D₂(√v)
    =(1/t²)∫_(l₁)^(h₁)∫_(l₂)^(h₂)[(1+x²v/t²)(1+y²v/t²)]⁻¹dxdy.

F is completely monotone: (−1)^kF^(k)≥0 for every k. Indeed each reciprocal factor has that property and Leibniz products preserve it, with differentiation under a bounded finite-domain integral. Formula (3) is

V_B(t)=(2t/π)∫₀¹F(v)dv/(2√v).

For this positive probability weight, n-node Gaussian quadrature is a LOWER bound, and n-node left Gauss–Radau (including v=0) is an UPPER bound. A proof avoids relying on a remembered sign convention: Hermite interpolation at all n Gaussian nodes has remainder F^(2n)(ξ) times a squared product, hence nonnegative; left-Radau interpolation has remainder F^(2n−1)(ξ)·v times a squared product, hence nonpositive. Integrate the interpolating polynomial, using exactness through degree 2n−1 and 2n−2 respectively.

With SciPy's convention weight (1−x)^α(1+x)^β on [−1,1], Gaussian nodes are `roots_jacobi(n,0,−.5)` mapped to [0,1], weights divided by 2√2. Left-Radau interior nodes come from `roots_jacobi(n−1,0,+.5)`; their mapped weights are w_i/(4√2 v_i), and the endpoint weight is 1−Σ_iλ_i. These are exact mathematical prescriptions; float nodes/weights are not directed bounds.

Whole-interval fixed n can be inefficient when max|endpoint|≫t: the v=0 boundary layer then has width about (t/max|endpoint|)². A tested improvement splits at ε=min(1,(t/max|endpoint|)²). On [0,ε], scale the Jacobi/Radau rule, including factor √ε. On successive geometric intervals [ε,4ε], …, ending at 1, use ordinary Gaussian and left-Radau rules for F(v)/(2√v). This integrand is itself completely monotone on v>0, so the same remainder signs apply. No claim of an implemented directed enclosure is made here.

## 3. The exact Jensen LP objective, including all omitted tails

Let α_j≥0 and W=Σ_jα_jW_j, W_j=U−V_(t_j)*σ_j, where each σ_j is a normalized centered even probability orbit. Then W is an integrable legal heavy-tail potential difference, its positive part has distributional C4, and

∫W=2πΣ_jα_jt_j=:S>0.

For disjoint bounded cells B_i let I_i=∫_(B_i)W=Σ_jα_jc_(i,j). On their union Jensen gives

Σ_i max(I_i,0)≤∫_(∪B_i)W_+.

There is a useful exact complement row even if that complement is unbounded:

c_(out,j)=2πt_j−Σ_i c_(i,j),
I_out=S−Σ_iI_i.

Hence the complete finite-partition objective

J_partition(α)=Σ_i max(I_i,0)+max(I_out,0)≤∫W_+.           (6)

Its half is a lower objective for kernel comparison cost when the D normalization W≥1 is valid. Epigraph variables for each finite cell and the complement preserve the LP. Cell integrals include their actual areas, and any D4 multiplicities must be applied to the finite sum before forming the complement. For a wedge cell copied eight times, its objective contributes eight times the cell positive mean; its complement signed coefficient subtracts eight times its signed mean.

Refining cells, or splitting the complement, only increases this exact lower objective. This is the triangle inequality max(a+b,0)≤max(a,0)+max(b,0). For a fixed integrable W, expanding partitions and vanishing mesh recover ∫W_+ via L¹ approximation by cell averages. With a fixed finite basis, t_min>0, and a bounded objective, S≤J_partition bounds Σα_j≤J_partition/(2πt_min); finite-basis L¹ approximation then makes convergence uniform on that coefficient set. This explains an eventual certification strategy, not a guarantee that a given coarse partition is accurate.

Quadrature must be separated from this theorem. A numerical c_ij matrix produces a floating Jensen approximation. It has no rigorous lower-bound direction until coefficient errors are enclosed. If the finite-cell errors are bounded by |I_i−Ĩ_i|≤ε_i, and the complement is defined from the SAME exact S, then the full finite-plus-complement Jensen error is ≤Σ_iε_i. Proof: the complete coefficient vector has unchanged signed sum, J=(Σ|I|+S)/2, and its L¹ error is at most Σ_iε_i+|Σ_i error_i|≤2Σε_i. This gives a linear correction budget ε_i≤Σ_jα_jδ_ij for nonnegative weights.

## 4. Quantifying the positive-part gap for a retained candidate

A cell can contain positive and negative W even when its signed mean is zero. For example W(x,y)=x on [−1,1]×[0,1] has Jensen value zero and positive integral 1/2. Thus exact cell means alone cannot certify an upper cost. Point quadrature has neither direction: midpoint quadrature gives zero in that example, but on W=1−4x² over the same rectangle it gives 2 whereas the true positive integral is 2/3.

For a candidate and a finite cell B, suppose m≤W≤M. If m≥0 or M≤0 there is no Jensen gap. If m<0<M, convexity's chord gives

max(I_B,0)≤∫_B W_+≤M(I_B−m|B|)/(M−m).                   (7)

In particular the gap is at most |B|(M−m)/4. This follows by maximizing the chord gap over the average μ=I_B/|B|; its maximum is −mM/(M−m)≤(M−m)/4. A bound using the actual known mean in (7) is sharper.

Away from the origin, r_min=min_Bρ>0, the R55 derivative estimate gives

|∂_iW|≤A/r_min²+CΣ_jα_j/t_j²,
A=Σ_jα_j,  C=3√3/(8π).

Therefore the cell oscillation is at most its sum of coordinate widths times this bound. Cells meeting the origin use the exact upper ∫_B W_+≤A U_B instead, since all subtracted V terms are nonnegative. Directed box ranges can sharpen both estimates.

Outside a region containing {ρ≤R}, use the generalized R55 full-tail bound. With L the largest source coordinate absolute value and S₀ the largest source ρ, R>max(L,S₀),

∫_(ρ>R)W_+≤8AL log(R/(R−S₀))+8Σ_jα_jt_j²/(R−S₀).       (8)

The proof applies the compact cross-strip bound to the total atomic position measure and adds each positive Q tail. The finite region must truly contain the stated diamond; a cross-shaped region omitting corners does not contain arbitrarily large diamonds simply because its arms extend to a large cutoff.

D constraints also need their own full-domain check: sampled W≥1 and a sample minimum do not certify essinf_DW. Tiny negative rounded weights can invalidate C4 and the positive bound W_+≤AU. Keep the retained weights exactly nonnegative and verify the entire D independently.

## 5. Actual bounded tests and evidence limits

`computation/r060_cauchy_rectangle_jensen.py` actually ran with Python3.13.18, SymPy1.14.0, mpmath1.3.0 at 90 decimal digits, NumPy2.5.4 and SciPy1.18.1, one process, five deterministic rectangles, no random inputs. The mixed derivative of (4) was exactly checked. On all five rectangles, independent positive integrals (1)/(2) satisfied U_B−V_B−Q_B below 10⁻⁷⁵ in the stated scale; the 256-step monotone Riemann bracketing and its telescoping gap identity passed numerically. Gaussian/Radau moment checks and positive weights passed within 10⁻¹² residual.

Whole-interval n=16 Gaussian/Radau relative width reached about .497 for a far cell and .617 for a thin cell. The geometric split with n=8 per piece reduced relative widths across these five examples to about 1.4e−8…1.3e−7. These floating widths support the algorithm design, not directed enclosures.

Direct Q-cell fixed GL16 had relative error +1.2325905e−4 for [−8,−7.5]×[0,.01], t=1, and −2.1915099e−7 for [.999,1.001]×[−.001,.001], t=.03. Other three selected cases agreed near floating precision. Thus fixed GL16 has no universal lower direction and can miss endpoint transition scales. For proposal generation it can still be practical, with explicit floating status and GL16/32/refined checks on retained candidates. Breaks at nonzero endpoint distances in s can also resolve those transition scales.

Saved output is `computation/results/r060-cauchy-rectangle-jensen/checks.json`. All code/inputs are deterministic and embedded; reproduce from the topic directory with:

    /Users/zzh/.venvs/main/bin/python3 computation/r060_cauchy_rectangle_jensen.py

R55 analytic source read earlier and preserved: `reviews/cauchy-smoothed-potential-independent-review.md`, SHA256 `97acd0c6c955e6de40f8372cb242d72d3c14f7e8125bc49f8100879a38b084b2`. This R60 note extends its scalar potential and tail formulas to exact rectangle means and joint-basis optimization. It does not repeat or alter the frozen R55 run. No fixed-ray search was restarted.
