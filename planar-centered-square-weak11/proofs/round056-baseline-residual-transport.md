# R56: local baseline cancellation and the spatial residual transport interface

Author: Sol recursive geometry helper, 2026-11-15. Scope: a constructive
reformulation and exact small cases, not a new bound for C₂. Work in R² with
Lebesgue receivers, f≥0 in L¹, λ>0, Q_r(x)=x+[-r,r]², and the strict all-scale
halo H={x:∃r>0, ∫_{Q_r(x)}f>4λr²}. The two axes always use the same r.

## 1. Constructive cancellation before any box charging

Define, with densities relative to Lebesgue measure,

    b=min(f,λ)1_H,       q=(λ−f)_+1_H,
    σ₁=f−b=(f−λ)_+1_H+f1_{Hᶜ}.

Thus λ1_H=b+q and f=b+σ₁ pointwise a.e. Send b(x)dx directly from source x
to receiver x. This is an explicit transport and is allowed by every positive
radius centered witness at x. It uses source capacity at most f and leaves
residual demand q dx and residual source σ₁ dx. In particular q and σ₁ are the
negative and positive parts of the signed spatial field d=f−λ1_H.

For every measurable S (all expressions finite for a bounded witness square),

    σ₁(S)−q(S)=∫_S f−λ|S∩H|.                         (1)

If Q=Q_r(x) is a genuine strict witness, then

    σ₁(Q)−q(Q)=∫_Q f−λ|Q∩H|>λ|Q\H|≥0.           (2)

This is exact cancellation, not a replacement of one witness radius by separate
coordinate or source radii. It is valid for arbitrary L¹ inputs, not only atoms,
periodic inputs, finite boxes, or an asymptotic approximation.

For target source capacity Kf, K≥1, keep the same explicit baseline and leave

    σ_K=Kf−b=σ₁+(K−1)f.                             (3)

The second term is the retained old-mass reserve. Removing it and claiming a
capacity-one residual allocation would incorrectly prove C₂≤1.

## 2. The exact cut identity, with all source overlaps retained

Fix one Borel strict witness r(x) for each x∈H; a choice exists by testing
positive rational radii in a fixed enumeration. Strictness and L¹ continuity of
square integrals allow this rational restriction. Put Q_x=Q_{r(x)}(x). For a
Borel receiver set A⊂H define N(A)=∪_{x∈A}Q_x. This is an analytic
projection and therefore Lebesgue measurable. Finite certified receiver/source
partitions have no measurability issue.

Every residual allocation along y∈Q_x, with source capacity σ_K, necessarily
satisfies q(A)≤σ_K(N(A)). Since A⊂N(A), this is exactly

    λ|A| + b(N(A)\A) ≤ K∫_{N(A)}f.                  (4)

Source neighborhoods are unions, not sums: the same atom or density patch in
many capturing boxes is counted once on the right. The extra cross term is the
baseline already spent on receivers in N(A)\A. Formula (4), rather than only
(2) for single squares, is the remaining geometric obligation.

For any finite certified incidence model with disjoint receiver/source cells,
with receiver
demands q_j and source capacities s_k=σ_K(T_k), the residual allocation exists
iff ∑_{j∈J}q_j≤∑_{k∈N(J)}s_k for every receiver-index subset J. Necessity is
summation. For sufficiency, connect source→receiver cells along the certified
incidences, add demand and capacity edges, and maximize throughput on the compact finite
flow polytope. A nonzero unfilled demand with no augmenting path identifies a
violated cut. This existence proof does not require finite termination of naive
augmentations with irrational capacities. This is a precise finite verification interface; it is not a
claim that single-square positivity (2) proves all cuts, nor a proposed LP
substitute for the missing geometric bound. For general measurable incidence
relations we assert necessity only; extra closure/tightness is needed before
invoking an infinite transport existence theorem.

The forced identity baseline is a substantive restriction. An unconstrained
transport may move baseline mass away from its own receiver. Therefore (4) can
be stronger than the original unmodified Hall inequality λ|A|≤Kμ(N(A)); no
claim is made that this particular baseline attains the optimal weak constant.

## 3. Why the logarithmic background obstruction cancels exactly

Let D=[−L,L]×[0,L], A₀=[−1,1]×[0,1], L≥16, and
f=1_D+1_{A₀}. On the wedge

    W={(s,t):4≤t≤L/4, |s|≤t/2},

Q_t((s,t)) lies in D, contains A₀, and has mass 4t²+2. It is a strict genuine
common-radius witness. Thus W⊂H and f=λ=1 on W. Consequently b=1 and q=0 on W.
All demand of this wedge is paid directly, with no residual charge to any source.

Also S=[−1,1]×[1,2]⊂H: Q_t((s,t)), now 1≤t≤2 and |s|≤1, lies in D and
intersects A₀ in positive area, so its average is strictly >1. Hence σ₁=0 on S,
and σ_K=(K−1) on S. The proportional whole-witness construction repeatedly
charged this same patch; the canceled formulation cannot repeat that charge at
K=1, and explicitly labels its available reserve at K>1. This resolves that
obstruction to this interface only; it does not solve the residual cut problem.

## 4. Exact minimal cases and an L¹ obstruction to residual Hall-one

One atom mδ_p has H=p+(−sqrt(m/(4λ)),sqrt(m/(4λ)))² and λ|H|=m. Sending all
receivers to that atom has exact load one. The atomic baseline is zero because
its source measure is singular; the L¹ formula above is not applied to δ_p.

For a fully rational genuine L¹ two-source case, set λ=1, ε=1/100, p₀=(0,0),
p₁=(12/5,0), T_i=p_i+[-ε,ε]² and f=10000(1_{T₀}+1_{T₁}). Each source square
has mass4 and total mass is8. The following three disjoint open rectangles lie
in H (their boundary lines have zero area):

    R₀=(−49/50,49/50)²,
    R₁=(71/50,169/50)×(−49/50,49/50),
    R₂=(101/100,139/100)×(−139/100,139/100).

For R₀ and R₁ the same radius 99/100 captures its whole source square and has
average 10000/9801>1. For R₂ radius7/5 captures both whole source squares and
has average50/49>1. Hence

    |H|≥2·2401/625+2641/2500=21849/2500>8.         (5)

Both source supports are inside H, so ∫b=|T₀|+|T₁|=1/1250, and exactly

    ∫σ₁=8−1/1250=19998/2500,
    ∫q=|H|−1/1250≥21847/2500.

Even allowing every source to serve every receiver, residual capacity-one fails
by at least1849/2500; any residual inflation factor must be ≥21847/19998.
For the original full source capacity K, the total-demand cut already requires
K≥21849/20000. This is a small exact witness that (2) does not imply all cuts.
It is not competitive with the project's established lower bound.

## 5. What must be remembered

For capacity one the signed spatial field d=f−λ1_H determines residual demand
and supply via its two Jordan parts. At K>1 one additionally needs the spatial
old-mass reserve (K−1)f. Given H and d, this reconstructs f, so it is an exact
reparameterization, not a finite-dimensional compression.

For a box-gluing implementation one needs these fields restricted to every
source subregion queried by a common-radius square, or their full two-dimensional
rectangle cumulative profiles. Total deficit/excess per box and the two separate
coordinate marginals do not determine arbitrary rectangle intersections. The
R32 exact-collar result already supplies the latter marginal obstruction; it is
not rerun here. Eligible source boxes must remain actual capturing boxes for the
same receiver and radius; source capacities are charged on their spatial unions.

Atomic inputs can instead use the exact E_B receiver rectangles with
s_B=sqrt(μ(B)/(4λ)) and E_B=(b−s_B,a+s_B)×(d−s_B,c+s_B). Their finite endpoint
arrangement gives receiver cells with constant eligible-box lists. Do not merge
source labels shared by multiple boxes into independent budgets. The present
baseline theorem concerns L¹ inputs; the atomic arrangement is an exact separate
finite incidence interface, not an atomic substitution into f−λ1_H.

## 6. Actual exact checks and remaining question

The companion Python Fraction script checks all rational rectangle areas,
containments at 12 vertices, the strict averages, total residual imbalance, and
36 wedge witnesses, and records Python version and source SHA. These are finite
exact auxiliary checks; (1)–(4) are analytic proofs, not numerically certified
general inequalities. No MATLAB job was needed. No new general C₂ upper bound.
Open: prove a useful K bound for the overlapping-neighborhood cuts (4), or
construct a better allocation which can release part of the identity baseline.
