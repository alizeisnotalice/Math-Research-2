# R40 covering-core review: sharp insertion reserve and ordered-envelope obstruction

Date: 2026-11-15 (Asia/Shanghai). Main repository and old T read only; no commit/push. This review concerns an auxiliary covering envelope, not an evaluated C₂.

## 1. Sharp nonvanishing old-mass reserve

For any finite nonnegative Borel measure μ of mass M and ν=εδ_p, ε>0, let A=H(μ+ν)\(Hμ∪Hν), with the fixed strict threshold 1 and all common-radius centered squares. The new statement in `proofs/point-insertion-shell-budget.md` is analytically verified:

    A⊂S_(M+ε)(p)\S_ε(p),   |A|≤M,
    S_t(p)={x:||x−p||∞<√t/2}.

Every witness at x∈A must capture p; otherwise x∈Hμ. Hence ||x−p||∞≤r<√(M+ε)/2. Removing the point halo S_ε yields a square shell of area M. For μ=Mδ_p and ε≥M, equality holds. No atomicity or bounded support is needed for μ. This is an actual repair of the R34 vanishing-budget failure, but it spends old mass and cannot be summed naively.

The explicit overpayment examples were checked. N distinct equal-mass points in a segment of width 1/10 have Σ(previous total mass)=(N−1)/2 under every order, while their final halo lies inside a rectangle of area at most 11/10. Also μ=δ_(2,0), ν=δ_0 has no mixed witness because a square capturing both has area at least4>2; an unrestricted total-mass shell still charges1. A valid joint account must retain both old coverage and common-capture geometry.

The auxiliary claim |Hμ|≥μ(R²) is false, as the two unit atoms at (±1/10,0) have complete halo equal to their two-source rectangle, of area 2−√2/5<2. The containment of both singleton halos follows from √2>7/5. Therefore source capacity P=M cannot be converted into a receiver-area lower bound this way.

## 2. Effective-mass envelope and its implication for a general upper bound

For a finite atomic input and a chosen insertion order, let μ_j be the old prefix and m_jδ_(p_j) the inserted atom. For each nonempty old source subset S define d(S,p_j) as the largest coordinate span of its bounding box together with p_j. Set

    E_j=max{m(S):d(S,p_j)²<m(S)+m_j},

with E_j=0 if there is no eligible nonempty subset. The candidate envelope square is Q_j=S_(E_j+m_j)(p_j).

For a newly activated center outside Hμ_j, its witness contains p_j and an old captured subset S, unless it already belongs to the inserted singleton halo. In the former case d(S,p_j)≤2r<√(m(S)+m_j), hence m(S)≤E_j, and the center lies in Q_j. The singleton halo is also contained in Q_j. Thus

    H(μ_j+m_jδ_(p_j))⊂Hμ_j∪Q_j,
    H(μ_final)⊂U_π:=∪_j Q_j.

Filling any candidate source box with all old atoms it contains preserves its geometry and increases its mass, so a maximizing set can be taken to be a filled coordinate box. This gives a finite geometric construction; it does not refer to the original halo-area supremum.

A uniform bound min_π|U_π|≤K∑m_j for every finite atomic input would imply the general L¹ weak bound with coefficient K. To check the transfer explicitly, first truncate a nonnegative L¹ f to a bounded square. Replace each grid cell's mass by an atom at its center. Every source point moves by at most δ in ∞ distance, so for any fixed strict witness Q(x,r), the atomic measure captures at least its original f mass in Q(x,r+δ). For sufficiently fine cells the witness remains strict. Hence H(f) is contained in the pointwise liminf of the atomic halos, and Fatou gives the same K bound; truncation is removed monotonically. This is a valid implication, not a claim that any proposed K has been proved.

The envelope loses geometry because Q_j is centered at the newly inserted source rather than retaining the source-box intersection rectangle. Four and six source examples demonstrate this loss exactly. The four-rectangle family has every-order value 13/8+d/2, tending to 13/8+√2/4≈1.97855339. The six-source all-order minimum is 5/3−√2/6+(24/25)/√6≈1.82288277. These are envelope lower bounds only.

## 3. Eight-source all-order certificate: analytic scope

The exact input consists of the eight points in {−49/100,0,49/100}² other than (0,0), each of mass1/8. The source order is the lexicographic point order

    0=(-a,-a), 1=(-a,0), 2=(-a,a), 3=(0,-a),
    4=(0,a), 5=(a,-a), 6=(a,0), 7=(a,a),  a=49/100.

D4 symmetries preserve the input, squared coordinate-span eligibility, and union area. Every permutation can therefore be mapped to one ending at corner0 or edge-midpoint1; all 7! permutations before each fixed endpoint suffice. This is a complete reduction of the 8! order problem, not a random sample.

The producer enumerated all10,080 representative orders and obtained348 distinct area-type signatures. Each signature records the square area numerator k in k/8 at every vertex. The effective choices were calculated using exact source geometry. For b=32, every half-side r=√(k/32) has the rigorous bracket

    floor(2^b√(k/32))/2^b ≤r≤ceil(2^b√(k/32))/2^b.

Both numerators come from integer `isqrt(k·2^(2b−5))`. Replacing each square by the inward bracket gives an exact lower bound on its union; replacing it by the outward bracket gives an upper bound. All coordinates are integers after multiplication by100·2^b. Therefore no floating-point area is used in the certificate.

The complete envelope minimum lies between

    378116406113191688339 / 180143985094819840000

and

    12099725004195137798823 / 5764607523034234880000.

The lower bound exceeds2 by exactly

    17828435923552008339 / 180143985094819840000 >0.

Thus every insertion order has envelope area greater than2, for a true positive rational eight-source input of mass1. This strictly disproves the proposed universal envelope bound2. It does not prove C₂>2, because Hμ is only a subset of U_π.

## 4. Exact radical value of the minimum

Signature237 is the unique possible minimizer among the348 signatures: every other signature has inward area at least

    777358489775631255821 / 360287970189639680000 >2.1576,

which is larger than signature237's outward area. A representative order is [0,2,5,3,7,4,6,1]. Its square area numerators are [1,8,1,2,2,1,3,1]. The four corners have half-side r=√2/8, the final left edge-midpoint has half-side1/2, the top/bottom midpoint squares have half-side1/4, and the right edge-midpoint has half-side u=√6/8.

There are exactly six positive-area pair intersections: (0,1),(1,2),(1,3),(1,4),(3,6),(4,6), and no positive-area triple intersection. The two global/pair intersections each have area(13/50)²; the two global/corner intersections sum to4r(r+1/100); the two right-triad/pair intersections each have area(u−6/25)². All other intersections are absent; in particular a>r+u excludes the right triad from the adjacent corner squares. Inclusion–exclusion therefore gives

    min_π |U_π|
      =19/8−2(13/50)²−4r(r+1/100)−2(u−6/25)²
      =18121/10000+3√6/25−√2/200
      ≈2.098967701322116.

The decimal is only a display. The exact root formula and the full-order integer certificate give the strict conclusion.

## 5. Independent actual verification

I inspected the producer's definitions and DP direction. Independently, `computation/review_eight_point_order_proxy.py` enumerates every submask directly (17,496 tests), rather than using its subset DP; it reconstructs every recorded order, checks all10,080 distinct orders and both factorial cardinalities, and verifies each signature. For geometry it uses the full two-dimensional endpoint cell decomposition, rather than the producer's one-dimensional sweep over interval unions. It recalculates696 exact inner/outer union areas for all348 signatures.

The first complete independent run passed in about0.135 seconds. After adding the unique-minimizer and pair-overlap checks, a fresh full run passed in about0.110 seconds, exit0. After the author files were staged byte-for-byte for integration, the explicit lexicographic vertex-order assertion was added and the full audit was run once more against the staged files: about0.084 seconds, exit0. This third run is the archived final result `computation/results/round040/eight-point-independent-review.json`. It verifies signature237's six pair intersections already have positive inward area, no additional pair has positive outward area, and every outward triple intersection has zero area. This independently supports the exact radical simplification in §4. Python3.15.0, standard-library integers and Fraction only, one process, no random seed or tolerance. Repeated runs are not counted as different inputs.

Reproduce after integration from the topic directory:

    python3 -B computation/review_eight_point_order_proxy.py --source-dir computation/results/eight-point-order-proxy --out computation/results/round040/eight-point-independent-review.json

## 6. Precise remaining issue

The sharp single-insertion reserve and effective-mass envelope are valid. The envelope coefficients7/4 and2 are false. The whole interface is not ruled out: no proof or counterexample here settles a general envelope bound K<3.590865312933. The available easy inclusion U_π⊂{Mμ>1/9} yields only9C₂ times total mass and does not improve the existing weak bound. Its pointwise factor is asymptotically sharp for a fixed unfavorable order, while a different order can substantially reduce that example.

Two distinct avenues remain: a global ordering/overlap theorem with coefficient above the exact eight-point obstruction but below the inherited upper bound, or a covering allocation that retains the source-box receiver rectangles instead of replacing them by full source-centered squares. Source-load constraints without area coverage are still insufficient. No general upper-bound impossibility, new C₂ lower bound, or exact C₂ value is claimed.
