# R73: cap-preserving atomization and a local full-halo obstruction

2026-11-16, Asia/Shanghai. All halos use strict threshold1, centered open axis-parallel squares with one common radius. The target remains actual deletion D=H(sigma+nu)\H(nu), with area Delta. This proof does not bound Delta by7/2 or violate that bound. It shows that residual-halo credit cannot be dropped even after fixing a unique maximum uniform-grid block in the connected nonterminal finite domain.

## 1. Objects and the exact-cap construction

For an integer n>=2 set e=1/n², h=(n+1)/n², w=1−e, t=n/(n+1). Let S_n={0,h,...,(n−1)h}², with mass1/n² at each point, and write sigma_n for this mass1 measure. Set J=[−e,1]² and B=[−1,2]². The R72 shared Borel outside measure is nu_n=t 1_(B\J) dx.

Choose 0<delta<=1/(2n), and let

    J_delta=J+[-delta,delta]²,
    c_delta=t/(1+2delta).

Use the square lattice cells delta([i,i+1]×[j,j+1]). Take **every closed cell completely contained in B and disjoint from the closed set J_delta**. At its center place outside mass c_delta delta². The resulting finite atomic measure is nu_delta. Boundary-touching cells are excluded. No subset selection or enumeration of all atom subsets is part of this definition.

**Theorem.** Every closed square Q of side R satisfies

    (sigma_n+nu_delta)(Q)<=max(1,R²).                  (1)

Consequently S_n is a maximum eligible block of mass1 for every delta. This is an exact statement at each mesh, stronger than merely having maximum mass tend to1. For any 0<lambda_delta<1, replacing nu_delta by lambda_delta nu_delta makes S_n the unique maximum block by the R72 outside-contraction theorem.

## 2. Enlarging cells without losing the square caps

Let Q^+ be the concentric square of side R+2delta. Every selected cell whose center belongs to Q is contained in Q^+. Since the selected cells are pairwise interior-disjoint and avoid J_delta,

    nu_delta(Q)<=c_delta |Q^+\J_delta|.               (2)

Suppose R>=1 and Q captures k grid coordinates on one axis and l on the other. If k,l>=1, put

    a_k=(k−1)h+e=[k(n+1)−n]/n²,
    a_l=(l−1)h+e.

The R72 interval-overlap lemma gives |I∩[−e,1]|>=a_k for an interval I of length R capturing k grid coordinates. Both intervals have nonempty intersection. Expanding both by delta commutes with their intersection in one dimension: the maximum of lower endpoints and the minimum of upper endpoints each move by delta. Thus the corresponding expanded overlap is at least a_k+2delta. The two coordinates are independent, giving

    (sigma_n+nu_delta)(Q)
      <=kl/n²+c_delta[(R+2delta)²−(a_k+2delta)(a_l+2delta)].  (3)

This bound uses the same actual outside measure for every Q. The cap argument uses no shape or size property of B beyond finiteness of the selected cells: the same construction and all caps hold with any bounded truncation region in place of B. The weak-limit and numerical statements below retain their stated B=[−1,2]² unless another region is explicitly specified.

For clarity the required cancellation can be written as an exact factorization. Let G=(n−k)(n−l)/n³ and c=c_delta. The excess of the right side of (3) over R² equals

    −(1−2n delta)G/(1+2delta)
    −(R−1)[(1−c)(R+1)−4c delta].                       (4)

The first term is nonpositive, and

    (1−c)(R+1)−4c delta
       =(1−c)(R−1)+2(1−t)>0.                          (5)

Thus (3) is at most R². One can also see the compensation through P=(1−a_k)(1−a_l):

    (R+2delta)²−(a_k+2delta)(a_l+2delta)
       =(1+2delta)(R²−a_k a_l)+2delta[P−(R−1)²],
    1−t−[kl/n²−t a_k a_l]=nP/(n+1)².

The density reduction exactly pays the expansion term whenever2n delta<=1.

If k=0 or l=0, the internal mass vanishes. By (2),

    nu_delta(Q)<=c(R+2delta)²
       <=c(1+2delta)²R²=t(1+2delta)R²<=R².

Finally any square of side R<=1 lies in a closed unit square, already controlled by the preceding argument. This gives mass<=1, not mass<=R², and proves(1) in all cases, including closed boundary ties.

## 3. Exact admissible atomization still fails to transfer deletion area

For fixed n, as delta tends to0, the complete selected cells exhaust B\J up to a boundary layer of area tending to0. Replacing uniform mass on each selected cell by its center changes integrals against compactly supported continuous functions by their modulus of continuity at distance at most delta. Since c_delta tends to t, nu_delta converges weakly to nu_n. The same holds after any common contraction0<lambda_delta<1 tending to1.

Every outside singleton has strict halo a square of side sqrt(lambda_delta c_delta)delta. It lies in its own cell because lambda_delta c_delta<1; these halo squares are pairwise disjoint. Let C=[−1/2,w+1/2]². Counting cells completely inside C and outside J_delta gives

    liminf |H(lambda_delta nu_delta)∩(C\J)|
       >=t|C\J|=t(3−6e).                             (6)

The omitted cell boundary layers have vanishing area. Unlike R72's unrepaired atomization, the maximum mass here is exactly1, so the actual deletion is contained in C at each mesh. Subtracting (6) yields

    limsup Delta_(sigma_n+lambda_delta nu_delta)(S_n)
       <= |C|−t|C\J|
       =(1+e)²+(3−6e)/(n+1)=U_n.                      (7)

For n=5 this is1927/1250=1.5416, whereas the R72 rigorous Borel-family lower bound is27886815467/10485760000>2.659. Hence deletion-area discontinuity occurs **within the exact maximum-mass-one finite class**, and even with S_n uniquely maximizing at every mesh. This concerns the difference of full and residual halos; it does not contradict the single-full-halo Fatou transfer used for general weak bounds.

The limit n->infinity in U_n is taken only after the fixed-n mesh limit. No conclusion about arbitrary finite outside measures follows from this particular approximating sequence.

## 4. The collar lies in the internal halo for n>=5

We record the independent R73 collar observation with its explicit witnesses. Around each grid coordinate jh take the closed band[jh−1/(2n),jh+1/(2n)]. Their intervening gaps have width e. The endpoint collars[−e,0] and[w,1] are within the first/last bands for n>=2.

If both coordinates of a point of J are in bands, the associated singleton square covers it except possibly band endpoints. At such endpoints one may instead use the adjacent pair or four-point square described below; alternatively regard the closed band boundaries as belonging to the gap category.

If one coordinate is in a gap and the other in a band, use the two adjacent sources in the gap direction. Their mass is2/n². Their maximum source distance is at most max((n+2)/(2n²),1/(2n)), strictly less than sqrt(2)/(2n) for n>=5, since n+2<sqrt(2)n. This gives a strict open witness. If both coordinates are in gaps, use the four grid corners: maximum distance at most(n+2)/(2n²)<1/n for n>2, while the mass is4/n². Points strictly inside both bands have a singleton strict witness. Thus

    J subset H(sigma_n),  n>=5.                       (8)

This proof also explains the small-n qualification; we do not assert(8) for n=3,4. Notice that nu(J)=0 does not imply H(nu)∩J is empty.

## 5. A finite exact input disproves a local full-halo budget

Equations(6),(8) imply

    liminf |H(sigma_n+lambda_delta nu_delta)∩C|
       >= |J|+t|C\J|.

For every n>=7 the right side is strictly greater than7/2. This already shows that a proposed bound on the local **full** halo cannot replace actual deletion. The following explicit finite input avoids relying only on existence through a limit.

Set

    n=9, delta=1/100, c=15/17, lambda=9999/10000.

Take the above full cell construction, then multiply every outside mass by lambda. The actual outside mass per cell is29997/340000000. Complete cells inside B have each axis index−100,...,199. Cells intersecting the closed J_delta have each axis index−3,...,101, so there are

    300²−105²=78975 outside atoms,
    78975+81=79056 total positive sources.

These are counts from the explicit rectangular index formula, not an executed enumeration of79056 sources. Every axis index−50,...,147 gives a complete cell inside C. Excluding the expanded collar leaves198²−105²=28179 cells whose singleton halos are inside C\J. Their halos are pairwise disjoint and disjoint from J; all of J is internally covered by(8). Therefore

    |H(sigma_9+lambda nu_delta)∩C|
       >=(82/81)²+28179·29997/340000000
       =7832077922743/2230740000000
       >7/2,                                         (9)

with exact excess24487922743/2230740000000, about0.01097749. This is a lower bound, not a computation of the full halo's exact area. It says nothing comparable about the actual deletion, because the large singleton contribution in(9) is residual halo and must be subtracted.

The example meets further structural qualifications. Exact caps(1) hold before contraction, and the R72 theorem makes S_9 unique after contraction. The outside grid is a rectangular frame with a rectangular hole strictly inside it; it is connected by horizontal and vertical nearest-neighbor edges. Each edge is eligible since its pair mass exceeds delta². Adjacent internal sources are eligible since2/81>(10/81)². The source at(0,0) connects to the outside cell center(−7/200,1/200), because its pair mass is greater than49/40000. These eligible edges establish connectedness of the entire actual input without enumerating all eligible subsets.

The extreme grid centers give actual span299/100. The total mass and nonterminal gap are

    m(A)=108360523/13600000,
    d_A²−(35/33)m(A)=8787943/17952000>0.                (10)

Thus a unique maximum uniform S_9, connectedness and nonterminal status do not rescue the local-full-halo shortcut.

## 6. Relation to a unit-window integral certificate

For any nonnegative integrable function f of unit-square centers z, write

    K_f(q)=integral f(z) 1_(|z−q|_infinity<=1/2) dz.

Integrating the same actual unit-square capacity gives

    integral K_f dnu <= integral f −integral K_f dsigma_n.

For f=1_C the continuous bounded kernel is a product of one-dimensional tents, with plateau1 on[0,w] and support[−1,w+1]. Its value is1 at every internal source, so the right side is|C|−1. The independent R73 unit-window proof evaluates the slack for the Borel family as

    s_n=|C|−1−t[|C|−(1+e−e²)²] ->0.

Since K_f is continuous and bounded, the exact-cap atomic measures above have slack tending to s_n for fixed n. A diagonal choice of meshes and contractions therefore gives finite atomic unique-maximum inputs with this averaged-cap slack tending to0 as n tends to infinity. No n-independent positive budget can be extracted from **the slack of this one uniform unit-window integral alone**. This does not rule out other weighted window certificates or inequalities using actual residual credit. The large-n sequence is not asserted nonterminal; the explicit n=9 obstruction(9) is nonterminal.

## 7. Actual computational scope

The accompanying standard-library Fraction program checks3092 positive k,l parameter pairs,9276 radius identity samples and42 zero-capture samples across n=2,3,5,7,9,17,33 and two admissible meshes each. These are exact sample checks of the algebra; formulas(2)–(5) prove all parameters. The same run computes the explicit n=9 rectangular counts, mass, connectivity margins, nonterminal gap and lower area(9) exactly. It does not enumerate the79056 actual sources, their source subsets, or any complete halo catalogue. Its result and process receipt record the executed program hash, Python version, normal exit and process absence.

The main remaining task is still a bound on Hfull minus Hresidual, with the same outside measure used throughout the edge and corner accounting. Neither the stronger exterior-gap conjecture nor the corner-only actual-deletion conjecture is settled here.
