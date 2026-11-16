# R79: two-sided local profiles for general maximum blocks, and a strict transfer obstruction

2026-11-16. Strict threshold1, centered axis-parallel squares with one common radius. The source measures are finite positive atomic. This note separates a general safe state interface from the unresolved area payment. It does not infer a general weak bound from a uniform-grid calculation.

## 1. The general deletion problem retains arbitrary internal mass profiles

Normalize a maximum eligible block S to mass1, with internal measure sigma and bbox [0,w]x[0,h], w,h<1. Let nu be the original common outside measure. The R63/R64 maximum-block criterion gives, for every closed axis-parallel rectangle B of maximum side length d_B,

    sigma(B)+nu(B)<=max(1,d_B²).                     (1)

This follows by enclosing B in a square of side d_B. All actual strict full and residual radii are below1/2. The deletion D=H(sigma+nu)\H(nu) lies in the R63 container

    I_S=(-1/2,w+1/2)x(-1/2,h+1/2).

The empty source collar and spatial localization from R63 remain exact constraints. The hard connected nonterminal branch still allows arbitrary numbers of atoms, arbitrary positive weights, and concentration inside a nearly square S. None of these facts supplies the grid spacing or the near-density estimate used in the separate R79 quadratic proof.

At a southwest receiver x=(-X,-Y), X,Y>0, any square capturing an old source has r>max(X,Y). Its old mass is exactly

    F_sigma(r-X,r-Y),
    F_sigma(u,v)=sigma({p_x<u,p_y<v}).                (2)

Thus an anchored cumulative profile is a natural state in this corner, even if S has no atom at its bbox corner. The radius is common in both arguments. Other parts of I_S require two-sided rectangle queries rather than just this anchored profile. Section2 gives a boundary-correct way to carry that information.

## 2. Four boundary conventions and safe rectangle intervals

For epsilon,eta in {<,<=}, define

    F^(epsilon,eta)(u,v)=sigma({p_x epsilon u,p_y eta v}).

Atomic boundary conventions matter. For an open rectangle Q=(l,u)x(b,v) with l<u and b<v, the exact identity is

    sigma(Q)=F^(<,<)(u,v)-F^(<=,<)(l,v)
               -F^(<,<=)(u,b)+F^(<=,<=)(l,b).        (3)

For the closed rectangle B=[l,u]x[b,v] with l<=u and b<=v, it is

    sigma(B)=F^(<=,<=)(u,v)-F^(<,<=)(l,v)
               -F^(<=,<)(u,b)+F^(<,<)(l,b).          (4)

These identities follow by multiplying the two one-dimensional open or closed interval indicators; they allow boundary atoms under the stated endpoint order. Zero- or negative-width open rectangles are empty and are assigned mass0 directly, not evaluated using(3). Closed rectangles may have zero width and then carry boundary mass; reversed closed intervals are empty and handled separately. Half-open variants follow the corresponding interval conventions, with empty intervals handled directly. In particular(3) is not extended to l=u: that difference of indicators would produce a negative boundary atom instead of the empty open interval.

Suppose every queried quadrant value has verified lower and upper bounds L_i<=F_i<=U_i. Write either(3) or(4) as sum_i s_i F_i, with two signs +1 and two signs -1. Safe rectangle bounds are

    ell(B)=max(0, sum_(s_i=+1)L_i-sum_(s_i=-1)U_i),
    u(B)=min(1, sum_(s_i=+1)U_i-sum_(s_i=-1)L_i).     (5)

They satisfy ell(B)<=sigma(B)<=u(B) for the original sigma. The same notation with the correct endpoint conventions is used for open squares. Subtracting upper quadrant bounds in an alleged rectangle upper bound is unsafe; the negative terms require lower bounds.

All four profiles come from one measure. Exact profiles are monotone, have the corresponding boundary limits, and have nonnegative rectangular increments. Arbitrary endpoint intervals need not admit a common realizing sigma, even when each separate interval looks plausible. For a safe upper relaxation it is permissible to allow additional unrealizable states, provided the actual state is contained. Such a relaxed feasible point is not a new source input or a lower certificate. A claimed input still requires one positive measure realizing all coordinates, weights, endpoint conventions and normalization together.

On a fixed finite endpoint arrangement, nonnegative masses on its full point/interval product strata can enforce common realizability. This is the existing finite arrangement principle, not a bounded-size discretization of all general blocks. The number of required queries and the admissible mass error are not uniformly controlled by this observation.

## 3. Full witnesses use an upper old bound; capacity rows use a lower old bound

Assume verified rectangle bounds ell<=sigma<=u have been obtained by(5), or by another valid method. For an actual full witness W=Q_open(x,r),

    nu(W)>4r²-sigma(W)>=4r²-u(W).                   (6)

This is a necessary upper-envelope condition. In contrast, the original capacity(1) implies only

    nu(B)<=max(1,d_B²)-sigma(B)
          <=max(1,d_B²)-ell(B).                     (7)

Therefore a safe outside-measure relaxation uses the OLD LOWER bound in capacity right sides. Replacing sigma(B) by its upper bound in(1) would strengthen the row and could discard a legal original outside measure. This is not a harmless approximation in a dual infeasibility argument. The R77 unit-window example, whose original mass339/400 becomes409/400 under the synthetic replacement, supplies an explicit already-qualified warning about that direction.

Retain the complete true residual H(nu), with the same nu at all centers and radii. Define

    O_u(nu)={x in I_S: some 0<r<1/2 has
                         nu(Q_open(x,r))+u(Q_open(x,r))>4r²}.

Every actual deletion point satisfies

    D subset O_u(nu)\H(nu).                        (8)

This interface applies to arbitrary many outside atoms; no two-atom assumption is needed. It may be extremely loose. In particular its outside mass in the receiving square is still the actual shared mass, not an independently selected amount at each receiver. Exact total mass, the old bbox and the zero outside collar should be retained as separate known constraints even if approximate profile bounds would otherwise weaken them.

## 4. A quantitative shared-cut direction, not a new finite LP theorem

The R71 geometric weak-duality certificate can use these interval rows. To make the direction explicit, choose witness squares W_j, upper capacity rectangles B_a, and residual squares V_b at centers assumed to lie outside H(nu). Their radii and coefficients are fixed. Suppose nonnegative coefficients lambda_j,alpha_a,beta_b obey a pointwise cover on the allowed outside-source domain,

    sum_j lambda_j 1_(W_j)
       <=sum_a alpha_a 1_(B_a)+sum_b beta_b 1_(V_b),
    sum_j lambda_j>0.                               (9)

Every hypothesized common outside measure then obeys the integrated upper bound

    sum_j lambda_j nu(W_j)
       <=sum_a alpha_a[max(1,d_(B_a)²)-ell(B_a)]
          +sum_b beta_b area(V_b).                  (10)

Residual square inequalities retain all original radii and closed boundary atoms, as in R71. Simultaneous strict full witnesses would force the left side strictly above sum_j lambda_j[4r_j²-u(W_j)]. Thus if the right side of(10) is no greater than that last quantity, the proposed common witnesses and residual exclusions are impossible.

For a fixed original sigma, the loss in the separating margin relative to exact old masses is exactly

    sum_j lambda_j[u(W_j)-sigma(W_j)]
       +sum_a alpha_a[sigma(B_a)-ell(B_a)]>=0.       (11)

Both errors have the same adverse sign. If every relevant error is at most epsilon, the loss is at most epsilon(sum lambda+sum alpha); it is not at most epsilon without a coefficient bound. Formula(11) identifies a concrete stability requirement for a compressed profile state. It does not assert that useful certificates with controlled weights exist.

This is a directional substitution into the already-established R71 shared-measure cut, rather than a replacement for that theory. It must be checked over all allowed source positions, including boundary strata. A finite candidate-column dual is not automatically such a cover. Nor does excluding finitely many selected receiver witnesses bound the area of the full continuum. A general area certificate must cover all relevant receiver/radius states or supply a valid integral analogue, while retaining the same nu.

## 5. Compatibility with a deletion tree still requires a universal area budget

The four exact cumulative profiles are additive for disjoint source measures: F_(sigma1+sigma2)=F_sigma1+F_sigma2 with each boundary convention unchanged. Accordingly, sums of verified lower and upper states safely bound a union of disjoint source blocks. Residual subtraction reverses the removed block's intervals: from ell_A<=F_A<=u_A and ell_S<=F_S<=u_S one obtains ell_A-u_S<=F_(A-S)<=u_A-ell_S, with nonnegativity imposed. Repeated interval operations can widen uncertainty; additivity by itself gives no stable area payment.

For a covering construction with overlapping source labels, even this simple addition must be replaced by the actual incidence multiplicities or a justified nonnegative cover as in(9). The R64 fact that actual strict source states are eligible unions preserves common-radius geometry; it does not authorize arbitrary hypergraph or independently optimized receiver states.

A sufficient route back to R63 is now precise. For every connected nonterminal finite input encountered by the deletion induction, select an actual maximum block S and construct a state containing its true sigma. Prove, uniformly for every common outside nu in the corresponding SAFE capacity relaxation, the area bound

    |O_u(nu)\H(nu)| <=(7/2)m(S),                   (12)

after scaling back from mass1. Then(8) gives the true block deletion budget, and the R63 disjoint-deletion telescoping, together with its inherited terminal-component payment, yields the conditional general weak bound. One can instead certify a smaller necessary envelope than O_u, but it must still contain the true deletion for every input in the stated domain.

What is missing is(12), including an admissible general profile construction and a receiver-area payment. A bound restricted to uniform S_n with two effective outside positions does not supply it: maximum blocks in general inputs need not be uniform, and arbitrary remaining sources need not have bounded corner multiplicity. An O(1/n) error in a grid's captured mass also does not bound the area of near-threshold receivers. Neither subtraction of residual halos nor an overlapping covering tree can be passed through such an approximation without an additional stability or direct integral argument.

## 6. A five-atom obstruction within the hard general-block domain

The distinction from uniform grids already occurs in the following exact finite input. Take old block S at

    (0,0), (9/10,0), (0,9/10), (9/10,9/10)

with masses97/100,1/100,1/100,1/100. Add the single outside atom

    q=(-3/10,0), mass1/200.                          (13)

The old block has mass1 and span9/10, hence is eligible. Any proper subset of S together with the outside atom has mass at most199/200<1. The entire five-source set has mass201/200 but span squared36/25>201/200, so it is not eligible. Therefore S is the unique maximum eligible block; all square caps follow from the exact maximum-block criterion.

Every small old corner atom connects to the heavy origin by an eligible pair, since98/100>(9/10)²; the outside atom connects to the origin since39/40>(3/10)². Hence the source input is eligible-connected. It is proper and nonterminal, with

    36/25-(35/33)(201/200)=823/2200>0.               (14)

Both bbox side lengths are9/10>3/4, so it is not removed by the elementary thin-bbox branch of R63.

For every0<u,v<=9/10, the strict anchored old cumulative profile is

    F_sigma(u,v)=97/100.                            (15)

Consequently no estimate F_sigma(u,v)<=Cuv with a fixed finite C can hold for this block at all small positive u,v. An upper model C(u+eta)(v+eta) with fixed C requires C eta²>=97/100 by taking u,v down to zero. A density1 model uv+epsilon requires epsilon>=97/100. Thus a mesh-error term tending to zero, of the uniform-grid type, is not implied by full maximum-block capacity, fat bbox, uniqueness, connectivity or nonterminal status. This does not preclude an adapted profile recording the heavy atom explicitly.

There is also a genuine regional obstruction to equal per-corner allocation. Let zeta=1-w=1/10 and take the southwest corner

    Gamma_S=(-1/2,-zeta)².

The residual singleton halo has vertical extent only sqrt(1/200)/2<1/10, so it misses this entire corner. The heavy old singleton already yields deletion area at least

    (sqrt(97/100)/2-1/10)²
       =(101-4sqrt(97))/400>61/400.                 (16)

For a square bbox of width w, the noncorner container area is3-3zeta². If one were to divide the remaining7/2 budget equally between its four corner regions, the quota would be

    [7/2-(3-3zeta²)]/4=1/8+3zeta²/4=53/400.        (17)

Thus this legal input's actual southwest deletion strictly exceeds that equal quota. This does not disprove a four-corner total budget: the other corners may cost much less. It proves that an equal corner template inherited from uniform symmetry is not a general-block theorem.

A rational whole-cell certificate avoids any need to integrate a full halo. Set

    E=[-49/100,-101/1000]², radius491/1000.

Every center in E has the heavy origin strictly inside this same open square. Its old label mass97/100 exceeds4r² by1419/250000. The outside atom is also captured, only increasing full mass. No residual witness exists because its minimum required radius is at least101/1000, whereas1/200<4(101/1000)². Hence E is contained in true deletion and

    |E|=151321/1000000>53/400,
    |E|-53/400=18821/1000000.                        (18)

The source fixture establishes profile concentration and unequal corner cost, not a violation of the complete7/2 deletion target or a new lower bound for C2.

## 7. Actual exact kernel and scope

`computation/round079/check_general_profile_obstruction.py` was run once with Python3.15.0/Fraction. It checks all31 nonempty subsets of the fixed five-source input, the resulting14 distinct bbox closures and their caps, all14 eligible subsets, unique maximum mask15, connectivity and(14). It verifies the whole-cell witness/residual exclusion in(18). It also checks200 open/closed rectangle identities with atomic boundary ties and their interval signs in(3)-(5). Its remaining eight n-value rows verify the scalar parity formulas of the separate thin-joint-band note; they are not a search over input measures.

PID17722 exited0 and was confirmed absent. Kernel time was about.005 seconds. No full halo union, producer tile computation, general profile optimization or random scan is part of this run. The proofs of the all-state directions and the five-source qualifications above are analytic; the finite kernel gives explicit reproducible verification of their constants and boundary accounting.
