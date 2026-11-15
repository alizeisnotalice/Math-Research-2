# R67: positive support padding does not repair the three-quarter envelope

2026-11-15. The R66 three-source counterexample used a maximum block of rank2, which the elementary support-rank bound pays directly. This note gives an entirely positive rank9 maximum block for which replacing the filled low bbox by the union of equal-radius source neighborhoods makes **no geometric change**. It refutes that proposed theta3/4 repair; it does not claim that every fixed theta or every mass-sensitive neighborhood construction fails.

## 1. Positive input and unique maximum block

Put

    d=1023/1024,  eta=2^(-20),  epsilon=1/1024,
    G={0,d/2,d}×{0,d/2,d}.

Let S be the nine grid sources in G. The origin has mass1/5, the corner(d,d) has mass h=4/5−7eta, and the other seven grid sources each have mass eta. Add one outside source

    q=(d+1/8,d),  m(q)=epsilon.

All ten masses are positive, m(S)=1 and m(A)=1+epsilon. S has both bbox widths d∈(3/4,1), hence is eligible and fat. Every proper subset of S has mass<1. A subset containing q but not the origin has mass at most4/5+epsilon<1. A subset containing both q and the origin has largest span d+1/8, whose square exceeds the mass of the entire input. Thus none of these crossing subsets can be eligible with mass≥1. Therefore L=1 and S is the unique maximum block, now of rank9.

The full S is an eligible edge, and the pair{(d,d),q} is eligible because its span squared1/64 is below h+epsilon. Thus all ten sources form a single eligible component. Its nonterminal margin is unchanged from the three-source example:

    d(A)²−(35/33)m(A)=6982433/34603008>0.                (1)

No zero-mass support point is used to enlarge a bbox or support neighborhood.

## 2. Equal-radius support neighborhoods still fill the old low container

At theta3/4 let r=sqrt3/4 and define the support-based low container

    Z=union_(p∈S)(p+(-r,r)²).

The source grid is Cartesian. Consequently Z is the product of the two one-dimensional unions of intervals centered at0,d/2,d. Consecutive centers are distance d/2 apart, with

    (d/2)²=1046529/4194304<3/4=(2r)².

The three open intervals overlap, and their union is(−r,d+r). Hence exactly

    Z=(−r,d+r)²=C_(3/4).                               (2)

This is a statement about the support container, independent of how small the seven positive padding masses are. Its use of the common maximal mass L gives those tiny sources the same low-neighborhood radius as the heavy source.

## 3. Strict budget failure survives the padding

The residual halo after removing S is the singleton square about q with half-side1/64. It lies inside Z, so its old-halo credit is exactly epsilon. The heavy singleton is an eligible closed high label because

    h=4194269/5242880>3/4.

Its receiver square has half-side sqrt(h)/2. As in R66, its part outside Z has area

    phi(h)=3h/4−3/16−sqrt(3h)/4.                        (3)

Since the residual halo is inside Z, this entire part is unpaid high area. Thus the exact envelope obtained with the support container Z satisfies

    B_Z≥(d+sqrt3/2)²−epsilon+phi(h).                    (4)

For t between h and4/5,

    phi'(t)=3/4−sqrt3/(8sqrt(t))∈[0,3/4].

Therefore phi(h)≥phi(4/5)−(3/4)(7eta). The already elementary R66 rational lower certificate, or direct substitution of sqrt3>433/250 and sqrt(3/20)<97/250, now gives

    B_Z−7/2 > 246397/131072000−21eta/4
               =982963/524288000>0.                    (5)

Thus the support-neighborhood repair fails in a connected nonterminal input whose unique maximum block is fat and has rank9. The additional low-rank reduction does not remove this counterexample. The failure is still a sufficient-budget failure, not a counterexample to the true deletion bound.

## 4. What the example does and does not rule out

The exact argument rules out replacing the fixed theta3/4 bbox container by the proposed union of source-centered squares with the prescribed common half-side sqrt(theta L)/2. It does not rule out every possible radius function depending on theta and L, containers or covers using the actual captured masses, or all other theta. No theta optimization is attempted: the exact filled-container budget has the monotonicity already proved in R66.

For the full erosion union W defined in the companion R67 note, each label uses m(T) and m(T\S). In this particular input every crossing label meets q and some grid point, so its horizontal span is at least1/8, while its residual mass is only epsilon with sqrt(epsilon)=1/32. Its same-bbox residual core is empty. Internal labels also have empty residual core. Since the q singleton halo lies within the heavy singleton halo, all of H(A) is covered by labels meeting S. Consequently, for this one input,

    W_A(S)=H(A),  |W_A(S)|=Delta_A(S)+epsilon.           (6)

This demonstrates both why W removes artificial filled-container holes and why it can still overpay actual residual halo. It is not a proof of the proposed general W selection bound. The proposed general target and its distinct obligations are stated in `r067-support-rank-and-full-erosion-union.md`.

## 5. One exact control, with separate scopes

The author's `computation/check_r067_positive_padding.py` actually checked all1023 nonempty subsets, all46 complete bbox closures,18 eligible closures, uniqueness of maximum mask511, positivity, connectivity, the strict fat and nonterminal conditions, support-grid interval overlap and the rational sign(5). Python3.15.0, runtime0.016917seconds, exit0. The author minimal input `computation/inputs/r067-positive-grid-padding.json` has SHA `68cfb1045685d6c51de1a3d8d6f320154454c34e8021012b56ce81f15700202b`; the output is `computation/results/round067/positive-padding-arithmetic.json`. This source/arithmetic check does not calculate a halo union.

The separate frozen production geometry also actually ran this one input: the inward and outward support-neighborhood unions agree with the corresponding filled-container rectangles by zero directed symmetric differences. Its complete halo bounds give B_lower approximately3.50308723161359>3.5 and Delta_upper approximately0.999688823837595<1. One child process exited0 and was verified gone. These are decimal summaries of stored directed rational intervals, not substitutes for the analytic sign certificate(5). The producer adapter and the author minimal input retain separate source hashes.
