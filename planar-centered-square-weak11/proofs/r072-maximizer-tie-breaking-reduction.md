# R72: outside-mass contraction removes maximum-block ties

2026-11-16, Asia/Shanghai. General finite-atomic reduction. This does not
prove the proposed 3.5 deletion budget or determine C2.

## 1. The contraction and the unique maximum

Let A be a finite set of distinct points with strictly positive masses.
A nonempty subset T is eligible when d(T)^2<m(T), where d is its larger
coordinate span. Let L be the largest eligible mass and let S be ANY
eligible subset with m(S)=L. Write sigma for the masses on S and nu for
the masses on A\S. For 0<lambda<1 set

    mu_lambda=sigma+lambda*nu.

All positions remain in the positive support. The subset S keeps its
mass L and remains eligible. No formerly ineligible subset can become
eligible, because its mass only decreases. Any eligible T different
from S has contracted mass strictly below L: if T is a proper subset
of S, positivity already gives m(T)<L; otherwise T includes outside
mass and its mass strictly decreases from a value at most L whenever
it was originally eligible. Therefore S is the UNIQUE maximum eligible
block of mu_lambda, with maximum mass exactly L, for every 0<lambda<1.
There is no need to assume an original gap between tied maximum masses.

## 2. Actual deletion area converges

At strict threshold 1, write

    H(mu)={x: mu(Q_open(x,r))>4r^2 for some r>0},
    Delta(mu,S)=|H(mu)\H(nu)|.

As lambda increases to 1,

    H(mu_lambda) increases to H(sigma+nu),
    H(lambda*nu) increases to H(nu).                     (1)

The inclusions follow from monotonicity of mass. Conversely, a point
in a limiting strict halo has one fixed witness radius with positive
mass excess. That same strict witness works for all sufficiently large
lambda<1, proving the equality of the increasing unions. An equality
at threshold does not enter the strict halo and is not used here.

Since H(lambda*nu) is contained in H(mu_lambda), deletion indicators are
exactly the difference of these two halo indicators. By (1) they converge
pointwise to the original deletion indicator. All halos lie in a common
bounded set: any strict witness has radius at most sqrt(m(A))/2 and
contains at least one point of the fixed finite support. Thus dominated
convergence yields

    Delta(mu_lambda,S) -> Delta(sigma+nu,S).              (2)

The deletion SETS need not be monotone in lambda; subtracting two nested
halo families does not justify that assertion. The argument proves
convergence of their indicator differences and areas. Equivalently, the
absolute area difference is bounded by the sum of the two monotone halo
area losses. No differentiability, generic endpoint order or fixed halo
combinatorial chamber is required.

## 3. Structural hypotheses that survive

Connectivity of the eligible hypergraph survives for lambda sufficiently
close to 1. Choose finitely many originally strict eligible hyperedges
that connect the positive support. For each such T involving outside
mass, its eligibility remains whenever

    (1-lambda)*nu(T) < m(T)-d(T)^2.                       (3)

The finitely many positive right sides give one common lambda_0<1.
Pure internal hyperedges do not change. In fact, by applying (3) to all
original eligible subsets, the complete eligible hypergraph is unchanged
for lambda sufficiently close to 1; no new eligible subsets can appear.
This does not mean their masses or maximum ties remain unchanged.

The whole positive-support bounding box remains exactly the same for
lambda>0. Consequently the strict nonterminal condition used in R63,

    d(A)^2 > (35/33)*m(A),

is preserved by this contraction, since total mass decreases. Any fixed
condition depending only on S and its unchanged masses or positions
also survives: for example prescribed uniform S, its two-axis spans,
its source rank and its support-tree cost. No claim is made that every
possible mass-dependent ancillary equality or inequality survives.

## 4. Existential and universal maximum-block budgets

Let C be a class of finite positive atomic inputs with the following
property: for every A in C and every maximum eligible S, the contractions
above belong to C for lambda sufficiently close to 1. The unrestricted
finite-input class has this property. So does the R63 connected strict
nonterminal domain, with any additional conditions solely on the chosen
unchanged S. For a fixed K, compare the assertions

    (E) every A in C has some maximum S with Delta(A,S)<=K*L;
    (U) every A in C and every maximum S satisfy Delta(A,S)<=K*L.

When a class is described with conditions on a chosen S, the quantifiers
refer to the corresponding admissible pairs (A,S); the existential
version must range over the admissible maximum blocks. Then (E) and (U)
are equivalent whenever the contracted pair remains in the same class.
(U) implies (E) immediately. To prove the reverse, fix any admissible
(A,S), contract its outside mass, and apply (E) to the contracted input.
Its unique maximum is S, so Delta(mu_lambda,S)<=K*L. Taking lambda to 1
and using (2) proves the desired bound for the original S.

In particular, if ONE specified maximum block in the open connected
nonterminal domain has Delta(A,S)>K*L with a strict positive gap, then
for all sufficiently large lambda<1 the contracted input has S as its
unique maximum, remains in that domain, and still violates the bound.
Thus that single-block strict violation can be converted into a genuine
counterexample to the maximum-block selection assertion. One need not
first make every original tie expensive. All geometric and numerical
claims about such a proposed input must still be certified.

The non-strict sign is essential. If the only assumed statement is a
strict inequality Delta<K*L, the limiting argument gives <=K*L and need
not preserve strictness. This note asserts equivalence for the non-strict
budget, including K=7/2, not an unqualified strict version.

## 5. Relation to the uniform-grid question

For S=S_n fixed as the equal-mass uniform grid, the contraction changes
only the actual outside masses. It preserves S_n, maximum mass 1 and
all of the common-source geometry. Therefore a strict counterexample
for a specified S_n can be made uniquely maximizing. If its original
input was connected and strictly nonterminal, those hypotheses survive.
This reduces a quantifier issue but supplies neither a uniform-grid
area upper bound nor a counterexample. Previous individual tie counts
remain true historical facts; this reduction changes their role in
future counterexample verification, not those archived counts.
