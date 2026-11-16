# R75: exact full-grid capacities for an orthogonal two-outside family

2026-11-16, Asia/Shanghai. All halo thresholds are strict, the basis is centered open axis-parallel squares with a common radius, and outside sources are finite positive atoms. This proof retains every internal grid point. It does not use the false general two-outside northeast-collapse payment. The general uniform-grid four-corner and7/2 deletion problems remain open.

## 1. Family and exact all-square capacity theorem

Let n>=3, e=1/n², h=(n+1)/n², and w=(n−1)h=1−e. Give each point of S_n={0,h,...,(n−1)h}² mass e and write sigma_n for the resulting mass1 measure. Consider

    nu_(a,b)=b delta_(-a,0)+b delta_(0,-a),
    0<a<=1, b>=0.

At b=0 the outside measure is empty. For positive b these are two distinct outside positions. Define

    B_n(a)=1/2 min_(1<=k<=n)
             {max(1,[a+(k−1)h]²)−k²/n²}.             (1)

**Theorem.** The complete measure sigma_n+nu_(a,b) satisfies every closed-square capacity

    mu(Q_R)<=max(1,R²)                                (2)

if and only if b<=B_n(a). By the archived R64 finite-source theorem, this is exactly the condition that S_n remain a maximum eligible block of mass1. This is an all-position, all-side-length theorem for the specified orthogonal equal-mass family, not a restriction to unit windows.

Necessity follows by testing Q_k=[−a,(k−1)h]². It captures both outside atoms and exactly k² internal atoms, including all boundary atoms. Its mass2b+k²/n² and side a+(k−1)h give each inequality in(1).

For sufficiency, take any closed square Q_R and divide according to its outside capture. If it captures no outside atom, its total mass is at most1. If it captures outside atoms but no old point, its mass is at most2b. The k=1 row gives

    2b+e<=max(1,a²)=1,

so this branch also obeys(2).

If both outside atoms and some old points are captured, the captured old points form an anchored k-by-l prefix: the square's lower coordinate endpoints are at most−a and its upper endpoints include0. Set j=max(k,l). The side R is at least a+(j−1)h. Hence

    mu(Q_R)=2b+kl/n²
       <=2b+j²/n²
       <=max(1,[a+(j−1)h]²)
       <=max(1,R²).

If exactly(-a,0) is captured together with old points, the captured old points again form an anchored k-by-l prefix, since the square's lower x endpoint is at most−a and lower y endpoint is at most0. Put

    C_k=max(1,[a+(k−1)h]²).

The side R is at least a+(k−1)h, and the k-th row of(1) gives b<=(C_k−k²/n²)/2. Since l<=n,

    2kl−k²<=2kn−k²=n²−(n−k)²<=n²<=n² C_k.

Consequently

    b+kl/n² <= (C_k−k²/n²)/2+kl/n² <= C_k
             <=max(1,R²).

The other single-outside case is its coordinate transpose. This exhausts all closed squares, including ties, and proves sufficiency. In particular the apparently additional single-outside rectangle constraints are redundant for this family; no independent assignment of its two masses was made.

The function B_n is nonnegative and nondecreasing, being the minimum of nonnegative nondecreasing functions. The k=n row shows B_n(a)=0 for0<a<=e. For a>e every row is positive, so B_n(a)>0. If0<b<B_n(a), every outside-containing set has strict capacity slack when its span is below1 and is ineligible when its span is at least1. Thus S_n is the unique maximum block. This uniqueness statement does not assert connectivity or a nonterminal condition throughout the parameter family.

### Only the two prefixes adjacent to unit side can bind

For e<a<=1 define

    K=1+floor((1−a)/h),   1<=K<=n−1.

The complete minimum in(1) reduces to just two rows:

    2B_n(a)=min{1−K²/n²,
                    (a+Kh)²−(K+1)²/n²}.             (1a)

Indeed a+(k−1)h<=1 exactly for k<=K, and the corresponding quantities1−k²/n² decrease with k. For k>=K+1 put R_k=a+(k−1)h>1. Adjacent differences on this branch are

    [R_(k+1)²−(k+1)²/n²]−[R_k²−k²/n²]
       =2hR_k+h²−(2k+1)/n²
       >2h+h²−(2n−1)/n²
       =3/n²+h²>0

whenever an adjacent index k+1<=n exists. Thus the first above-unit row is the smallest on that branch. At a=e the n-th unit row gives0; at exact unit-side ties the floor convention includes the tied row in K. Formula(1a) therefore locates both possible binding constraints: the largest prefix fitting in the unit-side regime and the first prefix beyond it. It is not valid to discard the latter merely because unit-square capacities have been checked.

## 2. The exact n=3 boundary

For n=3, h=4/9, and put

    alpha=(3sqrt(14)−8)/9,
    beta=2sqrt(3)/3−4/9.

Then1/9<alpha<5/9<beta<1, and the minimum in(1) is

    B_3(a)=0                              for0<a<=1/9;
           ((a+8/9)²−1)/2                 for1/9<=a<=alpha;
           5/18                           foralpha<=a<=5/9;
           ((a+4/9)²−4/9)/2               for5/9<=a<=beta;
           4/9                            forbeta<=a<=1.  (3)

The formulas agree at the stated endpoints. To verify(3), the three candidate functions are4/9, [max(1,(a+4/9)²)−4/9]/2, and [max(1,(a+8/9)²)−1]/2. The switch from the third candidate to the second candidate's plateau solves(a+8/9)²=14/9, giving alpha. The later switch from the second candidate to4/9 solves(a+4/9)²=4/3, giving beta. Above5/9 the third minus the second candidate is positive, because their doubled difference is8a/9+1/27. Belowalpha the third is smaller than the plateau, and betweenalpha and5/9 it is larger. These comparisons prove the formula on the full domain.

For the scaled no-cap counterexample position a=2/3, this gives

    B_3(2/3)=32/81 < 4/9.

The decisive k=2 square is[-2/3,4/9]². If both outside masses were4/9, its mass would be4/3 while its squared side is100/81. Thus the admissible family must reduce each mass by at least4/81 at this position. Testing only unit windows misses this restriction.

## 3. Exact three-source corner geometry for admissible n=3 parameters

Assume0<=b<=B_3(a),1/9<=a<=1. The southwest receiver corner is

    Gamma=[−1/2,−1/9]²,
    P=(−1,7/18)².

Only the origin is an internal source in P, of mass m=1/9. The R74 all-radius localization theorem applies because the full S_3 capacities have been established, and both full and residual strict witnesses have radius below1/2. At a=1 the outside atoms are on P's excluded boundary; the same geometric formulas below give no corner contribution from them, so including this endpoint is harmless.

Write

    r0=1/6, u=sqrt(m+b)/2, v=sqrt(m+2b)/2,
    s=sqrt(b)/2, t=sqrt(2b)/2.

Up to boundary sets of zero area, all full-minus-residual membership on Gamma is exactly

    (A0 union A1 union A2 union A12)
       minus (R1 union R2 union R12),                 (4)

where every rectangle is clipped to Gamma, empty rectangles are discarded, and

    A0  =[-r0,r0]²,
    A1  =[-u,u−a] x [-u,u],
    A2  =[-u,u] x [-u,u−a],
    A12 =[-v,v−a]²,
    R1  =[-a−s,−a+s] x [-s,s],
    R2  =[-s,s] x [-a−s,−a+s],
    R12 =[-t,t−a]².

These are the complete subset receiver rectangles for the three local sources. A source label has a strict witness precisely in the open rectangle with coordinates max_source−sqrt(mass)/2 and min_source+sqrt(mass)/2. An ineligible label has an empty rectangle. Labels consisting only of outside atoms occur in both full and residual halos and cancel before(4); omitting A12 would be incorrect. Formula(4) keeps the true complete residual halo and the actual corner clipping. It does not identify the entire southwest quadrant with Gamma.

## 4. Uniform parameter-box upper certificate interface

The following interface permits a finite certificate for a continuum of admissible parameters. It does not by itself assert that a requested numerical budget has passed.

Take a parameter rectangle [aL,aU]x[bL,bU] with1/9<=aL<=aU<=1 and0<=bL<=bU. Since B_3 is nondecreasing, there are no admissible parameters in it if bL>B_3(aU). Otherwise bU may be reduced to min(bU,B_3(aU)). This tests the full capacity boundary(1), not only the unit rows.

Let uU,vU be rational upper bounds on sqrt(m+bU)/2 and sqrt(m+2bU)/2, and let sL,tL be rational lower bounds on sqrt(bL)/2 and sqrt(2bL)/2. Uniformly over the remaining parameter box, the four A rectangles are contained in

    A0,
    [-uU,uU−aL] x [-uU,uU],
    [-uU,uU] x [-uU,uU−aL],
    [-vU,vU−aL]².                                     (5)

The three R rectangles uniformly contain

    [-aL−sL,−aU+sL] x [-sL,sL],
    [-sL,sL] x [-aL−sL,−aU+sL],
    [-tL,tL−aU]².                                     (6)

Clip(5) and(6) to Gamma and discard empty rectangles. Their exact union-difference area is an upper bound for every actual deletion area(4) in the parameter box. Indeed upper full rectangles enlarge all possible full labels, whereas(6) is a subset of each parameter's residual union; subtracting this common residual subset can only enlarge deletion.

The bounds on sL,tL can be checked by squaring rational endpoints; uU,vU are checked similarly. The rectangle areas and union operations then use rational arithmetic alone. The omitted0<a<1/9 range admits only b=0. On Gamma, all internal sources except the origin are outside P; the pure internal singleton then gives area(1/6−1/9)²=1/324. A finite partition of[1/9,1]x[0,4/9] into such parameter boxes, each either pruned or certified below a target, would establish that target for the entire remaining closed parameter domain. Midpoint tests or a list of successful sampled parameters are not such a proof. The finite producer certificate, if completed, must record the covering partition, every retained upper bound and the exact cap reductions.

## 5. Scope

The new general result proved here is the exact all-n capacity characterization(1), including the redundancy of single-outside rows and the explicit nonunit n=3 obstruction. It reduces a genuine two-outside family while retaining full S_n geometry. The separate no-cap three-source example refutes the attempted general northeast-collapse m/2 theorem; it does not violate(1), because it is outside this feasible family. A full n=3 family area certificate requires the parameter-cover verification in Section4 and is recorded separately if obtained. No claim about arbitrary outside positions, unequal outside masses, or general four-corner7/2 payment is made by(1) alone.
