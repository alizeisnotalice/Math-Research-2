# R75 addendum: the exact unequal-mass orthogonal capacity polytope

2026-11-16. This addendum concerns capacities only. Its area conclusions do not extend beyond the separately certified equal-mass slice. All internal grid positions are retained.

Use the definitions n>=3, e=1/n², h=(n+1)/n², w=1−e and sigma_n from `r075-orthogonal-two-outside-full-grid-caps.md`. Fix0<a<=1 and put

    nu=b1 delta_(-a,0)+b2 delta_(0,-a),   b1,b2>=0,
    C_k(a)=max(1,[a+(k−1)h]²),
    U_n(a)=min_(1<=k<=n){C_k(a)−k²/n²},
    V_n(a)=min_(1<=k<=n){C_k(a)−k/n}.                (1)

**Theorem.** Every closed square satisfies

    (sigma_n+nu)(Q_R)<=max(1,R²)                    (2)

if and only if

    b1+b2<=U_n(a),  b1<=V_n(a),  b2<=V_n(a).       (3)

Thus(3), together with nonnegativity, is the exact feasible mass polytope for these two fixed outside positions and the complete internal grid. No separate cornerwise assignment or symmetrization is used.

Necessity of the sum inequalities follows from the squares[-a,(k−1)h]², as in the equal-mass proof. For the b1 inequalities, consider the source subset consisting of(-a,0) and the internal k-by-n prefix. Its mass is b1+k/n, and its bounding-box span is max(a+(k−1)h,w). A closed square of that side contains this entire subset, possibly more mass. Since w<1, its capacity is C_k(a); hence b1+k/n<=C_k(a). The b2 case is the coordinate transpose. This argument applies also when one outside mass is zero: including its zero-mass position in the indicated test subset does not invalidate a necessary upper bound.

Conversely, a square containing no outside point has mass at most1. A square containing outside points but no internal point has mass at most b1+b2<=U_n(a)<=1−e<1, using the k=1 row and a<=1. If both outside points and an anchored k-by-l internal prefix are captured, set j=max(k,l). Its mass is at most

    b1+b2+j²/n²<=U_n(a)+j²/n²<=C_j(a),

and its side is at least a+(j−1)h. This gives(2). If exactly the left outside point and an anchored k-by-l prefix are captured, then

    b1+kl/n²<=b1+k/n<=V_n(a)+k/n<=C_k(a).

Again its side is at least a+(k−1)h. The bottom-only case is the transpose. These cases exhaust all square captures, including closed boundary ties, and prove sufficiency.

## Two decisive rows for both constraints

For0<a<=e, the k=n rows give U_n(a)=V_n(a)=0. For e<a<=1 let

    K=1+floor((1−a)/h),    1<=K<=n−1.

Then

    U_n(a)=min{1−K²/n², (a+Kh)²−(K+1)²/n²},
    V_n(a)=min{1−K/n,   (a+Kh)²−(K+1)/n}.        (4)

The U formula was proved in the equal-mass interface. For V, the unit-side rows1−k/n decrease with k. On the above-unit branch R_k=a+(k−1)h>1, successive differences are

    [R_(k+1)²−(k+1)/n]−[R_k²−k/n]
       =2hR_k+h²−1/n
       >2h+h²−1/n=(n+2)/n²+h²>0.

Thus the only candidates are again the largest unit-side prefix and the first nonunit prefix. Exact unit-side ties are included in K by the floor convention.

## Relation to the equal-mass result and scope

For every k,

    C_k−k²/n² <= 2(C_k−k/n),

because the difference of the right side and left side is
C_k−1+(1−k/n)²>=0. Taking k that minimizes V yields U_n<=2V_n. Therefore on the equal slice b1=b2=b, the individual inequalities in(3) are redundant and the exact condition is b<=U_n/2=B_n, consistently with the frozen interface.

Formula(3) retains the coupling b1+b2<=U_n alongside the individual capacities. It is an explicit all-n reduction of the full source constraints, but it is not an unequal-mass deletion-area theorem. Any parameter-box proof for b1=b2 must remain restricted to that slice. No finite scan, weak-limit transfer, or collapse of internal positions is used to prove this addendum.
