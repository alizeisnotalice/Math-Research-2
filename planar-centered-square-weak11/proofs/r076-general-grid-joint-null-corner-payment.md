# R76: general-grid corner payment without a joint mixed witness

2026-11-16. All thresholds are strict, the maximal basis consists of centered open axis-parallel squares with one common radius, and the outside measure is finite positive atomic. This proof concerns an actual common outside measure and keeps the complete internal grid in every capacity condition. No assumption is made about a fictitious joint receiver after collapsing internal mass.

## 1. Objects and the actual joint-null condition

Let n>=3, e=1/n², h=(n+1)/n², w=1−e, K=floor(n/2). The internal measure sigma_n places mass e at each point of S_n={0,h,...,(n−1)h}². Let the same finite outside measure nu satisfy all closed-square caps

    (sigma_n+nu)(Q_R)<=max(1,R²).                    (1)

Equivalently, S_n remains a maximum eligible block of mass1. Write

    D=H(sigma_n+nu)\H(nu),
    C=[−1/2,w+1/2]²,
    Gamma=[−1/2,−e]²,
    P=(−1,1/2−e)²,
    m_P=sigma_n(P)=K²/n².

The other corners and source windows are the D4 reflections about the internal square. Boundaries of receiver regions do not affect area. The R74 localization theorem implies that every strict full or residual witness has radius r<1/2, each corner witness lies strictly within its P, and corner deletion is unchanged by restricting both measures to P. Restricting outside mass while retaining the COMPLETE sigma_n preserves every cap in(1).

Assume P contains at most two distinct positive outside positions. Coincident positions are combined before counting. For two positions write the restriction as b1 delta_q1+b2 delta_q2, with b1,b2>0. The **actual joint-null condition** means that no receiver x in Gamma has a strict full witness Q(x,r) which captures any internal source and both q1,q2. It is a condition on the original source positions and masses, not on an auxiliary collapsed input. With zero or one effective outside position this condition is vacuous.

We prove the strict auxiliary half-local-mass bound

    |D intersect Gamma| <= g_n < m_P/2,             (2)

where

    A0=(n−2)²/(4n⁴),
    g_n=A0+(1−2/n²)m_P/2.                           (3)

The positive gap is explicit:

    m_P/2−g_n = (n−1)/n⁴             if n is even,
                 (2n−3)/(4n⁴)        if n is odd.    (4)

In particular g_n<1/8. Formula(2) is only asserted under the actual joint-null condition; active joint witnesses remain outside this proof.

## 2. Pure internal corner witnesses capture only the origin

For x in Gamma, any square capturing an internal source captures an anchored k-by-l prefix, because all internal coordinates are nonnegative and both coordinates of x are negative. Let j=max(k,l). Its captured internal mass is kl/n²<=j²/n². A strict pure-internal witness must therefore satisfy r<j/(2n). If j>=2, capturing the furthest prefix coordinate requires

    r>(j−1)h+e.

But

    2n²[(j−1)h+e−j/(2n)] = n(j−2)+2j >0,

a contradiction. Thus a pure internal corner witness captures only the origin. Its radius threshold is1/(2n), and its halo crop in Gamma has area

    (1/(2n)−e)²=A0.                                (5)

This argument uses actual grid locations. Replacing the entire internal measure by its total mass at the origin would give an incorrect value for this term.

## 3. A uniform clipped annulus estimate

Let0<=v<=u<=1/2, 0<e<1/2, and q be arbitrary. Then

    |Gamma intersect Q_0(u) intersect [Q_q(u)\Q_q(v)]|
       <=(u−e)_+²−(v−e)_+²
       <=(1−2e)(u²−v²).                            (6)

Here Q_p(t) is the open square of radius t about p; closed representatives may be used to calculate area. If u<=e the left side is zero. Otherwise Gamma intersect Q_0(u) is I² with I=[−u,−e] and length L=u−e. For either coordinate write A_t=|I intersect[q_j−t,q_j+t]|. Then

    0<=A_v<=A_u<=L,       A_u−A_v<=u−v.             (7)

When q_j lies outside I, overlap is min(L,(t−d)_+) for a fixed distance d>=0 and is1-Lipschitz in t. When q_j lies inside I, the radius-u interval covers I because L<=u, while radius v captures length at least min(v,L). Hence A_u−A_v<=max(L−v,0)<=u−v. This proves(7) in every position.

The two-coordinate product difference is maximized under(7) by taking both outer lengths equal to L and both inner lengths at least(L−(u−v))_+. One can see this by noting that xy−(x-d)_+(y-d)_+ is nondecreasing in x,y on[0,L] for d>=0. Thus the area difference is at most

    L²−(L−(u−v))_+²=(u−e)²−(v−e)_+².

Finally F(s)=(sqrt(s)−e)_+² has derivative0 for s<=e² and derivative1−e/sqrt(s) for s>=e². On[0,1/4] it is continuously differentiable with0<=F'(s)<=1−2e. Integrating from v² to u² gives the second inequality in(6). This is an all-position annulus estimate, not a cap assertion for an auxiliary source measure.

## 4. Truncate single-outside collapsed radii before using the annulus

Consider a point of D intersect Gamma and one of its strict full witnesses. It must capture some internal mass: otherwise the same square would be a strict residual witness. Under actual joint-null, it captures either no outside atom or exactly one of b_i delta_qi. The former case is accounted for by(5).

For the latter case, any captured internal point is northeast of the origin, so the witness also contains the origin. Its internal mass is at most m_P. Define

    u_i=min(1/2,sqrt(m_P+b_i)/2),
    v_i=sqrt(b_i)/2.                                (8)

The real witness radius obeys r<1/2 by(1) and r<sqrt(m_P+b_i)/2 by its strict mass inequality. Hence its center lies in Q_0(u_i) intersect Q_qi(u_i). Also D avoids the true singleton residual halo Q_qi(v_i), because that halo is contained in H(nu).

The singleton cap gives b_i<=1. Therefore

    0<=v_i<=u_i<=1/2,
    u_i²−v_i²=min((1−b_i)/4,m_P/4)<=m_P/4.          (9)

The cap b_i<=1 is obtained by applying(1) to a sufficiently small closed square containing the atom. It does NOT imply m_P+b_i<=1: the original internal mass in P and q_i need not fit in a unit square. The truncation in(8) is essential to avoid that unsupported inference. No global cap is imposed on the collapsed auxiliary measure.

By(6) and(9), the possible single-outside contribution has area at most(1−2e)m_P/4 for each i. There are at most two such choices, and the complete residual halo can only remove additional points. Combining them with the actual pure-internal area gives

    |D intersect Gamma|<=A0+2(1−2e)m_P/4=g_n.

Only this union bound is used; it does not independently assign different outside masses or positions to different receivers. It never tests a joint label in the auxiliary collapsed geometry. Equations(4) follow by expanding

    g_n−m_P/2=[(n−2)²/4−K²]/n⁴,

and both gaps are strictly positive for n>=3. This proves(2).

## 5. An exact K²-prefix test for the actual condition

For fixed distinct q1,q2 in P and positive masses b1,b2, define the actual internal prefixes

    T_(k,l)={0,h,...,(k−1)h}x{0,h,...,(l−1)h},
    1<=k,l<=K.

For the bbox [A,B]x[C,D] of T_(k,l) union{q1,q2}, put

    delta_(k,l)=max{(B−A)/2,(D−C)/2,B+e,D+e,
                       −1/2−A,−1/2−C,0}.          (10)

This is exactly the minimum, over the closed receiver square Gamma, of the largest distance to these actual sources. In one coordinate the containing-center interval[B−r,A+r] meets[−1/2,−e] if and only if r is at least the corresponding three terms in(10); the two coordinates give the maximum.

Under(1), actual joint-null is equivalent to the complete finite system

    kl/n²+b1+b2 <=4 delta_(k,l)²
                  for every1<=k,l<=K.              (11)

Indeed any strict full witness capturing old mass and both outside atoms captures an anchored prefix T_(k,l), so its exact captured mass violates(11). Conversely, a strict reverse inequality in one row makes the corresponding label's strict receiver intersect Gamma. To justify the required radius below1/2, let L=kl/n²+b1+b2 and let d be the label's bbox span. Then delta>=d/2. If L>4delta² and L>1, the label would violate its enclosing-square cap because d²<L and1<L. Thus L<=1, and delta<sqrt(L)/2<=1/2. A radius strictly between delta and sqrt(L)/2 exists; its center can be moved slightly into the open Gamma while preserving the strict margin. Any extra source mass captured only strengthens the full witness. Equality in(11) creates no strict witness.

All K² prefixes must be retained. Arbitrary outside positions do not justify replacing this menu by diagonal k=l rows. Equivalently, set

    J(q1,q2)=min_(1<=k,l<=K){4delta_(k,l)²−kl/n²}.

For an actual positive pair satisfying the global caps, joint-null holds precisely when b1+b2<=J. This criterion uses the same actual outside mass sum in every row. Failing the condition says that a joint witness exists, not that the area bound is violated.

## 6. Fixed-position full-grid mass compatibility

A useful sufficient filter can combine(11) with the complete-grid mass polytope. For fixed distinct outside positions q1,q2, enumerate the nonempty bbox closures of S_n union{q1,q2}. Each row has capacity max(1,d_bbox²) minus its captured internal mass, and outside coefficients among(0,0),(1,0),(0,1),(1,1). The zero-outside rows are already valid for sigma_n.

Let V1,V2 be the minimum right sides among the(1,0),(0,1) rows and U the minimum among the(1,1) rows. These three classes are nonempty: singleton outside positions and their joint bbox supply rows. The exact admissible nonnegative mass polytope for the reduced complete-grid input is

    0<=b1<=V1, 0<=b2<=V2, b1+b2<=U.                (12)

This is the R64 full-closure criterion grouped by its three coefficient types. The largest feasible total mass is exactly U_eff=min(U,V1+V2). If other outside sources exist in the original input, deleting them yields(12) as a safe outer domain for the actual pair; it is not permission to add their masses again independently elsewhere.

If

    U_eff+kl/n²<=4delta_(k,l)²
                 for every1<=k,l<=K,               (13)

then every admissible positive pair at these fixed positions is in the proven joint-null branch. Coincident outside positions should first be merged and treated by the single-position theorem. Formula(13) is a sufficient all-mass certificate for these positions, not a necessary condition for the area payment itself. The inequalities for active joint witnesses remain an unresolved branch.

## 7. Four corners of one globally feasible input

Suppose the same input satisfying(1) has at most two positive outside positions in each effective P and the actual joint-null condition holds in every two-position corner. Sources elsewhere are arbitrary. Applying(2) separately to the four restrictions is legitimate because the common global caps are retained and each corner's actual full and residual memberships localize exactly. It does not construct a new measure by pasting independent corner choices.

All deletion lies in C. Its noncorner area is3−3/n⁴. Hence the complete true deletion obeys

    Delta<=3−3/n⁴+4g_n
       =7/2−4/n³+1/n⁴                         if n is even,
        7/2−1/n+1/(2n²)−2/n³                  if n is odd. (14)

Both expressions are strictly below7/2. At n=3, g_3=5/108 and(14) is85/27, agreeing with the independently proved three-source branch. This is a general-n arbitrary-position joint-null structural class, not the full two-outside or general maximal-operator theorem. Neither arbitrary active joint witnesses nor more than two effective outside positions per corner are settled here.

## 8. Actual scalar checks and their limits

`computation/round076/check_general_joint_null_constants.py` was actually run with Python3.15.0 and standard-library Fraction. For n=3,...,129 it checks127 parity/payment rows, all8255 inequalities for pure-prefix maximum index j=2,...,n, and635 clipped-radius squared identities at five fixed outside masses including0 and1. It recovers g_3=5/108 and the complete n=3 bound85/27. PID50576 exited0 and was confirmed absent; script time was about0.052 seconds.

This finite regression checks the written constants and truncation algebra. The all-n and all-mass conclusions are proved in Sections2–4, not inferred from these examples. The run does not sample outside positions, certify any actual finite input, evaluate the K² activity menu, or integrate a halo area. The n3 geometric proof, its own small kernel and the separate fixed-position probe belong to different components and are not counted as executions of this script.
