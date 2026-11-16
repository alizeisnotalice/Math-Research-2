# R82: robust shared-side payment for arbitrary measures on eight short axis segments

2026-11-16. Centered axis-parallel squares, one common radius, strict halo threshold1. The theorem varies the balanced old square width and permits arbitrary finite positive measures on eight short axis segments, together with unrestricted additional outside mass away from the effective corner windows. It is a direct finite-Borel-measure theorem under explicit square caps, not an atomic approximation or a general C2 result.

## 1. The actual measures and theorem

Let3/4<w<1, e=1-w, and

    sigma=(delta_(0,0)+delta_(w,0)+delta_(w,w)+delta_(0,w))/4.

Put T(x,y)=(w-y,x), P_0=(-1,w-1/2)^2 and P_i=T^i(P_0). In the reference frame define compact segments

    N=[-179/500,-44/125]x{0},
    F={0}x[-43/50,-17/20].

For i=0,1,2,3 let alpha_i,beta_i be arbitrary finite nonnegative Borel measures supported on N,F respectively, with total masses a_i,b_i. Let rho be any finite nonnegative Borel measure carried by the complement of union_i P_i. Define the one common outside and full measures

    nu=sum_i (T^i)_#(alpha_i+beta_i)+rho,
    mu=sigma+nu.

No atom-count bound is imposed on any bundle; a bundle may be zero. Assume explicitly that for every closed axis-parallel square Q,

    mu(Q)<=max(1,side(Q)^2).                         (1)

Let Delta=H(mu)\H(nu), where H uses open squares and the strict mass/area threshold1. Let Gamma_i be the four rotated corner squares, starting with Gamma_0=[-1/2,-e]^2. Then

    sum_i |Delta intersect Gamma_i|<62/125=.496,      (2)
    |Delta|<437/125-3(1-w)^2<7/2.                    (3)

For finite atomic inputs this is an unrestricted-rank outside theorem in the stated support geometry: the old S has mass1 and span squared w²<1, and(1) makes S a maximum eligible block. Uniqueness, connectedness and nonterminal status are not needed. In the general Borel case the hypothesis is(1) itself; no undefined extension of the project's finite eligible-block terminology is used.

The statement does not assert that the few necessary mass inequalities below imply all global caps. It upper-bounds every actual input satisfying(1) by a larger mass domain. There is no independently optimized corner pasting and no finite-to-Borel deletion-area limit.

## 2. Shared capacities, localization and measurability

All strict full and residual witnesses have r<1/2 by(1), since nu<=mu. For the southwest corner, every such open witness square lies strictly inside P_0: x,y in[-1/2,w-1] and r<1/2 imply all its coordinate endpoints are strictly between-1 and w-1/2. The rotated windows are separated by coordinate gaps1-w>0. Each contains exactly its own old corner and the corresponding two local segment measures. Other rotated old corners or segment measures have a coordinate at least w in this frame, beyond w-1/2. The reference negative segments lie strictly between-1 and0. The measure rho, including any mass on the boundaries of the windows, cannot enter these open witnesses.

Consequently full and residual halos on Gamma_i equal those of the local measure delta_0/4+alpha_i+beta_i and its outside alpha_i+beta_i, respectively. No rank-two label enumeration is used for this assertion.

The fixed side unit window[-43/50,7/50]x[0,1] contains both old left endpoints and the entire near bundle alpha_0 and rotated far bundle beta_3. Their old mass is1/2. Its rotations imply

    a_i+b_(i-1)<=1/2.                               (4)

Any extra captured rho or other outside mass only strengthens this necessary inequality. Thus

    0<=a_i,b_i<=1/2,   sum_i B_i<=2,
    B_i=a_i+b_i.

Likewise the fixed corner unit window

    [-179/500,321/500]x[-43/50,7/50]

contains the reference old corner and both entire local bundles. Rotations give

    B_i<=3/4.                                       (5)

These inequalities concern the same original nu; the additional rho is not charged independently to different corners.

All sets whose areas are used below are measurable. For fixed r, x↦mu(Q_r(x)) is lower semicontinuous by Fatou's lemma for the indicator of an open square. Hence each strict superlevel set is open, and H(mu), the union over r>0 of such open sets, is open; the same holds for nu. Moreover, a deleted point must have a strict full witness capturing positive old mass: otherwise that witness already belongs to H(nu). It is therefore within distance<1/2 of an old corner. Thus Delta lies in C=(-1/2,w+1/2)^2, directly under(1), also for Borel nu.

For a corner upper bound enlarge the local receiver crop from C_e=[e,1/2]^2 to C_0=[0,1/2]^2 in positive coordinates X=-x,Y=-y. This is an inclusion for the same local deletion set. We do not identify the enlarged local halo with the global halo away from the original Gamma_i.

## 3. Partial near capture and an actual residual inner rectangle

Fix one local near bundle alpha, mass a<=1/2. Write

    d_-=44/125, d_+=179/500, d0=71/200, delta=3/1000,
    v=1/4, u=sqrt(1/4+a)/2, z=sqrt(a)/2.

Consider a strict local witness containing the old atom, positive alpha mass and no far mass. It has r<u. Positive alpha capture is impossible if X<=d_--r, and old capture requires X,Y<r. Thus its receiver is inside

    N_out=(d_--u,u)x(-u,u).

Together with old-only witnesses it is covered by O union N_out, where O=(-v,v)^2. This upper bound uses the total a even if the witness captures only part of alpha.

The *actual* near residual contains

    R_in=(d_+-z,d_-+z)x(-z,z),                        (6)

interpreted as empty when its horizontal width is nonpositive. For any point strictly in this rectangle, every position of the entire near bundle has distance<z from the receiver. A common r between that maximum distance and z captures all of alpha and satisfies a=4z²>4r². This also covers endpoint atoms and continuous Borel mass, because the bounds refer to the entire compact segment. For a=0 the residual inner rectangle is empty.

Use a single virtual midpoint atom only as a comparison device, with its receivers

    N0=(d0-u,u)x(-u,u),
    R0=(d0-z,d0+z)x(-z,z),
    G_e(a,d0)=|[(O union N0)\R0] intersect C_e|.

No capacity or maximum-block property is asserted for that virtual measure. In particular we do not replace the real bundle by it and then reuse the original caps. The companion [uniform near-halo lemma](r082-dual-near-halo-uniform-payment.md) proves directly for all0<=a<=1/2 that

    G_(1/1000)(a,d0)<699/10000=.0699.                (7)

Its higher-mass branch includes the case where N0 covers the old corner. When a=0, the zero-weight mixed receiver is redundant and contained in O, so the convention agrees with the actual singleton measure.

The rectangle inclusions N0 subset N_out and R_in subset R0 permit a one-sided area comparison. Inside C_0, the added full strip has area at most delta*u, and the missing residual core has area at most2delta*z. The latter bound also holds when z<=delta and R_in is empty: the residual comparison rectangle after cropping to C_0 has area at most2z²<=2delta*z. Clipping at0 or1/2 can only reduce these differences. Therefore

    |[(O union N_out)\R_in] intersect C_0|
      <=G_0(a,d0)+delta(u+2z).                       (8)

Uniformly for a<=1/2,

    u+2z<=sqrt(3)/4+1/sqrt(2)<6/5.                   (9)

For example sqrt(3)/4<7/16 and1/sqrt(2)<3/4 give a strict sum<19/16<6/5. This comparison requires no differentiability at a contact event. For a single actual near atom it also gives the usual distance-Lipschitz bound; for arbitrary bundles the full outer and residual inner are separately justified by the same original alpha.

Finally

    |C_0\C_(1/1000)|=1/1000-1/10^6<1/1000.

Equations(7)-(9) give the robust near-plus-old cost

    |[(O union N_out)\R_in] intersect C_0|
       <.0699+.001+(6/5)(.003)=.0745<.075.           (10)

The actual complete residual may cover more than R_in, which only improves this bound. No finite-atom approximation is used.

## 4. Positive-capture classes and safe far/joint outer rectangles

Fix the two local bundles alpha,beta, with masses a,b and B=a+b. Every deleted point has an old-bearing strict witness. Classify such a witness by whether its captured near mass and its captured far mass are zero or positive. There are four possibilities, regardless of the number of atoms or whether the measures are atomic. Old-only and near-only are handled by Section3. A witness capturing no old mass cannot contribute to deletion because it is already a residual witness.

A far-only witness has r<sqrt(1/4+b)/2 and Y>D_--r, where D_-=17/20, because positive far capture is impossible below that radius-distance condition. A witness capturing positive mass from both bundles has r<sqrt(1/4+B)/2 and the same necessary lower Y constraint. In both cases old capture gives X,Y<r. Let F_out and J_out denote the resulting far and joint outer rectangles. At the level of sets, after restricting to the positive crop,

    Delta_local subset
        [(O union N_out)\R_in] union F_out union J_out. (11)

Indeed every deleted point avoids the full actual residual, and hence avoids the guaranteed subset R_in. The three sets on the right need not be disjoint: the sum of their areas is only an upper bound. There is no use of R81's alternating identities a=m or b=m. The upper bounds on F_out and J_out may ignore residual credit, but the actual deletion on the left always subtracts the same complete alpha+beta halo.

If a or b is zero the corresponding positive-capture class is absent; retaining its stated outer rectangle is harmless. The old-only case at a=b=0 is included in O.

For the far outer bound put r=sqrt(1/4+b)/2<=sqrt(3)/4. Its area in C_0 is at most

    r(2r-17/20)_+<=sqrt(3)/4 (sqrt(3)/2-17/20)<7/1000. (12)

The function is increasing in r in its positive range. The last strict inequality is exact: its left side is3/8-17sqrt(3)/80, and sqrt(3)>736/425 follows from3*425²-736²=179>0.

For the joint outer bound put U=sqrt(1/4+B)/2<=1/2 using(5). Its crop is contained in

    [0,U]x[17/20-U,U],

with area U(2U-17/20)_+. If this rectangle is nonempty, U>17/40>d_+, so discarding the near-distance lower X restriction causes no additional issue. A witness need only capture part of either bundle; using total B gives a valid upper radius for all such possibilities.

Let

    B0=(17/20)^2-1/4=189/400,
    kappa=1/2-(17/20)/4=23/80.

The upper function is zero at and below B0; for B>B0 it has derivative

    1/2-(17/20)/(4sqrt(1/4+B))<=kappa,

since B<=3/4. This gives the hinge bound

    |J_out intersect C_0|<=kappa(B-B0)_+.             (13)

At equality B=B0 the outer rectangle has zero height and contributes zero area. When U=1/2 it is only a necessary upper radius: every actual witness still has r<1/2. We do not claim that an arbitrary point in J_out admits an actual witness, because this outer rectangle has discarded source-position constraints. Including its boundary lines affects no Lebesgue area.

## 5. The arbitrary-mass shared budget

Let k be the number of positive hinge terms. Equations(4)-(5) imply

    sum_i(B_i-B0)_+<=min(2,3k/4)-kB0.

For k=0,1,2,3,4 the respective upper values are

    0, 111/400, 111/200, 233/400, 11/100.

Therefore for every allowed common mass vector,

    sum_i(B_i-B0)_+<=233/400,
    sum_i |J_out,i intersect C_0|<=5359/32000<21/125=.168. (14)

This uses shared capacities and individual corner caps simultaneously. It is a safe relaxation of their true feasible mass set, not an assertion that each extremizing vector can be jointly realized. No alternating assignment or equality of corner masses is required. The maximum of this hinge relaxation occurs with three positive terms; summing four separate worst-case joint costs would be too large.

Combining the four near estimates(10), four far estimates(12), and(14) in(11) proves

    sum_i |Delta intersect Gamma_i|
       <4(.075)+4(.007)+.168=.496,

which is(2). The corner boxes have side1/2-e=w-1/2. In the full container C=(-1/2,w+1/2)^2 the noncorner complement has area3-3e². All actual deletion lies in C, so paying this complement by its whole area yields(3). No noncorner residual credit is needed for this restricted theorem.

## 6. Qualification, verification and general gap

This theorem concerns every finite positive Borel input in the stated support class satisfying all square caps. Within the four effective windows the near and far masses may have arbitrary distributions on their short axis segments, with no rank bound; outside those windows rho is unrestricted apart from finiteness, positivity and the same global caps. It does not manufacture globally feasible inputs from independently feasible corners or replace(1) by(4)-(5).

The width interval3/4<w<1 and all bundle positions are covered analytically. Finite atomic outside measures of arbitrary rank are direct special cases, not limits. The Borel extension was established by positive-capture classes and an actual whole-bundle residual core, rather than a finite subset menu or a deletion-area convergence claim.

The sole imported local analytic estimate is(7), with its own proof and scalar receipt. The producer's separate small exact geometry/constants certificate checks necessary window bounds, crop/Lipschitz constants, the far inequality and hinge arithmetic. It does not run a parameter area search. The measure-theoretic capture argument, one-sided outer/inner comparison and shared budget above carry the uniform theorem; software arithmetic alone does not establish its Borel scope.

The old block remains the balanced four-corner measure. Arbitrary old profiles, general two-dimensional external supports *inside* the corner windows and the associated maximum-block selection problem are not covered. No general C2 bound changes. The next missing step is a comparable local cost bound with shared capacities for general old profiles and positions, rather than an unbounded sequence of fixed-support variants.
