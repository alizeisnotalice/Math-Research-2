# R72: a sharp-density Borel family and an atomic-transfer obstruction

2026-11-16, Asia/Shanghai. This proof concerns actual shared outside
measures. It distinguishes finite atomic sources from a diffuse Borel
outside component throughout. It does not prove the uniform-grid 3.5
budget for arbitrary finite atomic outside inputs.

## 1. Grid, caps and the actual deletion question

Let n>=2 be an integer and set

    e=1/n^2,  h=(n+1)/n^2,  w=1-e,
    S_n={0,h,...,(n-1)h}^2,
    sigma_n=(1/n^2) sum_(p in S_n) delta_p.

Thus sigma_n has mass 1, support width w<1, and is an eligible block.
For finite atomic outside inputs, remaining a maximum eligible block
of mass 1 is equivalent to the R64 closed-square caps

    (sigma_n+nu)(Q)<=max(1,R^2)                          (1)

for every closed axis-parallel square Q of side R. The same formula is
an explicit admissibility condition for general finite Borel nu; we do
not identify a diffuse measure with a finite atomic block catalogue.
All halos below use strict threshold 1 and centered open common-radius
squares. Actual deletion is H(sigma_n+nu) minus H(nu).

Every strict full witness under (1) has radius r<1/2. If it survives
outside the residual halo it must capture at least one internal source.
Consequently actual deletion is contained in

    C_n=[-1/2,w+1/2]^2,  |C_n|=(2-e)^2.                 (2)

This alone pays n=2: |Delta|<=49/16<7/2. The single-source n=1 case has
|Delta|<=1, since a mixed eligible set containing the mass-1 source
would have mass greater than 1. For n=3 the container exceeds 7/2 by
11/162; for large n its excess tends to 1/2. Closing that gap for
arbitrary finite atomic outside inputs remains the main open question.

## 2. The complete internal halo is small at large n

Completing an internal source subset to its bounding-box closure leaves
its box unchanged and increases its mass, so these complete labels
suffice to cover the internal halo. Every bbox closure of S_n is a
k-by-l consecutive grid rectangle.
Its mass is kl/n^2 and its maximum span is (max(k,l)-1)h. Eligibility
therefore means

    (max(k,l)-1)^2*(n+1)^2 < kl*n^2.                    (3)

If, say, k>=l+2, then kl<=(k-1)^2-1, and (3) is impossible because
h>1/n. The eligible complete internal labels are exactly:

* every k-by-k square, 1<=k<=n;
* k-by-(k+1) rectangles and their transposes when k(2n+1)<n^2.

The first claim follows from (k-1)(n+1)<kn for k<=n. Expanding (3) for
the second gives precisely the displayed strict arithmetic condition.

The receiver rectangle of a k-by-l label can extend to the left of the
whole source bbox by at most sqrt(kl)/(2n)-(k-1)h. If k=1, eligibility
forces l<=2, so this is at most sqrt(2)/(2n). If k>=2, l<=k+1 and
sqrt(k(k+1))<=k+1 give

    2(k-1)(1+1/n)-sqrt(kl)>=k-3>=-1>-sqrt(2).

Reflection and interchange of coordinates yield

    H(sigma_n) subset [-sqrt(2)/(2n),w+sqrt(2)/(2n)]^2.   (4)

On the other hand, the n^2 singleton receiver squares each have side
1/n. Their centers have spacing h>1/n, so they are pairwise disjoint and
have total area 1. Therefore

    1<=|H(sigma_n)|<=(w+sqrt(2)/n)^2 ->1.                (5)

The weak limit of sigma_n is Lebesgue measure on [0,1]^2, whose strict
threshold-1 halo is empty. Thus even this simple internal sequence shows
that halo area is not continuous under weak convergence at the strict
threshold. Formula (5) concerns the actual internal halo, not R70's
source-only virtual outer regions whose areas tend to 4.

## 3. The shared diffuse outside family satisfies every square cap

Set

    J_n=[-e,1]^2,  t=n/(n+1),  B=[-1,2]^2,
    nu_n=t*1_(B\J_n) dx,  mu_n=sigma_n+nu_n.            (6)

The density of nu_n is at most t<1, so H(nu_n) is empty. Its entire
contribution to H(mu_n) is therefore actual Borel deletion, not a
receiver-dependent choice of virtual outside mass.

We first prove a one-dimensional overlap estimate. If an interval I of
length R>=1 captures k>=1 grid points, then

    |I intersect [-e,1]| >= a_k,
    a_k=(k-1)h+e=[k(n+1)-n]/n^2.                        (7)

If I is contained in [-e,1], its overlap has length R>=1>=a_k. If I
contains [-e,1], its overlap has length 1+e>=a_k. Otherwise it crosses
exactly one endpoint of [-e,1]. Its captured grid points span at least
(k-1)h, and between the extreme captured point and that crossed endpoint
there is at least e additional overlap. This proves (7), including a
source exactly on an interval endpoint. For k=n it gives a_n=1.

Now take a closed square Q of side R>=1. If one coordinate interval
captures no grid point then sigma_n(Q)=0 and mu_n(Q)<=tR^2<=R^2. Otherwise
let its coordinate counts be k,l>=1. Using (7),

    mu_n(Q)<=tR^2 + kl/n^2 - t*a_k*a_l.

For f(k,l)=kl/n^2-t*a_k*a_l, its forward difference in k is

    f(k+1,l)-f(k,l)=(1-l/n)/n^2>=0,

and similarly in l. Hence f(k,l)<=f(n,n)=1-t, and

    mu_n(Q)<=tR^2+1-t<=R^2,   R>=1.                    (8)

For R<=1, include Q in a closed unit square and use (8) to get
mu_n(Q)<=1. It would be incorrect to assert mu_n(Q)<=R^2 in this range.
Thus (1) holds for every square, with no grid scan or finite outside
position assumption. Truncating the diffuse density to B only reduces
mass, so the comparison with the full complement of J_n is valid.

## 4. The density t is sharp within this constant-density family

Replace t in (6) by a constant c>=0, keeping the same B and J_n. Let

    a=(n+2)/n^2=e+h,
    Q=[-a,1-a] times [0,1].

For n>=2 this is a unit square contained in B. It captures exactly n-1
internal columns and all n rows, so sigma_n(Q)=(n-1)/n. Its part outside
J_n has area a-e=h. The unit cap consequently forces

    (n-1)/n+c*h<=1,  hence c<=n/(n+1)=t.                (9)

Together with Section 3, this proves that t is the exact largest allowed
constant density in this specified B\J_n family. It is not an extremal
classification among arbitrary outside measures.

## 5. The Borel halo has an exact finite distance-event criterion

Under (1), every strict witness has radius below 1/2. Since nu_n alone
has no halo, a strict witness captures an internal grid point p. Every
point in that witness square is then within distance 2r<1 of p, so the
square lies in [-1,w+1]^2, which is contained in B. Therefore truncating
B instead to [-1,w+1]^2 does not change H(mu_n); both truncations are
identical on all relevant witnesses, and both residual halos are empty.
The area component of this round uses the explicit B=[-1,2]^2 in (6).

For fixed x, let d range over the distinct distances ||p-x||_infinity
of the internal grid, including all ties and d=0 when present. Let
K_x(d) be the number of internal points at distance at most d, and set

    V_x(d)=|Q_closed(x,d) intersect J_n|,
    E_x(d)=K_x(d)/n^2-4(1-t)d^2-t*V_x(d).                (10)

Then exactly

    x in H(mu_n) iff E_x(d)>0 for at least one d<1/2.   (11)

For a strict open witness r, take the largest distance d<r of its
captured internal sources. Between source-capture events, the internal
mass is constant while 4(1-t)r^2+t*V_x(r) is strictly increasing for
r>0, since V_x is nondecreasing. Thus E_x(d) exceeds the witness excess.
Conversely a positive value in (10) remains positive after a sufficiently
small increase to an open radius r>d: the diffuse overlap is continuous,
and all sources counted at d are then strictly captured. This proves
(11). A zero value is not a strict witness. Overlap breakpoints of J_n
do not create additional maxima between source events.

Formula (11) is an exact actual-measure area interface. Directed
classification of every receiver cell can use all guaranteed and
possible source events and rigorous overlap extrema. Finite midpoint
sampling alone is not an area certificate. The separate adaptive-area
proof and saved complete cell certificates record that verification.

## 6. Canonical fine-grid atomization has maximum mass tending to 1

Fix n for this section. For delta=1/m tending to zero, take ALL closed
grid cells of side delta whose closure is contained in B\J_n. At the
center of each such cell put mass t*delta^2, forming the finite atomic
measure nu_(n,delta). Omitting the boundary-touching cells changes only
a vanishing-area boundary layer. Thus these outside measures converge
weakly to nu_n. Put mu_(n,delta)=sigma_n+nu_(n,delta) and let L_delta be
its actual maximum eligible block mass. We do not assume that S_n is
exactly maximizing for every delta.

For any square Q_R of side R, all cells whose centers belong to Q_R are
contained in the concentric square of side R+2delta. Their diffuse mass
is exactly their atomic mass. Therefore

    mu_(n,delta)(Q_R)<=mu_n(Q_(R+2delta)).                (12)

The use of 2delta is conservative; it includes every closed cell and
all boundary ties. For any eligible finite source subset of span R,
its mass m is at most the left side of (12) for a containing square.
All source positions lie in B, so R<=3. If R>=1, (8), (12), and R^2<m
imply

    m<=t(R+2delta)^2+1-t
      <t*m+1-t+4t*delta*(3+delta),

and hence m<1+4n*delta*(3+delta). If R<1, the full cap gives
m<=max(1,(R+2delta)^2)<=(1+2delta)^2, which is also at most that bound.
Since S_n remains an eligible block of mass 1,

    1<=L_delta<=1+4n*delta*(3+delta) ->1.                (13)

This is a genuine finite-atomic assertion. It does not claim exact
preservation of the maximum-block caps or exact uniqueness of S_n.

## 7. The residual halo quantitatively obstructs deletion-area transfer

Each outside atom in Section 6 has a singleton receiver square of side
sqrt(t)*delta<delta, strictly inside its own grid cell. These singleton
halos are disjoint. Count only cells wholly contained in C_n\J_n.
Because this region has a finite rectangular boundary, their total cell
area tends to |C_n\J_n|. Consequently

    liminf_(delta->0) |H(nu_(n,delta)) intersect C_n|
        >=t*|C_n\J_n|.                                (14)

This lower bound is from actual residual singleton halos; it does not
use a guessed limiting residual maximal function.

Any point in actual deletion for mu_(n,delta) has a witness containing
an internal source, with radius below sqrt(L_delta)/2. Hence deletion
lies in

    C_delta=[-sqrt(L_delta)/2,w+sqrt(L_delta)/2]^2.

By (13), C_delta decreases in size asymptotically to C_n (monotonicity
in delta is not asserted), and C_n is contained in C_delta because
L_delta>=1. Subtracting (14) from this valid full container gives

    limsup_(delta->0) Delta(mu_(n,delta),S_n)
      <= |C_n|-t*|C_n\J_n|
       = (1+e)^2 + (3-6e)/(n+1) =: U_n.                (15)

Here n>=2 ensures J_n is contained in C_n, and
|C_n\J_n|=(2-e)^2-(1+e)^2=3-6e. In particular U_5=1927/1250=1.5416,
and U_n tends to 1 as n tends to infinity, with the atomization limit
always taken first at fixed n. The same conclusion holds if every
outside mass is additionally multiplied by factors 0<lambda_delta<=1
with lambda_delta->1: (13) remains an upper bound, and singleton halo
coverage in (14) has limiting density t. This includes any such
mass-contraction repair that happens to restore exact caps.

Thus canonical atomic approximation cannot automatically transfer a
large Borel deletion area from (6). Even though the full maximum mass
tends to 1, the residual halo develops a large area that was identically
zero for the diffuse outside measure. A certified Borel deletion lower
bound strictly above U_n, as supplied for n=5 by the separate area
certificate, yields a strict discontinuity for this very atomization.
It is not a counterexample to the finite-atomic uniform-grid 3.5 bound.

## 8. The remaining actual finite-source obligation

The explicit Borel family gives one admissible common outside measure
and a sharp density parameter within that family. The adaptive area
bounds concern that specified measure, not the supremum over arbitrary
outside measures. The atomization theorem explains why simply replacing
its density by many small atoms does not preserve deletion area.

For arbitrary finite atomic nu that keeps S_n maximizing, the proposed
bound Delta<=7/2 remains open for n>=3. A strict violation for any chosen
maximum S_n can be uniquely maximized by the separate outside-contraction
reduction; if initially connected and strictly nonterminal, those open
structural conditions persist. Neither this reduction nor the Borel
family supplies such a finite-atomic violation.
