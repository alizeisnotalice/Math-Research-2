# R79: a continuous quadratic envelope for actual discrete prefixes

2026-11-16. The threshold is 1, squares are centered, axis parallel, and use
one common radius. This is a safe upper bound for the original grid count.
It does not replace the old atoms by a continuous measure.

## 1. Actual input and radius domain

For n>=2 write

    e=1/n², h=(n+1)/n², c=e/h²=n²/(n+1)², K=floor(n/2),
    sigma_n=e sum_(i,j=0)^(n-1) delta_(ih,jh),
    Gamma=(-1/2,-e)², P=(-1,1/2-e)², m_P=e K².

Let nu be one finite positive atomic outside measure and assume the actual
mu=sigma_n+nu satisfies all closed square caps

    mu(Q)<=max(1,side(Q)²).                              (1)

All statements about residuals use this same nu. For clarity suppose nu|P
has at most two positive locations q1,q2 with masses b1,b2; the argument
works for any fixed finite number, with all its outside subsets retained.
Outside P there may be further sources. Define H(eta) by the existence of
an open square Q(x,r) with eta(Q(x,r))>4r², and Delta=H(mu)\H(nu).

By (1), every full or residual strict witness has r<1/2. If x in Gamma and
r<1/2, Q(x,r) is contained in P. Thus sources outside P do not contribute to
either halo in this corner. A witness containing any old source necessarily
contains the origin, since the source coordinates are nonnegative and the
receiver coordinates are negative. In positive receiver coordinates
x=(-X,-Y), it therefore has

    e<X,Y<1/2,  max(X,Y)<r<1/2.                         (2)

We also allow X=e or Y=e in estimates on the closed crop; that does not add
boundary points to the definition of the open corner. Strict closed-square
witnesses and strict open-square witnesses give the same halo: slightly
increase a strict closed witness radius while keeping its mass surplus and
r<1/2. The argument below always counts the open square at its stated radius.

## 2. Counting the original old atoms, including integer events

Put a=(r-X)/h>0 and b=(r-Y)/h>0. The captured old coordinates on the first
axis satisfy ih<r-X, including i=0. Their number is exactly ceil(a), even
when a is an integer: at a=j the strict inequality captures i=0,...,j-1.
There is no saturation at n because r-X<1/2-e<w=1-e. The same applies on the
second axis. Hence the actual captured old mass is

    M_old=e ceil(a)ceil(b).                             (3)

In particular the two axes use the very same radius. Since a<=ceil(a)<a+1,

    L_n=c(r-X)(r-Y) <= M_old
       <= U_n=c(r-X+h)(r-Y+h).                         (4)

The upper inequality is in fact strict in this domain. The weaker displayed
form is convenient for safe tests. To state an exact error identity, set
theta_x=ceil(a)-a and theta_y=ceil(b)-b, both in [0,1). Then

    U_n-M_old
      =e[a(1-theta_y)+b(1-theta_x)+1-theta_x theta_y]
      <=c h(2r-X-Y+h)=:E_n(X,Y,r)
      < epsilon_n=(n²+n-1)/(n²(n+1)).                  (5)

The last strict inequality uses r<1/2 and X,Y>=e. Indeed c h=1/(n+1), so
its uniform upper bound is (1-2e)/(n+1)+e, the stated epsilon_n. This is a
uniform error in captured mass, not an error in deletion area.

## 3. A quadratic necessary test with one shared radius

For a subset J of the actual outside locations in P, including J=empty, put

    B_J=sum_(i in J)b_i,
    R_J(X,Y)=max(X,Y, ||q_i-(-X,-Y)||_infinity : i in J).

The maximum for J=empty is simply max(X,Y). All outside coordinates and
masses remain those of the original input. If a full strict witness captures
old atoms and precisely the outside subset J, (2) gives r>R_J. By (4),

    B_J > F_n(X,Y,r) :=4r²-U_n
        =(4-c)r²+c(X+Y-2h)r-c(X-h)(Y-h),
    r in (R_J,1/2).                                    (6)

Define O_J to be the set of receivers in Gamma for which such an r exists.
An O_J test does not require J to be the exact captured subset: allowing
additional outside sources only enlarges a necessary upper envelope. The
forward proof can always choose the actual captured J, which establishes
the required inclusion. Conversely (6) alone is not an actual halo witness,
since U_n is an upper bound for the discrete old mass.

Let R_nu=H(nu) intersect Gamma, evaluated with the complete residual. Then

    H(mu) intersect Gamma subset R_nu union (union_J O_J),
    Delta intersect Gamma subset (union_J O_J)\R_nu.    (7)

Pure outside witnesses are already in R_nu. The J=empty test represents
pure old witnesses and must not be discarded merely because B_J=0. With two
positive outside locations, R_nu is exactly the union of their two singleton
receivers and pair receiver, cropped by Gamma. For each nonempty label T its
receiver is

    (max_T p_x-r_T,min_T p_x+r_T)
      x(max_T p_y-r_T,min_T p_y+r_T),
    r_T=min(1/2,sqrt(nu(T))/2).

An interval of zero or negative width is empty. Extra sources captured at a
label's common radius only help its true residual witness. This is precisely
the actual-menu argument in `r078-actual-prefix-menu-and-joint-payment.md`;
no residual label or radius is replaced by the quadratic old-count bound.

## 4. Exact convex minimization, ties, and rational cell tests

The leading coefficient A=4-c is positive. Define

    r0=c(2h-X-Y)/(2(4-c)).

If R_J>=1/2, O_J is empty, irrespective of any endpoint polynomial value.
If R_J<1/2, define r_star=clamp(r0,R_J,1/2) and phi_J=F_n(X,Y,r_star).
Continuity and convexity give the exact equivalence

    (X,Y) in O_J  iff  R_J<1/2 and B_J>phi_J.           (8)

Although the minimization uses the closed interval, its strict comparison
is equivalent to existence in the open interval: a strict surplus at an
excluded endpoint persists after a sufficiently small inward perturbation.
If B_J=phi_J, there is no strict witness in the open interval. This includes
an interior minimizer, a boundary minimizer, and a zero surplus. In particular
B_J=0 is accepted exactly when phi_J<0, not when phi_J<=0.

For the actual receiver domain X,Y>=e,

    r0<=c(h-e)/(4-c)=n/(3n²+8n+4)<1/2.

Thus the top clamp is never needed there: r_star=max(R_J,r0). The general
clamp formula remains useful for implementation. If r0>=R_J, phi_J equals

    -c(X-h)(Y-h)-c²(X+Y-2h)²/[4(4-c)];

otherwise it is F_n(X,Y,R_J). Every quantity is rational for rational input
and receiver coordinates. No square root approximation is needed in (8).

For a closed rational receiver cell D=[XL,XU]x[YL,YU] in the positive crop,
write D_Gamma=D intersect (e,1/2)². The cell assertions below concern
D_Gamma; a crop-boundary point is not thereby added to the open corner.
let Rlo=min_D R_J and Rhi=max_D R_J, both computed from the same source bbox
with the origin included. Explicitly, let [a,b]x[d,f] be the bbox of the
reflected locations {0} union {-q_i:i in J}. Then

    Rlo=max((b-a)/2,(f-d)/2,b-XU,XL-a,f-YU,YL-d,0),
    Rhi=max(XU-a,b-XL,YU-d,f-YL,0).

The one-dimensional farthest-distance function is max(z-a,b-z); its minimum
on an interval is the maximum of the half span and the two endpoint
separation terms. Taking the maximum of the independent axis minima gives
Rlo; taking all farthest endpoints gives Rhi. On r>=max(X,Y),

    partial_X F_n=c(r-Y+h)>0,
    partial_Y F_n=c(r-X+h)>0.                           (9)

Consequently the following are safe cell statements:

* D_Gamma is disjoint from O_J if Rlo>=1/2, or, when Rlo<1/2,
  min_[Rlo,1/2] F_n(XL,YL,r)>=B_J.
* D_Gamma is contained in O_J if Rhi<1/2 and
  min_[Rhi,1/2] F_n(XU,YU,r)<B_J.

For exclusion, any actual feasible cell point and radius lies above the
lower-coordinate polynomial on the larger interval. For inclusion, use one
common radius strictly between Rhi and 1/2 after the endpoint perturbation;
(9) makes its surplus valid throughout D. Crucially this inclusion certifies
membership in the safety envelope O_J, not in H(mu). No actual-full inner
area certificate follows from it. The lower bound L_n in (4), or the robust
margin test below, is needed for an actual-full witness.

## 5. Quantitative scope, and improvement over total local mass

At a fixed receiver and common radius, write the soft surplus as
S_soft=B_J+U_n-4r². If r>R_J and S_soft>E_n, then (5) proves
B_J+M_old>4r², so this very radius is an actual full witness. Capturing an
outside atom not listed in J only strengthens the conclusion. In the other
direction every actual old-bearing witness passes (6). Thus the remaining
uncertainty is an explicit mass-margin strip of width at most epsilon_n.
This does not imply that a correspondingly small receiver area is uncertain:
one would additionally need a quantitative control on near-threshold level
sets. In particular no deletion-area continuity or Borel-to-atomic transfer
is asserted here.

For completeness the limiting polynomial at fixed valid X,Y,r is

    F_infinity=3r²+(X+Y)r-XY.

On the whole domain (2) a direct estimate gives

    0<=U_n-(r-X)(r-Y)
       <=(1-c)/4+1/(n+1)+1/n².                        (10)

Indeed U_n=c(r-X)(r-Y)+c h(2r-X-Y)+e, the positive term
c h(2r-X-Y) dominates the possible loss (1-c)(r-X)(r-Y):
if a=r-X,b=r-Y in (0,1/2), then ab<=(a+b)/4 and
c h=1/(n+1)>=(1-c)/4. This also proves the lower sign in (10).
Equation (10) and (5) are uniform state/mass estimates of order 1/n, with
no inference about deletion-area limits.

The quadratic bound can improve the collapsed total-mass proxy strictly.
Use the qualified R77/R78 actual S20 input with

    e=1/400, h=21/400, c=400/441, m_P=1/4,
    q1=(-71/200,0), b1=1/4,
    q2=(0,-171/200), b2=1/2.

At (X,Y)=(2/5,39/100), R_empty=R_{1}=2/5 and
R_{2}=R_{1,2}=93/200. The quadratic is increasing at these radii and above,
and its two minima are

    phi_empty=phi_{1}=16/25-1/336,
    phi_{2}=phi_{1,2}=8649/10000-799/58800.

These exceed respectively 1/4 and 3/4, so every O_J rejects the point.
The complete residual also rejects it: the singleton squared distances are
(39/100)²>1/16 and (93/200)²>1/8, and the pair squared distance
(93/200)²>3/16. The total-local-mass proxy accepts it via the joint label,
since R_{1,2}=93/200<1/2=sqrt(m_P+b1+b2)/2. All comparisons have strict
margins, so the separation also holds on some positive-area neighborhood.
This is a fixed-input local improvement, not a claim of sharpness.

The raw U_n is not everywhere below m_P: already for n=20, X=Y=e and
r=49/100 it exceeds 1/4. To enforce pointwise domination of both mass bounds,
one may use U_clip=min(m_P,U_n), or equivalently the radius test

    B_J>max(4r²-m_P,F_n(X,Y,r)) for one common r.        (11)

This remains safe and is no weaker than either original test. It is no
longer a single quadratic test; separately minimizing its two terms and
combining the resulting inequalities would lose the common-radius condition.

## 6. Actual scalar verification and remaining problem

`computation/round079-kernel/check_quadratic_envelope.py` uses only Python
integer/Fraction arithmetic. Its minimal input is
`computation/inputs/r079-kernel/quadratic-fixtures.json`. It verifies the
corrected S20 raw-source hash and atoms, exact open-prefix counts at rational
events (including integer ties), (4)-(5), convex-minimum strict/equality and
B=0 cases, and the strict comparison point and U_n>m_P example above.
It does not evaluate an area, rebuild full prefix menus, rerun global caps,
or search outside parameters. The analytic argument carries the all-n and
all-receiver claims. Software version, executed command, actual cwd, script
and input hashes, checks and run time are saved in the result receipt.

The original discrete menu remains the exact full-halo representation.
This continuous quadratic envelope is a cheaper necessary test with an
explicit mass error and the original exact residual. A general joint-active
corner area payment and the unrestricted target Delta<=7/2 remain open.
