# R74: a true corner payment with one effective outside atom

2026-11-16. Strict threshold 1, Lebesgue area in R², centered axis-parallel
squares with one common radius. All halo differences below use the actual
residual halo at every radius, not a full-halo proxy or a unit-window slack.

## 1. A one-dimensional overlap inequality

Let I be an interval of length u>0, let 0<=v<=u, and let q be any real number.
Write

    L_s=length(I intersect [q-s,q+s]), s in {v,u}.

Then

    0<=L_v<=L_u<=u,    L_u-L_v<=u-v.                  (1)

Translate I to [0,u]. If q belongs to I, its radius-u interval contains I,
so L_u=u. Its radius-v interval covers length at least v inside I: explicitly
L_v=min(v,q)+min(v,u-q)>=v. If q lies to the left of I, put d=-q>=0. Then
L_u=(u-d)_+ and L_v=(v-d)_+, whose difference is at most u-v. The case q to
the right follows by reflection. These cases prove (1), including endpoints.

For a square A=I_x times I_y with side u and an arbitrary point q in R²,
coordinate products give

    |A intersect (q+[-u,u]²)|
      -|A intersect (q+[-v,v]²)| <=u²-v².               (2)

Indeed put d=u-v. By (1), each inner overlap length is at least
(L_u-d)_+. For 0<=a,b<=u, the function

    ab-(a-d)_+(b-d)_+

is nondecreasing in each variable (verify the cases a<=d, b<=d and both
greater than d). It is at most its value at (u,u), namely u²-v². This proves
(2). It is a bound on the genuine annular intersection with a fixed square;
no independent choice of an inner and outer source measure has occurred.

## 2. The true two-atom payment

Let alpha=m delta_p, beta=M delta_q with m>0, M>0, and arbitrary positions
p,q. Let H denote the centered-square strict halo, and set

    D=H(alpha+beta)\H(beta),
    O_p={x:x_1<p_1, x_2<p_2}.

Then

    |D intersect O_p|<=m/2.                            (3)

Translation allows p=0. Put

    r_m=sqrt(m)/2, v=sqrt(M)/2, u=sqrt(m+M)/2.

For a single atom the complete all-radius halo is exactly its open square
of this radius. For two atoms the complete full halo is the union of the
two singleton receiver squares and the mixed receiver rectangle

    E=(p+(-u,u)²) intersect (q+(-u,u)²).                (4)

To see completeness, a strict witness captures just the first atom, just
the second, or both. In each case its radius is smaller than the corresponding
square-root mass radius. Conversely a point in one of those open receivers
admits a radius slightly larger than its farthest captured distance and
strictly below that square-root radius. Thus (4) is not a chosen-radius
relaxation. It remains valid if a receiver rectangle is empty or the atom
positions coincide.

The beta singleton contributes nothing after subtracting H(beta). The alpha
singleton has area m/4 in O_0. For the mixed part, let A_u=(-u,0)². Since
v<=u, its true residual subtraction has area

    |(E\H(beta)) intersect O_0|
      =|A_u intersect (q+[-u,u]²)|
          -|A_u intersect (q+[-v,v]²)|
      <=u²-v²=m/4,                                    (5)

by (2). Open/closed receiver boundaries are finite unions of line segments
and have zero Lebesgue area; the underlying strict halo description remains
(4). The same outside mass M determines both u and v in (5). Adding the
alpha singleton's possible m/4 proves (3). Their overlap is simply counted
twice in this upper bound, which is harmless. No sharpness of m/2 is claimed.
For M=0 the simpler bound m/4 suffices.

## 3. Collapsing arbitrary northeast internal mass

Now let sigma be any finite positive atomic measure supported in
[p_1,infinity) times [p_2,infinity), of total mass m. Keep the same single
outside atom beta=M delta_q. Then

    |[H(sigma+beta)\H(beta)] intersect O_p|<=m/2.       (6)

For a center x in O_p, if a square centered at x captures any internal source
s in the northeast support, it also captures p: each coordinate distance to p
is no larger than its distance to s. Consequently for every centered square
at such x,

    sigma(Q)<=m delta_p(Q).

This pointwise domination of captured masses on the specified receiver
orthant implies

    H(sigma+beta) intersect O_p
       subset H(m delta_p+beta) intersect O_p.

Subtract the **same** H(beta), then use (3). Collapsing internal mass may
destroy other global constraints, but none is needed for (3), so this is a
valid upper bound. We do not claim domination of measures on arbitrary
squares or any analogous result for several outside atoms.

## 4. Uniform grid: one effective outside atom per corner

For n>=2 use

    e=1/n², h=(n+1)/n², w=1-e,
    sigma_n=(1/n²) sum_(i,j=0)^(n-1) delta_(ih,jh),
    J=[-e,1]², C=[-1/2,w+1/2]².

Let nu be a finite positive atomic outside measure, and assume the full
measure mu=sigma_n+nu satisfies all closed axis-parallel square caps

    mu(Q)<=max{1, side(Q)²}.                            (7)

Let actual deletion be Delta=H(mu)\H(nu). Define

    Gamma_SW=(-1/2,-e)²,
    P_SW=(-1,1/2-e)².

Every strict full witness has radius r<1/2: for r>=1/2, (7) bounds its mass
by 4r². For a center in Gamma_SW, every source in such a witness lies strictly
inside P_SW, since x_j-r>-1 and x_j+r<1/2-e. Thus this localization retains
every possible strict full radius, not only a critical-radius menu.

There are exactly k=floor(n/2) internal coordinates in (-1,1/2-e) on each
axis. Indeed i(n+1)<(n²-2)/2 has precisely the integer solutions
i=0,...,k-1. Hence local internal mass is

    m_SW=k²/n².

Assume nu restricted to P_SW has at most one positive-mass location, denoted
beta=M delta_q (or beta=0). Sources outside P_SW may be arbitrarily numerous.
For a point in Delta intersect Gamma_SW, its strict full witness uses only
sigma_n restricted to P_SW and beta. Also H(beta) is contained in H(nu), so
the point is outside H(beta). Therefore

    Delta intersect Gamma_SW
       subset [H(sigma_n|P_SW+beta)\H(beta)] intersect Gamma_SW.

The local internal sources lie northeast of the origin, while Gamma_SW lies
in its southwest orthant. Applying (6) proves the stronger single-corner
payment

    |Delta intersect Gamma_SW|
        <=floor(n/2)²/(2n²)<=1/8.                      (8)

In particular (8) implies the proposed 1/8+3/(4n⁴) corner bound in this class.
For n=3 it gives 1/18. Other outside sources were not discarded from the true
residual halo: replacing it by H(beta) only reduces residual credit and is
used in the safe direction for an upper bound.

Reflect Gamma_SW and P_SW independently about the grid midlines x=w/2 and
y=w/2 to obtain the other three corner receiver regions and effective source
windows. On either axis, the source intervals are

    (-1,1/2-e),       (1/2,w+1),

which are disjoint. Thus the four effective windows are pairwise disjoint.
If each contains at most one positive outside location, (8) applies to each
using the **same actual global nu and caps**. The corner receiver regions are
also disjoint, so

    sum_(four corners) |Delta intersect Gamma|
       <=2floor(n/2)²/n²<=1/2.                         (9)

No assumption that square caps act independently on the four source windows
is made. A cap crossing several windows is part of the original hypothesis.

## 5. A complete 3.5 payment for this structural class

Every actual deletion point has a strict full witness capturing an internal
source; otherwise that same witness would put it in H(nu). With r<1/2 this
places Delta inside C. The four corner receiver regions have total area
(1-2e)², while |C|=(2-e)². Their complement in C therefore has area

    (2-e)²-(1-2e)²=3-3e².

Combining this elementary full container bound with the **true deletion**
corner bound (9) gives

    |Delta|<=3-3/n⁴+2floor(n/2)²/n²
             <=7/2-3/n⁴<7/2.                          (10)

This proves a general structural class, not just a fixed input: any n>=2,
any finite outside positions and masses satisfying (7), and at most one
positive outside location in each of the four effective corner windows.
Outside those windows the measure is unrestricted apart from (7). We do
not require S_n to be the unique maximum block. The case of multiple effective
outside atoms in a corner remains unresolved; neither (6) nor (10) is claimed
for it. In particular the argument does not settle the unrestricted uniform-
grid target or the general sharp weak constant.

## 6. Actual exact checks and the next obstruction

The minimal fixture is
`computation/inputs/r074/two-atom-corner-payment.json`, with runner
`computation/round074/check_two_atom_corner_payment.py`. The actual Fraction
run checked 255 one-dimensional and 4335 two-dimensional rational geometric
cases. These include v=0, v=u, boundary contacts, disjoint intervals and both
signs of the outside center, and retain the complete residual square of
radius v. They are regression checks of the continuum proof above, not a
finite-grid proof or a candidate search. Equation (2)'s annulus bound is
attained in some regressions; no sharpness of the complete m/2 bound follows.

The same run checks a genuine n=3 input: all nine grid atoms of mass 1/9 plus
one outside atom at (-2/5,-2/5), of mass 1/4. Literal coordinate-box membership
enumerates all 46 nonempty closures, with 28 eligible closures. It confirms
all square-cap inequalities, the unique maximum internal block, one eligible
component, and positive nonterminal margin 29891/89100. P_SW contains just
the internal origin, so (8) certifies its corner deletion is at most 1/18.
No complete halo area was computed or inferred from the regressions.

The actual run used Python 3.15.0 standard library, exit 0, about 0.058 seconds
internally. Input/script/version hashes and results are saved in
`computation/results/round074-two-atom/payment-checks.json`. From the topic
directory:

    python3 -B computation/round074/check_two_atom_corner_payment.py \
      --input computation/inputs/r074/two-atom-corner-payment.json \
      --output /Volumes/data/CodexTemp/T2/r074-two-atom-recheck.json

The next unhandled case is genuinely multiple effective outside labels in
one corner. Applying (5) separately to every such label would charge the
same internal mass repeatedly and is not justified. A useful continuation
must couple those labels through their actual residual union, or prove a
local atomic residual-area credit. Single-atom location refinement cannot
produce a counterexample to (8), so no such scan was run. This component
creates no subagents, runs no MATLAB, and makes no commit or push.
