# R73: shared unit-window integral kernels and a non-CDF obstruction

2026-11-16. Threshold 1, centered axis-parallel common-radius squares. Let
sigma_n be the uniform n-by-n grid of mass 1, n>=2, with

    e=1/n², h=(n+1)/n², w=1-e,
    S_n={0,h,...,(n-1)h}²,
    J=[-e,1]², C=[-1/2,w+1/2]².

Assume the same nonnegative finite outside measure nu satisfies all closed
square caps for mu=sigma_n+nu:

    mu(Q)<=max{1, side(Q)²}.                            (1)

For finite atomic inputs, the separate R64/R72 result identifies this with
S_n remaining a maximum eligible mass-1 block. All formulas below keep one
common nu; edge and corner masses are never selected separately by receiver.

## 1. The zero collar and the exact doubly shifted L

The four closed unit squares whose coordinate intervals are chosen from
[-e,w] and [0,1] each contain S_n. Each has sigma mass 1, hence outside mass
zero by (1). Their union is J, so

    nu(J)=0, including its boundary.                    (2)

For e<=a,b<=1 define

    Q_(a,b)=[-a,1-a] x [-b,1-b],
    k(a)=#{i in {0,...,n-1}: ih<=1-a}.

Its exact unit cap is

    nu(Q_(a,b)) <= 1-k(a)k(b)/n².                      (3)

Using (2), its outside part can be written as the disjoint L-shaped union

    [-a,-e) x [-b,1-b]
       union [-e,1-a] x [-b,-e).                       (4)

The lower-left corner belongs to the first piece and is not counted twice.
The half-open conventions in (4) retain all possible atoms outside J; only
points in J have been removed using their proved zero mass. Increasing a
or b moves the opposite end of the unit window as well. Thus these L sets
are not a nested two-parameter family. The actual function nu(Q_(a,b)) need
not even be coordinatewise increasing; (3) cannot be treated as a two-variable
cumulative distribution constraint.

## 2. A general positive continuous dual kernel

Write W_z=z+[-1/2,1/2]² for a closed unit square. For any nonnegative
integrable measurable weight f on receiver centers, define

    K_f(q)=integral f(z) 1_(q in W_z) dz.

Integrating the valid unit caps against the **same** f and using Tonelli gives

    integral K_f(q) dnu(q)
       <= integral f(z) dz - (1/n²) sum_(p in S_n) K_f(p).  (5)

This is an explicit all-position nonnegative dual, valid for finite atomic
and finite Borel outside measures. No candidate source grid is involved.
Closed/open window boundaries in the z integration are null for each fixed
source q, so the kernel identity does not discard atomic nu mass on any
source-space boundary. If preferred, the general window parameter can carry
any nonnegative finite measure rather than f(z)dz; the corresponding kernel
is its exact window-incidence integral.

For f=1_C, the kernel has the closed form

    K_C(q)=kappa(q_x) kappa(q_y),
    kappa(s)=0,             s outside [-1,w+1];
             1+s,          -1<=s<=0;
             1,             0<=s<=w;
             1+w-s,         w<=s<=w+1.

Indeed kappa(s)=length([-1/2,w+1/2] intersect [s-1/2,s+1/2]). Every internal
grid coordinate lies in [0,w], so K_C(p)=1 for every p in S_n. Formula (5)
therefore becomes the exact continuous weighted edge/corner budget

    integral kappa(q_x)kappa(q_y) dnu(q)
       <= |C|-1 = 3-4e+e².                            (6)

By (2), the source integral may be restricted to the complement of J. The
kernel is positive on the interior of [-1,w+1]² and weights external sides
and corners by their actual incidence with unit windows. It is not a sum
of independently optimized outside masses. This is a mass budget, not yet
an area payment for actual deletion or its missing outer ring.

There is also a useful accounting distinction for four reflected corner
families. Integrate Q_(a,b), a,b in [0,1], and its three reflections about
the grid midlines. In each family the one-axis average internal count divided
by n is (1+e)/2. Their summed right-hand side is

    4[1-((1+e)/2)²]=3-2e-e².                          (7)

The corresponding center weight is the product of the indicators of
[-1/2,1/2] and [w-1/2,w+1/2], summed on each axis. These two intervals overlap
in a strip of width e; hence this weight overcounts those receiver centers
relative to 1_C. The difference between (7) and (6) is exactly

    2e(1-e)>0.

Summing corner-window families is valid, but failing to account for their
overlap loses this amount of budget. It does not manufacture additional
independent choices of nu.

## 3. A legal finite atomic input violating the CDF interpretation

Take n=5 and add one outside atom

    q=(-1/5,4/5),   nu({q})=1/100.

This is a genuine admissible input, with all masses positive. S25 has mass 1
and span 24/25, so it is eligible. Any proper subset of S25 together with q
has mass at most 24/25+1/100=97/100<1. The whole 26-atom input has mass 101/100
and span 29/25; its span squared 841/625 exceeds its mass. Thus S25 is the
unique maximum eligible block, and every square cap (1) holds.

The pair consisting of q and the internal atom (0,18/25) is eligible: its
span squared is 1/25 and its mass is 1/20. Since S25 itself is eligible, the
eligible-source graph is connected. Its nonterminal margin is strictly

    (29/25)²-(35/33)(101/100)=22637/82500>0.

Nevertheless

    nu(Q_(3/10,1/10))=1/100,
    nu(Q_(2/5,3/10))=0,

although both shift parameters increase. The first square is
[-3/10,7/10]x[-1/10,9/10] and the second is
[-2/5,3/5]x[-3/10,7/10]. Their internal masses are 12/25 and 9/25; their full
unit-cap slacks are 51/100 and 16/25. Thus the nonmonotonicity is present in
the actual shared measure of a connected nonterminal unique-maximum input,
not merely in arbitrary unrelated L sets. It refutes the CDF substitution,
not any genuine unit-cap inequality or the target 3.5 deletion bound.

## 4. The uniform integral has no uniform positive slack in the Borel class

Use the actual R72 Borel family

    nu_n=t 1_(B\J) dx,   t=n/(n+1), B=[-1,2]².

The separate R72 proof establishes all square caps (1). Its final source is
`proofs/r072-uniform-grid-borel-family-and-atomic-discontinuity.md`, SHA-256
`7b49ddee130e5deee6bbf17641554fbfe122dae9925f73c865574e7073eddf53`.
The support of K_C lies in B. Direct one-dimensional integration gives

    integral kappa(s) ds = 2-e,
    integral_[-e,1] kappa(s) ds = 1+e-e².

Hence the actual left side of (6) is

    t[(2-e)²-(1+e-e²)²].                               (8)

Its nonnegative slack is exactly

    S_n=[(3-4e+e²)+n((1+e-e²)²-1)]/(n+1).              (9)

For 0<e<=1/4, (3-4e+e²)<=3 and
(1+e-e²)²-1=2e-e²-2e³+e⁴<=2e. Consequently

    0<=S_n<=(3+2/n)/(n+1)<=4/(n+1),
    n S_n -> 3.

Therefore there is no n-independent positive lower bound on the slack of
this integrated unit-cap budget over the admissible Borel class. A strategy
that first demands such a fixed positive unit-cap integral deficit cannot
pay the target outer-area gap |C|-7/2, which tends to 1/2. This is a precise
obstruction to that intermediate demand; it does not say that (6), finer
weights, common-radius information or actual residual-halo credit cannot
participate in a successful proof.

No finite-atomic area conclusion is inferred by approximating this diffuse
measure: R72 proves actual deletion can be discontinuous. The next section
instead proves a uniform **unit-window mass** comparison, which transfers
the slack obstruction to actual finite inputs in the remaining domain. It
does not transfer deletion area.

## 5. Uniform slack obstruction in the unique, connected, nonterminal finite class

For every n>=5 set

    delta=delta_n=1/(100n²), lambda=1-delta²,
    B_*=[-1,3] x [-1,2], J_delta=J+[-delta,delta]²,
    c=t/(1+2delta).

Take every closed delta-lattice cell contained in B_* and disjoint from the
closed set J_delta. At its center put mass lambda c delta²; call the resulting
finite outside measure eta_n. This definition is an implicit finite index
formula, not a request to enumerate its very large atom list. The exact-cap
atomization theorem of Section 2 of
`proofs/r073-cap-preserving-atomization-and-local-full-halo-obstruction.md`
applies because delta<=1/(2n), and its proof permits any bounded truncation
B_*. Thus sigma_n plus the uncontracted outside masses c delta² satisfies all
square caps (1); contracting outside masses by lambda preserves those caps.

S_n remains eligible of mass 1. For any subset eligible after contraction,
its uncontracted version is also eligible and has mass at most 1. If it
contains outside mass, its contracted mass is strictly below 1. A proper
internal subset also has mass below 1. Therefore S_n is the unique maximum
eligible block.

Here are explicit checks of the remaining domain properties. Put m=1/delta=
100n². The full lattice rectangle has 4m by 3m cells. Because e=100delta,
the cells intersecting closed J_delta have each axis index, relative to the
origin -1, in the single integer interval

    [m-102,2m+1].

It is strictly inside the outer lattice rectangle. Its complement is a
four-neighbor connected rectangular grid ring. For n>=5 and delta<=1/100,
lambda c>=66/85>1/2 (use t>=4/5, lambda>=99/100 and 1+2delta<=51/50).
Every neighboring pair of outside atoms is consequently eligible, since its
mass 2lambda c delta² exceeds its span squared delta². The selected cell

    [-e-3delta,-e-2delta] x [0,delta]

has center q=(-e-5delta/2,delta/2). Its distance to the internal source (0,0)
is e+5delta/2=41/(40n²)<1/n. Their pair is eligible even from the internal
mass alone. This joins the outside ring to S_n, whose whole mass-1 label is
eligible. The complete eligible-source graph is therefore connected.

The source span is exactly 4-delta, realized by the leftmost and rightmost
outside centers. Total mass is at most 1+12lambda c<=13. Hence its nonterminal
margin is bounded below uniformly by

    (399/100)²-(35/33)*13=703633/330000>0.              (10)

The rectangle B_* is essential here: using the old square truncation and
letting n tend to infinity would not preserve this nonterminal conclusion.

We next prove a uniform comparison for **all closed unit windows**, including
boundary ties. Let E=B_*\J and let E_delta be the union of the selected
lattice cells before replacing their mass by center atoms. B_* is lattice
aligned, and every missing cell outside J_delta lies within distance delta
of J_delta. Consequently, up to zero-area boundaries,

    E_delta subset E,
    E\E_delta subset (J+[-2delta,2delta]²)\J.

For n>=2 and delta<=1/4 the latter area is at most

    8(1+e)delta+16delta² <=14delta.                    (11)

For a closed unit square Q, the union of selected cells whose centers lie
in Q differs from Q intersect E_delta only inside the centered boundary strip
between square sides 1-delta and 1+delta. Its area is 4delta, even when a cell
center is exactly on a Q edge. As lambda c<=1,

    |eta_n(Q)-lambda c |Q intersect E_delta||<=4delta.

Also t-lambda c<=2delta+delta²<=3delta. Comparing to the diffuse measure
nu_n^*=t 1_E dx and using unit area and (11) therefore gives

    |eta_n(Q)-nu_n^*(Q)|<=4delta+14delta+3delta
                              =21delta.               (12)

This is a uniform closed-window mass estimate, not a weak-convergence claim
for discontinuous indicators. It is independent of the number of atoms.

To control the Borel window deficits, let F_n be the one-dimensional grid
empirical CDF and F the CDF of Lebesgue measure on [0,1]. On a gap
[kh,(k+1)h), F_n=(k+1)/n. Its difference from F lies between
-(k+1)/n² and (n-k)/n²; left limits obey the same 1/n bound. On [w,1] the
difference is at most e. A unit interval cannot have both endpoints strictly
inside (0,1). Its closed grid count k/n therefore differs from its overlap
with [0,1] by at most 1/n, including endpoint ties. Replacing [0,1] by [-e,1]
changes overlap by at most e. If the two axis overlaps with J are L_x,L_y
and the grid counts are k,l, all quantities are at most 1, so

    L_x L_y-kl/n² <=2(1/n+e).

For every center z in C_infinity=[-1/2,3/2]², its unit square W_z lies inside
the original B=[-1,2]² and B_*. Thus its full Borel mass is exactly
kl/n²+t(1-L_x L_y), and its nonnegative unit-cap deficit obeys

    0<=1-(sigma_n+nu_n^*)(W_z)
       <=1/(n+1)+2/n+2/n² <=4/n.                      (13)

Combining (12), (13) and the finite exact square caps proves the final
pointwise theorem:

    0<=g_n(z):=1-(sigma_n+eta_n)(W_z)
       <=4/n+21/(100n²), for EVERY z in C_infinity.     (14)

For any finite positive measure rho_n on window centers supported in
C_infinity, define K_(rho_n)(q)=integral 1_(q in W_z) d rho_n(z). Tonelli,
using the one common finite outside measure eta_n, now gives

    D_n(rho_n):=rho_n(C_infinity)
        -(1/n²) sum_(p in S_n) K_(rho_n)(p)
        -integral K_(rho_n)(q) d eta_n(q)
      =integral g_n(z) d rho_n(z),
    0<=D_n(rho_n)
      <=[4/n+21/(100n²)] rho_n(C_infinity).             (15)

Thus for every fixed M<infinity, all choices rho_n(C_infinity)<=M have
D_n(rho_n)->0 simultaneously. The weights may depend arbitrarily on n and
the input. This includes L1 density weights and atomic window-parameter
weights; no transfer of their possibly discontinuous kernels is needed,
because (14) holds at every window center before integration.

In particular there is no uniform positive unit-cap slack obtainable by
choosing bounded-total-mass positive window weights, even within finite
inputs having unique S_n, one eligible component, and strictly positive
nonterminal margin. To demand slack at least a fixed a>0 along this family
would require weight mass at least a/[4/n+21/(100n²)]. This is an obstruction
to that specific bounded-weight intermediate goal, not to every possible
use of unit caps. It does not control deletion area or prove its target bound.

## 6. Actual exact checks and remaining task

`computation/round073/shared_unit_window_kernel.py` reads the minimal fixture
`computation/inputs/r073/noncdf-one-outside.json`. It independently enumerates
all 271 nonempty bbox closures of the 26-atom input by literal coordinate-box
membership, finds 124 eligible closures, verifies every cap, the unique
maximum S25, one component and the exact nonterminal margin. It checks the
two actual unit windows, verifies K_C at 144 rational test positions including
every kernel kink and the outside atom's coordinates, and independently
integrates kappa by exact trapezoids. For the n=5 Borel fixture it obtains

    |C|-1=1776/625,
    S_5=42096/78125,
    integral K_C dnu_one_atom=1/125.

The same finite n=5 run also partitions the whole continuous unit-window
center domain C_infinity at every source entry/exit and overlap kink. It
checks all 441 point/open-interval product cells, evaluating all 961 relevant
corner limits with fixed cell counts. This certifies the n=5 deficit enclosure
[0,9/25] and checks the bound 97/150 from the first expression in (13), not
just sampled centers. The saved `maximum_deficit` is this corner-limit upper
bound; no attainment of its endpoint by an actual window center is claimed.
Closed ties and one-sided cell limits are retained.

Finally it checks the scalar geometry of the implicit n=5 atomization with
m=2500: its 68219184 outside atoms are specified by index formula and **not
enumerated**. It checks the forbidden index interval, the bridge distance,
positive neighbor-pair margin, coarse nonterminal bound and constants in
(11)-(12). The general square caps and all-n construction are justified by
the separate atomization theorem and the proof above, not inferred from this
scalar fixture. No large atom array or new halo calculation is created.

The Python 3.15.0 Fraction run actually passed, exit 0, about 0.038 seconds.
Exact output and
input/script/version receipt are in
`computation/results/round073-kernel/shared-unit-window-kernel.json`.
Only this fixed fixture was computed; the all-n identities, bound and limit
are proved above rather than inferred from a numerical scan. No MATLAB,
full halo geometry or new Borel area scan was run.

The missing step remains an inequality converting compatible unit-window
mass information and residual-halo exclusion into a fixed positive amount of
actual missing outer-ring area. Neither (6) alone nor the rejected CDF
substitution supplies that step. The target actual deletion bound 7/2 and
the general sharp weak constant remain unresolved.
