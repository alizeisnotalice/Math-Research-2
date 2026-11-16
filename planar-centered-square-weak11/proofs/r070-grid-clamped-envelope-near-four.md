# R70: internal-capture clamps and all remaining prefixes can tend to four

Independent derivation, 2026-11-15. Finite distinct positive atomic sources,
strict centered axis-parallel square receivers, threshold one. All areas are
planar Lebesgue areas. This concerns a necessary-condition outer envelope;
it gives no lower bound for an actual mixed halo or deletion area.

## 1. The envelope and an explicit sub-envelope

Normalize the maximum eligible source-block mass to L=1. For nonempty
proper U subset S set a=m(U), r_U(x)=max_(p in U)||x-p||_infinity, and
c_S(x)=max_(p in S)||x-p||_infinity. A mixed label with internal set U and
total mass m must satisfy

    a<m<=1, 4r_U(x)^2<m,
    c_S(x)>f_a(m):=sqrt(1+m-a)-sqrt(m)/2.              (1)

Let F(S) be the union of Halo(S) and all x satisfying (1) for at least
one U and at least one real m. This is the full S-only internal-capture
clamp, with the outside mass treated as a free scalar. Allowing equality
at minimization endpoints gives the same interior carrier, since a strict
gap in the last inequality and r_U<1/2 permit an upward perturbation of m
when needed. No actual outside sources or outside-state constraints enter
this definition.

Put

    U_x=S intersection (x+(-1/2,1/2)^2), A_S(x)=m(U_x),
    E(S)={x:0<A_S(x)<1,
              c_S(x)>sqrt(2-A_S(x))-1/2}.            (2)

Then E(S) subset F(S): U_x is nonempty and proper, every one of its
finitely many points has distance strictly below 1/2, and m=1 is feasible
in (1). In fact one may choose m<1 sufficiently close to 1, by continuity.
Thus the inclusion still holds if mass-one mixed ties are excluded. U_x
is a bbox closure within S: every S-source in bbox(U_x) lies in the same
open centered square. The construction respects an actual common-radius
internal capture; it does not assert an actual mixed source state exists.

## 2. A sequence with |E(S_n)| tending to 4

For odd n>=3 let w_n=1-1/n^2 and place n^2 old sources at

    p_ij=(i*w_n/(n-1),j*w_n/(n-1)), 0<=i,j<n,
    m(p_ij)=1/n^2.                                    (3)

The block S_n has mass 1 and bbox side w_n<1, hence is strictly eligible.
Its side-one source squares cover the entire square
C_n=(-1/2,w_n+1/2)^2. In particular F(S_n) subset C_n, so |F(S_n)|<4.

The empirical source measures converge to uniform area measure on [0,1]^2.
For each fixed x, the boundary of the open side-one square centered at x
has zero limiting area, so A_(S_n)(x) converges to

    A(x)=length([0,1] intersection (x_1-1/2,x_1+1/2))
         *length([0,1] intersection (x_2-1/2,x_2+1/2)).

Also c_(S_n)(x) converges to c(x)=max(|x_1|,|x_1-1|,|x_2|,|x_2-1|).
For x in [0,1]^2 write delta=max(|x_1-1/2|,|x_2-1/2|). If delta>0,

    c=1/2+delta, A>0, A<1, A>=(1-delta)^2,
    (c+1/2)^2-(2-A)>=2delta^2>0.                      (4)

Thus c>sqrt(2-A)-1/2. For x in C=(-1/2,3/2)^2 outside [0,1]^2,
c>1 and 0<A<1, whereas sqrt(2-A)-1/2<sqrt(2)-1/2<1.
Consequently every x in C except its center satisfies the strict limit
conditions in (2). Convergence and strictness imply x in E(S_n) for all
sufficiently large n. Fatou's lemma, applied on the bounded square C,
now gives liminf |E(S_n)|>=|C|=4. Together with

    E(S_n) subset F(S_n) subset C_n, |C_n|=(1+w_n)^2 ->4,

this proves

    |E(S_n)| ->4 and |F(S_n)| ->4.                    (5)

This proof does not enumerate the 2^(n^2) possible internal subsets and
does not identify a realizable mixed label for each point of the envelope.
It proves that these necessary S-only conditions alone cannot give a
uniform bound below 4, even with their full internal-radius clamp.

## 3. A legal ambient input with a unique maximum block

Adjoin the single source

    q_n=(33/32,w_n/2), m(q_n)=epsilon=1/1000.          (6)

All masses are positive and all positions distinct. We verify the ambient
conditions without an exhaustive subset claim. Any q_n-containing source
subset T of mass at least 1 contains old sources. If its leftmost old
column has index i, then its old mass is at most (n-i)/n. Hence

    m(T)>=1 implies i/n<=epsilon,
    min old x =i*w_n/(n-1)<=epsilon*n*w_n/(n-1)<=2epsilon.

Its horizontal span is therefore at least 33/32-2epsilon, and

    (33/32-2epsilon)^2>1+epsilon>=m(T).

It cannot be eligible. Every subset not containing q_n has mass at most 1,
with equality only for the entire S_n. Thus S_n is the unique maximum
eligible block of the full input, not just of the old sources.

The input is connected in the eligible-source incidence sense: S_n itself
is eligible, and q_n forms an eligible pair with (w_n,w_n/2), which exists
because n is odd. Indeed

    epsilon+1/n^2-(1/32+1/n^2)^2
      =epsilon-1/1024+15/(16n^2)-1/n^4>0, n>=3.

The full input has bbox span 33/32 and mass 1001/1000. Its nonterminal
margin, in the round-53 terminal criterion, is exactly

    (33/32)^2-(35/33)*(1001/1000)=1529/844800>0.

Both S_n widths exceed 3/4. These are fat blocks of rank n^2>=9. Their
side-one support-square area is (1+w_n)^2. They also fail the earlier
support-tree sufficient test: writing h=w_n/(n-1), every nonzero edge
between grid sources has cost dx+dy-dx*dy>=h, because 0<=dx,dy<1 and
each nonzero coordinate increment is at least h. A nearest-neighbor grid
Hamiltonian path attains this lower bound on every edge. Thus the exact
minimum tree cost is (n^2-1)h=(n+1)w_n>5/2.

The construction does not claim failure of every heavy-core certificate
or of the actual 3.5 deletion selection problem. Section 3 alone makes no
all-prefix area claim; section 5 proves the strengthened envelope
obstruction. The ambient input is a legal domain embedding for these
S-only envelope obstructions.

## 4. Independent finite exact checks

The standalone standard-library script
`computation/review_r070_grid_clamp_obstruction.py` reads
`computation/inputs/r070-independent-grid-clamp.json`. It imports no
producer or prior-round geometry code. For n=3,5,9,17,33,65 it partitions
C_n at every source coordinate plus/minus 1/2. A_S is constant on each
open rectangular cell. For each proper positive mass A, failure of (2)
is precisely the intersection with

    K_h=[w_n-h,h]^2, h=sqrt(2-A)-1/2.

Integer-isqrt bounds at denominator 2^80 enclose sqrt(2-A). The lower
area subtracts the intersection with K_(h_upper); the upper area uses
K_(h_lower). Fractions carry all coordinates, areas and domain margins.
Zero-area cell boundaries do not affect the result. All full-mass cells
are excluded from E by definition. The computation checks exact scalar
margins in section 3; it does not enumerate all subsets of the large
ambient input.

The actual finite lower areas (decimals for display only) are
3.242330569895 for n=3, 3.707650300417 for n=5,
3.874437996709 for n=9, 3.952014089891 for n=17,
3.981712671134 for n=33, and 3.993218376951 for n=65.
The archived result retains the exact rational endpoints and script/input
hashes. Sections 1-3 supply the general theorem; the finite checks are
independent, directed finite certificates.

Portable command from the topic archive (output must be chosen explicitly):

    python3 -B computation/review_r070_grid_clamp_obstruction.py --input computation/inputs/r070-independent-grid-clamp.json --output /Volumes/data/CodexTemp/T2/r070-grid-replay/result.json

No MATLAB process, full halo, LP, random search or source-area claim is
part of this calculation. Root owns publication and remote verification.

## 5. The all-remaining-subset strengthening also tends to four

This extension was proposed jointly by root and core and independently
derived here. It retains more internal geometry than sections 1-2.

For a fixed proposed internal capture U, total mass m and receiver x, every
nonempty H subset S\U with M=m(H) and M+m>1 must satisfy

    c_H(x)>sqrt(M+m)-sqrt(m)/2.                        (7)

Otherwise H union T would be eligible with mass M+m>1. Its H-H spans are
strictly below 1<sqrt(M+m), its T-T spans strictly below sqrt(m), and its
cross spans are strictly below sqrt(M+m) by the strict receiver radius.
H need not be eligible. For fixed x, any failed H can be enlarged to all
remaining sources of distance at most c_H(x). This preserves c_H and
increases M, so also fails. Thus it suffices to check the at most |S\U|
nonempty distance prefixes, each including all ties. A zero-distance
prefix must not be omitted.

Define G(S) as Halo(S) union the set of x admitting a nonempty proper
U subset S and one common m satisfying (1) AND (7) for every such H.
The pure internal halo is retained without a mixed-capture requirement.
Then, for the grids (3),

    |G(S_n)| ->4.                                     (8)

Here U and m may depend on x. The full ambient input (6) still provides
the legal unique-max/fat/connected/nonterminal domain embedding; G(S_n)
does not encode that input's actual outside masses or positions.

To prove (8), fix x in C=(-1/2,3/2)^2 with
delta_i=|x_i-1/2| in (0,1), i=1,2. Excluding the two coordinate midlines
costs zero area. Put

    epsilon_n=n^(-1/2), m_n=1-epsilon_n,
    r_n=sqrt(m_n)/2,
    U_n=S_n intersection (x+(-r_n,r_n)^2), a_n=m(U_n).

For all large n, U_n is nonempty, proper, a_n<m_n and r_(U_n)(x)<r_n.
Indeed a_n converges to A(x) in section 2, which lies strictly between
0 and 1. The use of moving radii is justified by the uniform discrepancy
bound below, together with continuity of the limiting rectangle area.

Write A_x(R) for the continuous uniform mass in the closed radius-R
square centered at x. Eventually r_n>|delta_i-1/2| for both coordinates.
For every R>=r_n we then have the exact expression

    A_x(R)=product_(i=1,2) min(1,R+1/2-delta_i).       (9)

The positivity and boundary-crossing condition preceding (9) matters:
the expression need not hold at arbitrary small R. Its almost-everywhere
derivative obeys both A_x'(R)<=2 and
A_x'(R)<=2R+1-delta_1-delta_2. Hence, writing r=r_n and
M_cont(R)=A_x(R)-A_x(r),

    M_cont(R)<=2(R-r),
    M_cont(R)<=(R-r)(R+r+1-delta_1-delta_2).           (10)

For completeness, the needed uniform discrepancy is elementary. The 1D
grid positions are i*w_n/(n-1). For 0<=t<=w_n its open or closed empirical
CDF lies between t*n/(n+1) and t*n/(n+1)+1/n, so its discrepancy from the
uniform [0,1] CDF is at most 1/n. The same bound holds outside [0,w_n].
Taking the difference of two CDFs gives discrepancy at most 2/n for any
open or closed interval. The Cartesian product has discrepancy at most
4/n for every open or closed axis-parallel rectangle. We use the weaker
8/n bound for each rectangle below.

For any distance prefix of S_n\U_n, let R be its farthest distance and
M_n(R) its mass. Necessarily R>=r. Open capture and closed prefix give
exactly

    M_n(R)=empirical_mass(closed Q(x,R))
             -empirical_mass(open Q(x,r)),
    |M_n(R)-M_cont(R)|<=16/n.                         (11)

Only active prefixes M_n(R)>epsilon_n need (7). By (10)-(11), such a
prefix has R-r>(epsilon_n-16/n)/2. The cap corresponding to (7) is

    cap(R)=(R+r)^2-m_n=(R-r)(R+3r).

With gamma_n=2r-1+delta_1+delta_2, equations (10)-(11) yield

    cap(R)-M_n(R)>=(R-r)gamma_n-16/n.

Eventually gamma_n>=(delta_1+delta_2)/2>0. Uniformly over ALL active
prefixes, therefore,

    cap(R)-M_n(R)
      >(delta_1+delta_2)(epsilon_n-16/n)/4-16/n>0      (12)

for all sufficiently large n. The last inequality follows from
epsilon_n=n^(-1/2), for the fixed positive delta_1+delta_2. Thus every
active prefix simultaneously has M_n(R)<cap(R), exactly the strict
inequality R>sqrt(M_n(R)+m_n)-r. All other H are covered by prefix
enlargement. In particular H=S_n\U_n is active because a_n<m_n, and
its c_H<=c_S proves the original condition (1) as well. This constructs
a witness (U_n,m_n) for x in G(S_n), for all large n.

Fatou's lemma and G(S_n) subset C_n finish (8), as in section 2. There
is no finite-n all-prefix area certificate asserted in section 4: those
saved finite areas concern E(S_n), the fixed-total sub-envelope only.
The all-prefix result (8) is the analytic uniform-error theorem proved
in this section. It still leaves actual outside-source compatibility,
and actual W/Wprime or deletion budgets, unresolved.

## 6. One explicit finite all-prefix area certificate

The following finite instance uses the uniform error estimates of section 5
directly. It does not enumerate the large source grid. Set

    n=400001, w=1-1/n^2, m=99/100, r=sqrt(99/100)/2,
    D=(-49/100,49/100) union (51/100,149/100),
    Omega=D cross D.

Use exactly the grid S_n in (3), with n^2=160000800001 positive sources,
each of mass 1/n^2. Then

    Omega subset G(S_n), |Omega|=(49/25)^2=2401/625>7/2.       (13)

For x in Omega, each delta_i=|x_i-1/2| lies strictly between 1/100
and 99/100. Since sqrt(99/100)>99/100, r>99/200>49/100,
and therefore r>|delta_i-1/2|. Formula (9) holds for every R>=r.
Its two factors at R=r are greater than 1/200 and less than 99/100.
For the open capture U=S_n intersection (x+(-r,r)^2), write a=m(U).
The section-5 open-rectangle discrepancy gives

    a>1/40000-8/n>0,
    a<(99/100)^2+8/n<99/100=m.                              (14)

Thus U is nonempty and proper, a<m<1, and its actual maximum distance
is strictly less than r. Also

    gamma=2r-1+delta_1+delta_2>1/100.

For every remaining complete distance prefix of mass M>1/100, (10)-(11)
give R-r>(1/100-16/n)/2>0 and

    (R-r)(R+3r)-M
      >=(R-r)gamma-16/n
      >(1/100-16/n)/200-16/n
       =1/20000-402/(25n)>0.                               (15)

Every inequality is uniform in x in Omega and in the active prefix.
Prefixes of mass M<=1/100 satisfy the inactive alternative M+m<=1.
The same fixed m=99/100 therefore works for all prefixes, and prefix
enlargement gives (7) for all H. The full complement S_n\U is active
because a<m; applying its inequality and c_H<=c_S also gives (1).
This proves Omega subset G(S_n). There is no incompatibility with
G(S_n) subset C_n: the right edge of Omega is below 149/100,
while w+1/2>149/100 because 1/n^2<1/100; its left edge exceeds -1/2.

The exact scalar checker `computation/review_r070_finite_prefix_area.py`
with input `computation/inputs/r070-finite-prefix-area.json` was actually
run using Fraction arithmetic. It checks (14)-(15), the squared radical
comparison, the area and C_n containment margins. This is an analytic
certificate for a finite formula-specified input, not enumeration of its
160000800001 sources or of any subsets, and not a calculation of the
actual halo, W/Wprime or deletion area. The ambient extension of section 3
applies unchanged. No general bound on C2 is changed.
