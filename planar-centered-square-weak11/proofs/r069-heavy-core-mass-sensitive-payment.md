# R69: a mass-sensitive heavy-core payment class

2026-11-15. Finite positive atomic inputs, threshold1, centered axis-parallel
common-radius squares, strict witness inequality. This theorem gives a new
sufficient class for the full erosion-union W and hence Wprime. It does not
prove the unrestricted3.5 maximum-block selection proposition.

## 1. Definitions and an explicit finite payment functional

Let S be a maximum-mass eligible source block, L=m(S)>0 and s=sqrt(L).
Its bbox widths w,h satisfy w,h<s. Define

    C = bbox(S)+[-s/2,s/2]^2,
    Z_s(H)=union_(p in H) (p+(-s/2,s/2)^2).

Choose ANY nonempty core H subset S and put

    M=m(H)>0, E=L-M>=0, V=S\H, delta=sqrt(E),
    B={q in A\S:dist_infinity(q,V)<delta}, beta=m(B),
    nu=the original atomic measure restricted to V union B.

If V is empty, define B empty. The original source IDs are counted once,
so nu has total mass exactly E+beta. B is the entire indicated finite near-
outside set, not a selected label or a sum of reused captured masses.

For a radius r>=0 satisfying4r^2<=3M, set

    K_H(r)=intersection_(p in H) (p+[-r,r]^2).

K can be empty. No assumption H is eligible or K subset C is made. Define

    Gamma_r(H)=|C\K|+|Z_s(H) intersection K|
              =|C|-|C intersection K|+|Z_s(H) intersection K|. (1)

This is finite rectangle-union geometry using S,H,L,r. It does not require
W, Wprime, the full remaining halo, or the true deletion region.

Let I_A(S) be the union of all outer receiver rectangles E_T(m(T)) over
eligible source subsets T meeting S. The R67 W and Wprime are contained
in this touching-receiver union. The main containment is

    I_A(S) subset (C\K) union (Z_s(H) intersection K)
                      union Halo(nu).                       (2)

In particular,

    |Wprime_A(S)| <= |W_A(S)|
       <= Gamma_r(H)+9(E+beta).                            (3)

Thus existence of a maximum block S and a nonempty H with right side at
most7L/2 is a general3.5 payment certificate. There is no restriction on
source count or the number of eligible labels.

## 2. The three-label proof of the containment

Fix any eligible T meeting S, write m=m(T). Maximality gives0<m<=L.
For every x in its open outer receiver E_T(m), all q in T obey

    ||x-q||_infinity<sqrt(m)/2.                            (4)

This already gives x in Z_s(S) subset C. Divide into three cases.

### A. T meets the core

Choose p in T intersection H. Equation(4) gives x in Z_s(H). If x is inK,
it belongs to the second set in(2); otherwise it belongs to C\K.

### B. T misses the core and m<=E

T contains a point p in V because it meets S and misses H. For every outside
source q in T\S, equation(4) at the two source points gives

    ||q-p||_infinity<sqrt(m)<=sqrt(E)=delta.

Therefore q belongs to B, and T subset V union B. The FULL diameter is
sqrt(m), not its half-radius; this is why delta is sqrt(E). The same strict
witness capturing T is a witness for nu, whose mass on the square is at
least m. Hence E_T(m) subset Halo(nu).

If E=0 this case is impossible for a positive-mass T, so the empty-V/B
convention causes no exception.

### C. T misses the core and m>E

The source sets H,T are disjoint. Consequently

    M+m>L.                                               (5)

Suppose x in E_T(m) intersection K. We prove that H union T would have
bbox span strictly less than sqrt(M+m), hence would be eligible with mass
greater than L, contradicting the choice of S.

There are three possible types of endpoint pairs controlling a coordinate
span of H union T:

1. H-to-H spans are at most d(S)<sqrt(L)<sqrt(M+m), because H subset S and
   S itself is strictly eligible. This argument does NOT require H eligible.
2. T-to-T spans are less than sqrt(m)<sqrt(M+m), by eligibility of T or(4).
3. For p in H,q in T, closed membership x in K and the STRICT bound(4) give

       ||p-q||_infinity < r+sqrt(m)/2 <= sqrt(M+m).          (6)

   The last inequality follows from the already established scalar minimum

       sqrt(M+m)-sqrt(m)/2 >= sqrt(3M)/2 >= r.

   It is also the elementary inequality
   M+3m>=2sqrt(3Mm). It holds for every M,m>0. The strict slack in(6) comes
   from the open receiver, so the boundary of closed K is also excluded.

All coordinate endpoint pairs are strictly below the required span.
The eligible union H union T has mass M+m>L, the desired contradiction.
Thus high-mass labels missing H have NO outer receiver in K.

These three cases prove(2). Bbox closure of a label is not needed anywhere;
the argument holds for every eligible T, and therefore for the closed
subfamily defining the R67 W/Wprime unions. No residual-core area or full
remaining halo is subtracted.

## 3. A self-contained weak9 bound for the small input nu

For any finite positive atomic measure mu of total massN,

    |Halo(mu)|<=9N.                                       (7)

If N=0 the assertion is immediate. Otherwise take a compact set F inside
its strict halo. For each x in F choose a CLOSED centered witness square
Q_x of positive half-side with mu(Q_x)>|Q_x|. Its interior contains x,
so compactness supplies a finite subcollection whose interiors cover F.

Order this finite collection by nonincreasing radius. Retain a square only
if it is disjoint, INCLUDING boundaries, from all previously retained
squares. Every discarded Q intersects a retained square R with radius at
least its own; therefore Q subset 3R by the L-infinity triangle inequality.
The retained closed squares are genuinely pairwise disjoint, so no atomic
boundary mass is charged twice. Consequently

    |F| <= sum_retained |3R|=9 sum_retained |R|
         <9 sum_retained mu(R)<=9N.

The finite-input strict halo is measurable (indeed open and bounded);
compact exhaustion proves(7). Applying it to nu, with N=E+beta, proves(3).
This is a single source-mass charge for a finite subinput. It is not a
per-label sum of the same a_T masses.

## 4. Two exact refinements and a finite choice, without a radius scan

Case B also gives x in Z_delta(V):=union_(p in V)(p+(-delta/2,delta/2)^2).
Thus the low branch is in BOTH Halo(nu) and Z_delta(V), and the last term
of(3) may be replaced by

    min{9(E+beta), |Z_delta(V)|},

or the simpler min{9(E+beta), |V|*E}. This can be useful when V is small;
the mass bound remains useful for arbitrarily many tiny positive sources.
For E=0 take this low union empty. This refinement does not create any
new label family or change the near-source definition.

For fixed H, K_H(r) grows as r increases, and

    Gamma_r(H)=|C|-|(C\Z_s(H)) intersection K_H(r)|.

Therefore Gamma_r is nonincreasing in r. The maximal permitted radius
r_star=sqrt(3M)/2 is always the best such geometric bound; the low branch
cost is independent of r. Rational smaller radii may be used for simple
exact certificates by verifying4r^2<=3M. No radius parameter scan is needed.
The resulting sufficient test is a finite choice of nonempty H subset S,
not a restatement involving W or the complete remaining halo.

## 5. An explicit arbitrary-rank dominant-corner class

A simple uniform class follows without computing an arbitrary core union.
Assume w,h>3s/4 and there is a source p at a bbox corner of S. Let H={p},
M=m(p), and define E,beta as above. If

    E+beta<=L/32,                                        (8)

then

    |Wprime|<=|W| < (2767/800)*L < (7/2)*L.               (9)

Proof. Here M>=31L/32. Choose r=17s/20, which is permitted since

    3M-4r^2 >= (93/32-289/100)L=(13/800)L>0.

After coordinate reflections put p at the upper-right corner. Fat widths
ensure each side of C intersection K has length r+s/2=27s/20. The single
square Z_s(H) has areaL and is inside K. Thus

    Gamma_r(H)=|C|-(729/400)L+L=|C|-(329/400)L.

Since |C|=(w+s)(h+s)<4L, combining(3) and(8) gives

    |W| < [4-329/400+9/32]L=(2767/800)L.

This class permits arbitrary rank and arbitrarily many tiny positive
filling sources. If no outside source lies within sqrt(E) of V, beta=0,
so E<=L/32 suffices. Condition(8) can also tolerate a positive nearby
outside mass. It is not a support-only tree condition.

## 6. Interpretation and scope

The method makes every eligible label use its OWN total captured mass to
choose one of three analytic branches, but never adds the corresponding
label masses. Meeting the core is paid by a global equal-radius core union;
low total mass is paid by one small-mass input or one small-radius support
union; high total mass missing the core is excluded from a protected
intersection by maximum-block geometry. These jointly differ from a
fixed-theta filled bbox or an equal-radius support-neighborhood repair.

Gamma must be implemented with C intersection K. Even a valid H can have
K extending outside C; replacing |C intersection K| by |K| is wrong.
The square-root near test is exact on rational data through
(dist_infinity(q,V))^2<E, with strict inequality.

The theorem pays a new general structural class, not every maximum block.
No claim of a general heavy/diffuse dichotomy, unrestricted Wprime3.5
selection, improved general C2 bound or matching extremizing sequence is
made here. Controlled certificates for existing fixed inputs and their
actual Python/SymPy checks are archived separately. Frozen R68 files and
other workspaces remain untouched; no LP, MATLAB or extra Agent is used.
