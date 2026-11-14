# Finite-source receiver dual: sparse completeness and precise obstructions

Date: 2026-11-15. Research note; no new numerical upper bound for C₂.

## 1. Complete finite-source dual and an at-most-N receiver theorem

Fix distinct p₁,…,p_N in R² and positive m_j. Set μ=Σm_jδ_pj, λ=1,
Q(x,r)=x+[-r,r]², and H={x: sup_r μ(Q(x,r))/(4r²)>1}.
All incidences use closed squares; all witness inequalities are strict.
For each nonempty S⊆{1,…,N}, let

    W_S = {(x,r): r>0, {j:p_j∈Q(x,r)}=S, 4r²<Σ_{j∈S}m_j}.
    c_S = sup{4r²:(x,r)∈W_S}.

Discard empty W_S. Every remaining c_S is positive and ≤m(S).
There are at most 2^N−1 labels, irrespective of the continuum of centers/radii.
Define the finite LP

    P(μ)=max Σ_S c_S z_S,
          z_S≥0,  Σ_{S∋j}z_S≤1 (j=1,…,N).                 (1)

**Theorem.** P(μ) is exactly the supremum of B_ν=∫4r²dν over all
nonnegative receiver measures on the labeled space of valid (x,r), whose
loads are ≤1 at every p_j. For every 0<ε<1 there is such a certificate
supported on at most N actual strict witnesses, with

    (1−ε)P(μ) ≤ B_ν ≤ P(μ) ≤ μ(R²).                      (2)

Proof. Push any ν to the finite incidence labels: z_S=ν(W_S). Each
W_S contains some source index, so z_S≤1; the total ν mass is ≤N.
Loads give exactly the constraints in (1), and B_ν≤Σc_Sz_S.
Conversely, (1) is a compact polytope (each variable is ≤1), hence has
an optimal extreme point. At most N coordinates of an extreme point
are positive: otherwise their N-dimensional incidence columns are
linearly dependent and a sufficiently small two-sided perturbation
preserves loads and nonnegativity. For each active label choose an
actual witness with 4r²≥(1−ε)c_S, possible by the definition of a
positive supremum. Put its z_S mass on that labeled witness. Its source
loads equal those of z and its B is at least (1−ε)P. Finally original
masses are primal feasible since m(S)≥c_S, or use Tonelli. □

Finite LP strong duality additionally gives

    P(μ)=min{Σ_j w_j: w_j≥0, Σ_{j∈S}w_j≥c_S for every valid S}. (3)

This theorem resolves the infinite receiver-space issue for fixed finite
support. It does not say that the suprema c_S are attained, or that
computing them and all valid incidences is cheap. The finite source LP
from the legacy source underlies (1); aggregation over *all* valid
witnesses and sparse ε-realization are the additions here.

## 2. Exact missing all-N geometric lemma

For any coordinate bounding box B of source points, including degenerate
boxes, put M_B=μ(B), s_B=√M_B/2 and

    E_B=(max_x(B)−s_B,min_x(B)+s_B)
        ×(max_y(B)−s_B,min_y(B)+s_B),

discarding empty rectangles. The finite-source version of the complete
halo proof gives H=∪_B E_B without any scale truncation. In the forward
direction take the bounding box of all atoms captured by a strict square;
in the reverse direction any radius between the maximum distance and
s_B is a strict witness.

The sharpness-capable target is the following geometric statement, with
one uniform K independent of N, source locations, positive masses, and
coordinate-box degeneracies:

    |∪_B E_B| ≤ K P(μ).                                  (G_K)

If proved, (2) implies |H|≤Kμ(R²), with arbitrarily small relative
certificate slack and at most N receiver atoms. No numerical mesh,
bounded source number, or radius menu is built into (G_K). For an actual
receiver certificate the ε form is B_ν≥(1−ε)|H|/K; a strict endpoint
B_ν≥|H|/K is not promised because c_S need not be attained.

This is not an evaluated formula for C₂. Its value is that all missing
work has become a definite geometry-versus-finite-LP inequality. It is
sharp in strength: if K is any valid global finite-measure weak constant,
take a minimizer w in (3). Each original strict witness with incidence S
has area ≤c_S≤w(S). For any η>0, the input (1+η)Σw_jδ_pj thus makes
every original halo center strict; hence |H|≤K(1+η)P(μ). Let η↓0.
Thus the ordinary sharp weak inequality implies (G_K) as well. This
equivalence is a reduction, not a resolution of its constant.

**Route decision.** Source-dependent load checks are sufficient for a
general theorem if (G_K), or its direct certificate version, is proved
uniformly for every μ. Requiring all-y loads is an optional stronger
strategy, and the next obstruction explains why that stronger strategy
may have a substantially smaller feasible objective on a chosen family.

## 3. Unbounded gap between support and all-space loads on a fixed family

For arbitrary N≥1 let t_i=i/(N+1), x_i=(t_i,1−t_i), r_i=1,
p_i=(1+t_i,2−t_i), and μ=5Σδ_pi. Every Q_i contains the origin.
Also p_j∈Q_i implies t_j≤t_i and t_j≥t_i, hence j=i. Therefore all
Q_i are strict witnesses, μ(Q_i)=5>4. On this selected family the source
load LP is max 4Σz_i with 0≤z_i≤1, of value 4N. The all-space load LP
additionally has Σz_i≤1 at the origin, and its value is exactly 4.

This is an arbitrarily large gap between two finite certificate problems.
It is **not** a gap between weak constants, and the finite set {x_i} has
zero Lebesgue area. In particular it does not refute μ-dependent
support-only proofs. It refutes transferring a fixed support certificate
unchanged into an all-y certificate. The all-y obstruction occurs at an
empty source location and persists even with equal radii.

## 4. Two precisely limited simplification failures

**Singleton halo reconstruction.** For μ=δ_{(-1/2,0)}+δ_{(1/2,0)}, the
two singleton halo squares have area 2 in union. The two-source box adds
the rectangle of width √2−1 and height √2. Its part outside the singleton
union has area (√2−1)²=3−2√2>0. Thus the complete halo has area
5−2√2. Omitting nonsingleton source boxes does not reconstruct even this
two-source halo. This does not prove that all certificates supported on
singleton witnesses fail to prove any proposed K; that stronger statement
would require a separate argument.

**Unweighted intersection-nerve rules.** Specify the class: a rule sees
only which subfamilies of labeled witness squares have nonempty common
intersection, is permutation-equivariant, and assigns receiver weights
without seeing areas, radii or source masses. On a complete intersection
nerve it must give every square the same weight, at most 1/N. Take N
concentric squares, one with radius 1, the others radius 1/N (duplicates
are valid labeled witnesses), and source μ=5δ_0. All are strict witnesses.
The true receiver LP value is 4, while the rule has objective at most

    4(1+(N−1)/N²)/N.

Its approximation factor deteriorates at least as
N/(1+(N−1)/N²), without bound. The same example defeats the particular
local rule z_i=1/max_{y∈Q_i}#{j:y∈Q_j}. The example does not rule out
all conceivable rules using overlap counts together with geometry,
areas, or source masses. The term “local states” without a specified
information model is too broad to admit a general impossibility claim.

## 5. Does P(μ) actually improve on total mass?

No such improvement has been established here. The second program
`computation/general_dual_mass_probe.py` computes the c_S using continuous
radius incidence, rather than replacing them by m(S). All changes in the
one-dimensional order of the center-boundary lines p_j,k±r occur at
r=|p_j,k−p_l,k|/2. Between successive such critical radii a<b, the family
of realizable exact source incidences is constant. It is enumerated at
(a+b)/2 using all coordinate endpoints and intervening intervals. For
each incidence S in this open radius layer, valid witnesses exist iff
4a²<m(S), and this layer contributes

    min(m(S),4b²)

to c_S. Critical radii themselves are separately enumerated, including
zero-area center cells; a contribution 4r² requires 4r²<m(S). A final
cap beyond all arrangement changes and beyond √μ(total)/2 terminates
the search. These statements follow directly from the constant strict
ordering between critical radii and the exact-capture inequalities.

Actual run: seed 20261115, 100 five-atom inputs with distinct points
sampled from {0,…,4}² and positive masses in {1/4,…,39/4}. Python
3.13.18 with installed SciPy1.18.1 was used. HiGHS located dual solutions;
every solution was reconstructed as rational numbers and checked for
nonnegativity, source loads ≤1, and objective **exactly equal** to the
original total mass. All 100 cases passed. Since original masses are
primal feasible, these are 100 exact P(μ)=μ(total) certificates, conditional
only on the explicitly proved coefficient enumeration above. Saved run:
`computation/results/general-dual/mass-probe.json`. Reproduction:

    /Users/zzh/.venvs/main/bin/python3 computation/general_dual_mass_probe.py

**Certificate-archiving rerun, 2026-11-15 00:43:25 Asia/Shanghai.** The
earlier saved output contained only the 100-case summary. After review,
the same seed and all 100 inputs were actually rerun once to save every
case's exact points, positive masses, positive dual support S, exact c_S,
exact z_S, source loads, and objective. The replacement `mass-probe.json`
is 113,533 bytes and contains all 100 certificates; it is not an expansion
of the old summary without recomputation. Each coefficient can be
recomputed from its saved input using the complete radius-arrangement
algorithm above. No floating LP matrices are stored.

At 00:43:38 a separate read-only verification run using Python 3.15.0
and only standard-library Fraction arithmetic recomputed all capacities
and checked all 100 saved certificates. It used no LP solver. Its output
is `computation/results/general-dual/mass-probe-verification.json`.
This is an arithmetic reread using the same coefficient algorithm, not
an independently implemented geometric enumeration. Reproduction:

    python3 computation/general_dual_mass_probe.py --verify computation/results/general-dual/mass-probe.json

Here **P(μ)=μ(total) concerns the LP with supremum capacities c_S**.
It does not supply an actual strict receiver measure with B_ν=μ(total).
Equation (2) realizes arbitrary relative ε approximation. Indeed any
nonzero strict receiver certificate has positive integrated witness
slack, so its B_ν is strictly below μ(total).

A separate targeted numerical solve with central mass 8 at (0,0), masses
1/8 at the other eight points of {−1,0,1}², returned P=9=total mass with
22 valid incidence labels. This targeted result was not subjected to the
100-case rational dual check and is only a numerical observation.

**Decision:** no example P<μ(total) was found, and no theorem P=μ(total)
has been proved. In particular (G_K) currently offers no demonstrated
smaller input quantity and no all-N control beyond a reformulation of
the original sharp inequality. The finite-source sparsity theorem is an
interface result; it should not be presented as a new C₂ upper-bound
method with established quantitative gain. The concrete failures in §4
identify precisely which information cannot simply be omitted.

## 6. Actual exact computation and source record

`computation/general_dual_obstructions.py` uses Python standard-library
Fraction arithmetic, exhaustive finite LP vertex enumeration, and exact
coordinate cell enumeration including endpoints to obtain every global
incidence mask. Proper-subset masks are removed only because receiver
weights are nonnegative. Python 3.15.0 was actually run, without random
input or floating-point tolerances. N=5 gave:

* support LP=20, all-space LP=4;
* nested-square LP=4, equal-weight objective=116/125, ratio=125/29;
* exact quadratic-field extra singleton-missed area=3−2√2>0.

Output: `computation/results/general-dual/exact-checks.json`.
Run from topic directory: `python3 computation/general_dual_obstructions.py`.
No MATLAB session was needed. SymPy was unavailable in default Python;
the final program needs no third-party library. LP checks certify these
finite examples, not (G_K) for arbitrary input.

Source dependencies read in this work:

* Current `proofs/periodic-complete-halo.md`, §2, SHA256
  `7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e`.
  The finite nonperiodic halo formula above adapts its bounding-box proof.
* Read-only legacy `/Volumes/data/CodexProjects/T/Math-Research-planar-centered-square-weak11/planar-centered-square-weak11/proofs/attainment.md`,
  §§5–6 (lines215–326), SHA256
  `94fd2370609f8d3ec321ac7a89158c8e28dc4207f13bc0fe69cede82f0b11348`.
  Legacy claims finite incidence strong duality and Tonelli/slack identity;
  those sections were read and the finite interface was checked here.
  This note does not recertify unrelated attainment claims.

No candidate numerical value of K has been proved or tested here. In
particular K=2 is not assumed. All-N (G_K), and hence the exact value of
C₂, remains unresolved. No commit or push was made by this delegated task.
