# R62 coupled whole-cell source residual and adaptive exact certification

2026-11-15. This is a constructive continuation of the positive occupation
interface, not a computed improvement of the general constant. The owner is
`pcw-sol-kernel-sources`; MATLAB/Python batch execution belongs to root.

## Objects and inherited scope

Keep the R60 normalization: rho(x)=|x1|+|x2|, U=1/rho,
A=−|D1|−|D2|, AU=−2*pi*delta0 and
Qh=U−PhU=2*pi*integral_0^h p_s ds≥0. Here

    k_h(t)=h/[pi*(h²+t²)],  p_h(x)=k_h(x1) k_h(x2).

The source support is a finite union of axis-parallel square cells. We use
D4 orbits only to reduce unknowns; each coefficient a_j≥0 is a density per
physical source cell in orbit j. Let

    phi_j(x)=sum_{C in orbit j} integral_C p_h(x−y)dy,
    g(x)=lambda*p_h(x)+sum_j a_j*phi_j(x),  lambda≥0,
    R_i(x)=g(x)−a_i  on physical target cell C_i.

The target value a_i belongs to the orbit of C_i. The mathematical condition
is R_i≥0 throughout each cell (a.e. suffices). Outside their union f=0, hence
g−f≥0 automatically. By R60, certified source positivity gives

    K=Qh*(lambda*delta0+f)≥0,
    AK=2*pi*[(lambda*p_h+Phf−f)−lambda*delta0],
    (1/2)*integral K=pi*h*(lambda+sum_j a_j*|orbit j|*delta²).

Source positivity alone does not give a new weak-type constant: a certified
lower bound m>0 for K on the prescribed diamond D is still necessary.
The eventual general bound is the displayed halfmass divided by m. A sampled
D minimum or ordinary Gaussian quadrature does not certify m.

## 1. A coupled rectangle residual lemma

Let B=[l1,r1]×[l2,r2] lie inside a target cell. Put ell_k=r_k−l_k.
Suppose g is C² on B and

    partial_11 g≤Mx,  partial_22 g≤My  throughout B,
    Mx,My≥0.

Then

    inf_B R_i ≥ min_{v vertex of B} R_i(v)
                −(Mx*ell1²+My*ell2²)/8.                 (1)

Only the full coupled values g(v) occur in (1). We never replace the sum of
positive basis functions by the sum of their separate infima. There is no
mixed-derivative term.

Proof. Let I1,I2 be endpoint linear interpolation in coordinates 1 and 2.
For a one-variable C² function F on an interval of length ell, the positive
Green-kernel representation gives

    I F(t)−F(t)=integral G(t,s)*F''(s) ds,
    G≥0, integral G=(t−l)*(r−t)/2≤ell²/8.

Consequently IF−F≤M*ell²/8 when F''≤M and M≥0. Since I2 has nonnegative
weights,

    I2 I1 g−g=I2(I1 g−g)+(I2 g−g)
             ≤Mx*ell1²/8+My*ell2²/8.

Bilinear interpolation I2 I1 g is a convex combination of the four vertex
values. Subtract the constant a_i to obtain (1). Endpoints shared by target
cells cause no difficulty: check each closed target cell with its own
coefficient; this is stronger than the required a.e. condition.

## 2. Linear proposal rows without entrywise-inf loss

For the basis phi_0=p_h, obtain upper bounds U_j,11 and U_j,22 for their
second derivatives on B. The bounds may be signed. For each basis define

    E_j(B)=[ell1²*max(0,U_j,11)+ell2²*max(0,U_j,22)]/8.

Nonnegative coefficients imply the sufficient linear rows, for every vertex v,

    lambda*[phi_0(v)−E_0(B)]
      +sum_j a_j*[phi_j(v)−E_j(B)]−a_i ≥ 0.           (2)

Indeed partial_kk g is bounded above by the corresponding coefficient
combination of max(0,U_j,kk), which can be used in (1). Signed derivative
intervals are first summed over all physical cells of a D4 basis, before the
max operation. This preserves internal orbit cancellations. Formula (2)
still need not exploit cancellation between distinct orbit coefficients;
that is recovered by the postchecker in Section 3.

A slightly stronger LP can introduce two nonnegative auxiliaries per B,
eta_k≥lambda*U_0,kk+sum_j a_j*U_j,kk, then replace the E sum by
(ell1²*eta_1+ell2²*eta_2)/8. This exactly retains the coefficient combination
of these signed upper bounds. The present bounded pilot deliberately avoids
those many extra columns. Both formulations have O(ell²) curvature loss at
fixed h; the old entrywise-inf construction generally has O(ell) loss.

`computation/cauchy_occupation_coupled_lp.m` implements (2) in ordinary
MATLAB doubles, with optional multicorner proposal mode that omits E_j.
It always labels outputs PROPOSAL. Its source rows, LP optimum, Gaussian
quadrature and D samples do not constitute a directed certificate.

### 2.1 Version 2: tight one-dimensional extrema and signed auxiliaries

The first root-owned version-1 pilot at h=.15, delta=.125, subdiv=2 had
floating cost about 6.97777284373. Its natural derivative interval evaluation
and separate positive curvature contributions were too costly. Version 1
is frozen; `cauchy_occupation_coupled_aux_lp.m` is a distinct version 2.

For k'' the derivative is

    k'''(x)=24*h*x*(h²−x²)/[pi*(h²+x²)^4].

Therefore the tight range on a real interval uses its endpoints and the
critical points 0,±h that lie in that interval. For k', use endpoints and
±h/sqrt(3) when present. Their exact values are respectively

    k'(h/sqrt(3))=−9/[8*sqrt(3)*pi*h²],
    k'(−h/sqrt(3))=+9/[8*sqrt(3)*pi*h²].

The v2 exact checker detects critical-point membership by rational squared
comparisons with h²/3, including the correct positive/negative half-line;
no floating root location is used. Sqrt(3) itself is bracketed by
floor(sqrt(3*2^320))/2^160 and the next lattice point, using integer isqrt.
Endpoint derivative values are rational multiples of 1/pi. Hence these
tight ranges remain genuinely directed in the exact checker. The MATLAB
proposal uses the corresponding ordinary floating formulas.

Version 2 implements the two signed-curvature auxiliaries from Section 2:
eta_x≥sum_j v_j U_j,11, eta_y≥sum_j v_j U_j,22, eta≥0. Four source rows per
subrectangle require

    sum_j v_j phi_j(vertex)−a_i
       ≥ell²*(eta_x+eta_y)/8.

There are B=J*subdiv² subrectangles, J+1+2B variables and 6B source/curvature
rows. Eta columns are sparse. This preserves cancellation across distinct
orbit coefficients in the LP itself. Upper derivative intervals for each
orbit still safely enclose the sum over its physical cells; tighter
one-dimensional ranges reduce this remaining interval loss.

Version-2 code is separate from the frozen version-1 actual runs:
`r062_coupled_source_verify_v2.py` has the same CLI and coverage/repair
semantics, replacing only the derivative range oracle. Neither version may
claim a new general constant without a source PASS and full-D certificate.

## 3. Actual-coefficient interval verification

Once rational lambda,a are fixed, evaluate derivative intervals for the
**actual complete g**. If I_j,kk encloses partial_kk phi_j(B), interval
arithmetic gives

    J_kk=lambda*I_0,kk+sum_j a_j*I_j,kk,
    Mk=max(0,upper J_kk).

These quantities may be substantially smaller than the LP proposal's sum of
positive upper bounds. Compute a directed interval for every vertex R_i(v),
put L_B=min_v lower R_i(v)−(Mx*ell1²+My*ell2²)/8, rounding the error upward.
If L_B≥0, the complete rectangle is certified by (1). Otherwise bisect the
longer side and repeat. Every complete initial cell partition and every
subsequent bisection is recorded; a PASS requires no pending or unresolved
rectangle. A budget expiration is INCOMPLETE, never PASS.

The one-dimensional derivatives used are

    k'_h(t)=−2*h*t/[pi*(h²+t²)²],
    k''_h(t)=2*h*(3*t²−h²)/[pi*(h²+t²)³],
    T_C(x)=integral_[c−delta/2,c+delta/2] k_h(x−y)dy,
    T'_C(x)=k_h(x−c+delta/2)−k_h(x−c−delta/2),
    T''_C(x)=k'_h(x−c+delta/2)−k'_h(x−c−delta/2).

Thus for a physical rectangle phi(x)=T_1(x1)T_2(x2),
phi_11=T''_1*T_2 and phi_22=T_1*T''_2. Probability T_C is even and decreases
with |x−c|: its derivative is negative for x>c because k is even and strictly
decreases on the positive half-line. Its interval range on B therefore uses
the minimum and maximum absolute distance exactly. The derivatives use
rational interval evaluation of the displayed formulas; their natural
interval dependency widens a bound but never invalidates it.

The standalone verifier `computation/r062_coupled_source_verify.py` uses
only Python integers and Fraction input parsing. Every number is enclosed
on the common lattice 2^−160. Signed products use the minimum and maximum
of all four endpoint products, with outward integer floor/ceil; reciprocals
are used only for strictly positive denominator intervals. Rational input
coefficients are scaled with outward integer rounding.

For arctan, arguments are reduced to [0,1/2] by oddness, reciprocal and
pi/4 identities. The 88-term alternating series is enclosed together with
the next-term absolute remainder. Pi is enclosed using
16*atan(1/5)−4*atan(1/239); the tangent identity is also checked exactly.
These are genuine directed bounds, not high-precision floating estimates.
The inherited R60 checker used a related construction; this implementation
adds signed interval multiplication, derivative enclosures and the explicit
whole-support partition proof.

A strict negative upper interval at a closed-cell corner gives a real
invalidity witness: continuity of g−a_i implies negativity in a one-sided
positive-area neighborhood inside that target cell. This is a rejection of
that candidate's source condition, not a general lower bound on the whole
occupation family. Such a witness does not require full support coverage.

## 4. Certified repair and its exact mass cost

Suppose every target cell is covered by a finite terminal partition and
R_i≥L_B there. Let p_B>0 be any certified lower bound for p_h on B. Then

    Delta=max_B max(0,−L_B/p_B)                       (3)

is a sufficient additional origin amplitude. Indeed R_new=R_old+Delta*p_h
≥L_B+Delta*p_B≥0 on each B. Formula (3) uses the old residual bounds;
there is no need to recompute their curvature after adding Delta.
The exact halfmass penalty is pi*h*Delta. Since Qh≥0, a previously certified
D minimum does not decrease. If the original D minimum was sampled only,
repairing the source does not fix that separate gap.

There is also a termination statement. If a fixed finite candidate has
R_i≥epsilon>0 on every closed target cell, then the exact adaptive test
terminates after finitely many longest-side bisections when its arithmetic
precision and the arctan truncation degree are chosen sufficiently high together. All finitely many derivative basis
functions are continuous on the compact support; their interval enclosures
remain bounded, the interpolation error tends to zero with mesh squared,
and vertex evaluation widths tend to zero when BOTH the fixed-point precision
and the arctan series degree increase. Increasing BITS alone with a fixed
88-term series does not remove its truncation floor. Choose precision, series
degree and mesh so both errors are less than epsilon/2.
If a legal source has zero slack, adding any Delta>0 produces uniform slack
at least Delta*min_support p_h>0. Therefore arbitrarily small explicit origin
repair removes the zero-slack certification obstruction, with arbitrarily
small exact halfmass penalty. This proves a constructive certification
route; it does not give an efficient complexity bound for a particular
candidate. The implementation fixes 160 bits for the bounded pilot and
reports unresolved cases instead of invoking this existence argument as a
computed certificate.

The verifier outputs a sufficient extra repair only when its accepted and
unresolved terminal leaves cover every initial cell. It explicitly outputs
null if the time/box budget leaves any unchecked rectangle. `--repair q`
first adds an exact rational q to lambda, and then verifies that new source.

## 5. Bounded root-owned pilot and decision rule

No old midpoint/cellinf candidate is rerun. Inputs are specified in
`computation/inputs/r062-coupled-source/pilot.json`.
A first useful pilot has h=.15, delta=.125, N=16, subdiv=2 (153 D4 unknowns
plus lambda, 2448 source rows); subdiv=4 has 9792 rows. Each MATLAB process
uses one numerical thread and is owned, time-bounded and explicitly exited
by root. The exact checker has a separate 45-second and 12000-box budget.

Compare robust-row cost, strict source result and the necessary/sufficient
repair penalty. A source-certified cost still above 3.590865312933 does not
improve the inherited general upper bound. A lower proposal is only worth
full-D certification if the directed repair penalty leaves a real margin.
This is a general upper-bound kernel interface; it does not substitute a
narrow-family exclusion or a negative search for the main objective.

## 6. A positive-integral full-D certificate interface

If an actually source-certified candidate reaches a competitive cost, full-D
verification can use the same integer/arctan engine rather than certify the
ordinary Gaussian quadrature. Fix 0<epsilon=s_0<...<s_M=h. For a test box X,
let d_k=max_{x_k in X_k}|x_k|. For a source-cell coordinate centered at c,
let d_C=max_{x_k in X_k}|x_k−c|. Define

    b_k=min(k_{s_m}(d_k),k_{s_{m+1}}(d_k)),
    t_C=min(T_{s_m,C}(d_C),T_{s_{m+1},C}(d_C)).

For each fixed distance, k_s(d) has at most one interior maximum as a
function of s, and its minimum on a closed positive time interval is at an
endpoint. The same is true for T_s(d): direct differentiation gives

    partial_s T_s(d)
      =−delta*[s²−(d²−delta²/4)]
        /[pi*(s²+(d−delta/2)²)*(s²+(d+delta/2)²)].

Both are decreasing in absolute spatial distance. Thus on X×[s_m,s_{m+1}],
p_s≥b_1*b_2 and each physical rectangle probability is at least t_C1*t_C2.
These lower bounds use the same shared time interval; replacing a factor by
its own minimum is valid by positivity, and improves under time refinement.
Consequently

    inf_X K ≥ 2*pi*sum_m (s_{m+1}−s_m)
      *[lambda*b_1*b_2+sum_physical_cells a_C*t_C1*t_C2].     (4)

The omitted interval [0,epsilon] contributes a nonnegative amount and may
be dropped. Rational interval lower endpoints and a lower pi bound make
(4) directed. It covers the singular origin without evaluating U there.
Cover D by boxes, discarding boxes disjoint from D and allowing arbitrarily
small outward overhang at its boundary; D4 permits one fundamental wedge.
A finite box partition on which (4) exceeds a target m certifies K≥m on D.
This is a ready analytic interface, not an implemented or run D certificate
in the present source pilot. It gives a direct route from a competitive
source candidate to a genuine general upper bound.

## Evidence and remaining work

Root reported an actual successful 160-bit exact tiny self-test and one
MATLAB pilot (h=.15, delta=.125, subdiv=2), cost about 6.97777284373 with
floating source-row minimum about −1.95e−14. Root owns the saved run receipts
and process cleanup; this author has not independently rerun those jobs.
The pilot exposes substantial curvature enclosure loss and is not a new
upper bound. The source checker subsequently returned PASS on all 612 subrectangles with
no repair (1.647 seconds in the saved output). The exact halfmass divided by
pi is 355375052752634298393/160000000000000000000. The candidate, full source
certificate and root execution receipt are copied byte-for-byte to
`computation/results/r062-coupled-source/v1-h015-d125-s2/`. No new C2 bound,
full-D certificate, commit or push is claimed here.

Further actual V2/V3 comparisons and the interface-jump proof are recorded
in `r062-jump-telescoped-source-curvature.md`; this main proof retains the
frozen V1 derivation and its original executed certificate.

Source files actually read: R60 independent positive-cell occupation proof,
R60 MATLAB `cauchy_occupation_cell_lp.m`, the kernel R60 handoff and the R60
independent occupation audit. Their SHA256 provenance is recorded in the
handoff manifest, distinct from the new code's execution status.
