# R65 nonuniform tensor rectangles: positive occupation source interface

2026-11-15. This is a new R65 increment. R62's published source, proofs,
certificates and manifests remain frozen (root reports remote-verified commit
39a92c77ed8ad365dcea221b9d716dc8a7d74114). The occupation/C4 theorem and the
coupled whole-cell interpolation lemma are inherited from R60/R62; the new
work is an exact nonuniform tensor geometry, its correct orbit/area/jump
representation, and a separate proposal/checker implementation.

## 1. Objects, geometry and density units

Keep A=−|D1|−|D2|, U=1/(|x1|+|x2|), Q_h=U−P_hU≥0,
p_h(x,y)=k_h(x)k_h(y), k_h(t)=h/[pi(h²+t²)]. Let

    e_0<e_1<...<e_m, e_k=−e_(m−k), m odd,
    I_i=[e_i,e_(i+1)], w_i=e_(i+1)−e_i,
    c_i=(e_i+e_(i+1))/2.

All endpoints are exact rational numbers. The central interval crosses zero;
the same one-dimensional partition is used on both axes. Physical cells are
C_ij=I_i×I_j. Let f=sum a_ij*1_Cij with a_ij≥0 a density, not a cell mass.
D4 maps a rectangle by reflecting its center and by simultaneously swapping
its two center coordinates AND its widths. Thus symmetry requires

    a_ij=a_(m−1−i,j)=a_(i,m−1−j)=a_ji.

Taking one representative for each orbit is legitimate. Reflected widths
agree by exact endpoint symmetry, and swapping coordinates preserves area.
If orbit q has size n_q, representative widths w_qx,w_qy and density a_q,

    integral f = sum_q n_q*a_q*w_qx*w_qy.            (1)

There is no common delta² factor. A coarse/fine rectangle and its transpose
belong to the same orbit but their width coordinates must swap.

Set T_i(x)=integral_Ii k_h(x−z)dz. The orbit probability basis is the sum of
T_i(x)T_j(y) over every physical rectangle in that orbit. Then for lambda≥0,

    g=lambda*p_h+P_h f,
    K=Q_h*(lambda*delta0+f),
    R_ij=g−a_ij on C_ij.

If all complete target cells have R_ij≥0, exterior source positivity is
automatic because f=0 there. R60 gives K≥0, K in L1, the legal source
lambda*p_h+P_hf−f≥0, and

    halfmass(K)=pi*h*[lambda+sum_q n_q*a_q*w_qx*w_qy]. (2)

This is a general positive-source kernel construction. A certified lower
bound for K on the entire diamond D remains necessary before converting
(2) into a general C2 bound.

## 2. Tensor-interface jumps retain exact cancellation

For a vertical interface x=e_k in source row j, define

    d^x_(k,j)=a_(k,j)−a_(k−1,j),  k=0,...,m,

with exterior densities zero. Since

    T_i''(x)=k_h'(x−e_i)−k_h'(x−e_(i+1)),

finite summation gives the exact continuous identity

    (P_h f)_11(x,y)
      =sum_(k,j) d^x_(k,j)*k_h'(x−e_k)*T_j(y).     (3)

The probability T_j in (3) has the ACTUAL width w_j of that source row.
There is no division by either adjacent x width. Different x widths do not
break the cancellation: adjacent intervals share the same endpoint, and
their density values multiply identical interface functions. Likewise

    (P_h f)_22(x,y)
      =sum_(i,k) [a_(i,k)−a_(i,k−1)]*T_i(x)*k_h'(y−e_k). (4)

Here T_i has the actual column width w_i. A constant-density tensor block
has zero internal jumps despite nonuniform widths. The outer zero-density
terms in (3)–(4) must be retained.

The tensor condition matters: every rectangle in one source row shares
that row's y interval. An arbitrary quadtree with unmatched row intervals
would require splitting interface segments first; the present code does
not silently apply (3) to such a partition.

## 3. Coupled anisotropic whole-cell constraints

For a target subrectangle B with side lengths ell_x and ell_y, suppose
eta_x,eta_y≥0 bound g_11,g_22 from above. R62's one-sided interpolation
lemma gives

    inf_B(g−a_target) >= min_four_vertices(g−a_target)
                          −ell_x²*eta_x/8−ell_y²*eta_y/8. (5)

The two errors generally differ; sharing one delta² coefficient would be
incorrect for coarse/fine rectangles.

For an interface function psi in interval[L,U], split its signed jump as
d=dplus−dminus, dplus,dminus≥0. Then

    d*psi <= dplus*U−dminus*L.                     (6)

As in R62, complementarity is unnecessary. Let ox,oy be signed upper
bounds for the origin p_h curvatures. Linear constraints

    eta_x >= lambda*ox+sum_vertical(dplus*U−dminus*L),
    eta_y >= lambda*oy+sum_horizontal(dplus*U−dminus*L),
    eta_x,eta_y>=0

together with four vertex rows using the COMPLETE original positive
convolution and separate errors from(5), are sufficient whole-cell source
constraints. All coefficient units remain densities.

## 4. Exact pilot geometry and implementation mapping

The prepared axis has fine intervals of width1/32 on[-17/64,17/64], including
the central cell[-1/64,1/64]. Fourteen coarse intervals of width1/8 extend
each side to outer boundary129/64. There are45 axis intervals,2025 physical
rectangles and276 D4 density orbits. The positive-edge list is

    1/64,3/64,...,17/64,
    17/64+1/8,...,17/64+14/8=129/64,

and its reflected negative list completes the partition. Every edge is
listed explicitly in
`computation/inputs/r065-tensor-occupation/central-fine-outer-coarse.json`.
This is the minimum exact geometry input, rather than an inferred grid
from printed floating centers. It is constructed with Python Fraction and
integer comparisons; no LP/kernel experiment was run by this author.

For a common subdivision count s, J=276, B=J*s² target subrectangles,
E=2*m*(m+1)=4140 physical interfaces. The proposal has

    J+1+2E+2B variables,
    6B source/curvature inequality rows,
    E jump equality rows,

plus D sampling rows. With s=2, B1104 and10765 variables. Fine global
resolution1/32 out to radius about2 would require many more source cells;
the present input concentrates fine density representation near the origin.
This is a cost choice, not a proof that the chosen mesh is accurate enough.

`computation/cauchy_occupation_tensor_jump_lp.m` has signature

    (h, grid_input, subdiv, output_root).

It reads exact rational strings to doubles for proposal construction,
groups orbit descriptors(center_x,center_y,width_x,width_y), and builds an
explicit physical-cell orbit map. JD has right-minus-left and top-minus-
bottom rows, including both exterior boundaries. Variable order remains
(base lambda,a), all dplus, all dminus, eta_x, eta_y, with equality
[JD,−I,+I,0]. Curvature edge ranges are vectorized, using each row/column's
actual width. The source matrix has distinct EX and EY entries
ell_x²/8 and ell_y²/8. The mass cost uses orbit_size*actual_cell_area.

The source output embeds `axis_edges_exact` from the input. Center, width
and area metadata are exact dyadic values for this pilot. An arbitrary
nondyadic rational partition can still define the mathematical theorem;
this prototype's optional exact metadata consistency assertions require
exact string metadata or omission of those floating diagnostic fields
if such a partition is later introduced. The current prepared dyadic input
has no such ambiguity.

D samples and Q_h*rectangle integrals use ordinary Gauss24, now with the
separate widths of each source rectangle. This remains a numerical
proposal; no full-D directed bound is provided by the MATLAB code.

## 5. Separate exact whole-support source checker

`computation/r065_tensor_source_verify.py` parses embedded edge strings with
Fraction, asserts exact symmetry/strict ordering/central crossing, and
independently rebuilds the D4 orbits and physical cell map. It recomputes
all areas, every actual signed density jump and row/column probabilities.
It uses original coupled vertex values, never LP jump parts or eta values.
Its directional interpolation error is computed separately for the actual
rational x/y side lengths.

Arithmetic is inherited from the reviewed R62 V3 checker:160-bit outward
signed integer intervals,88-term arctan series with exact alternating
remainder, Machin pi bounds and exact detection of k' critical points via
squared rational comparisons/integer sqrt3 bracket. Source coverage is an
explicit initial partition followed by recorded longest-side bisections.
PASS requires every source cell covered, no pending/unresolved rectangles.
A strict negative upper vertex value gives a one-sided interior witness;
a budget exhaustion reports INCOMPLETE. Complete terminal coverage can
supply a sufficient extra lambda repair, with exact pi*h*Delta cost.

As in corrected R62, arbitrary strict-source termination is an existence
statement that jointly increases precision AND series degree and refines
boxes. Fixed160/88 does not guarantee all small-slack candidates terminate.

The prepared tiny self-test uses a genuinely nonuniform three-interval axis
[-2,−1/2,1/2,2], h=4 and constant density1/10000. It checks actual total
area16, central/outer widths1 and3/2, transposed rectangle widths, exact
mass1+16/10000, disappearance of every interior constant-density jump,
and whole-cell source positivity. Root must actually execute it before
recording a computational PASS. Author validation so far is Python AST
parsing and exact input geometry construction only.

## 6. Bounded root-owned first test and current scope

Root chooses the first h=.075 or.15, with the SINGLE prepared grid and s=2;
no extra mesh family or automatic batch is requested. MATLAB should be
launched with requested thread controls and a bounded lifecycle, then exit
promptly; actual receipt and process absence belong to root. The exact
checker has an initial60-second/16000-box budget, independently adjustable
by root without turning a timeout into a certificate.

At preparation no optimized R65 candidate, source certificate, full-D
certificate or improved general bound exists. First compare saved source
cost and strict source status with R62, then decide whether another test
is worthwhile. The code and input are the minimum new interface for that
decision; no new Agent, MATLAB process, commit or push was started here.
