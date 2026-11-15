# R65 V3: canonical orbit-density-difference groups

2026-11-15. New independent V3 increment only. The original nonuniform
proposal V1, diagnostic V2, exact grid and source checker remain frozen.
This increment preserves the exact continuous kernel family and groups
identical density-difference linear forms before interval curvature bounds.
It does not assert equality of the old and new floating sufficient LPs.

## 1. Why the previous failure does not imply kernel infeasibility

The original V1 run at h=.075,subdiv2 saved only an assertion after linprog;
its exitflag is unknown. The actual diagnostic V2 rerun returned exitflag0,
algorithm dual-simplex-highs,iterations20000 and a MaxIterations message,
with no returned vector or objective. Thus the V2 failure is an actual
iteration-cap event. It cannot retrospectively identify V1's flag.

Root's separate lambda-only diagnostic returned positive minimum unit
source-row residual3.42737969744517e−5 over4416rows, sampled D minimum
qmin=.0013946309069233043 and finite normalization lambda717.036310493.
This is strong floating finite-formulation feasibility evidence, not an
exact finite-LP or full-D certificate. Mathematical lambda*Q_h feasibility
is independently explained by the positive time-integral argument in
`r065-solver-diagnostic-scope.md`.

## 2. Exact finite density-difference identity

Use the nonuniform tensor geometry and D4 density orbit IDs1,...,J. Give the
exterior density the special ID0 with a_0=0. A physical interface r has
oriented neighboring orbit IDs q_plus and q_minus: right/left for vertical
interfaces, top/bottom for horizontal interfaces. Its jump is

    d_r=a_(q_plus)−a_(q_minus).

If these IDs are equal, d_r is identically zero and the interface may be
omitted from the curvature sum. Otherwise define its canonical unordered
pair (lo,hi),lo<hi, and canonical jump

    d_g=a_hi−a_lo,
    s_r=+1 if q_plus=hi, −1 otherwise.

Then d_r=s_r*d_g exactly, including exterior pairs(0,hi). No spatial-
distance or cell-width identification is made. Interfaces share a group
only because their D4 density-difference linear forms agree up to sign.
Different physical widths and positions remain inside their functions.

Let JD be the original physical jump matrix, GD the canonical group
matrix, and SG the sparse physical-to-group sign matrix. Each nonzero
physical row of SG has one entry+1 or−1; zero physical rows are zero.
Then

    JD = SG*GD.                                    (1)

All entries are integers0,±1. The new MATLAB code explicitly asserts the
zero difference of these matrices. The author also ran an exact Python
integer/Fraction geometry ledger that verifies(1) row by row for the
prepared input, without running an LP or evaluating source residuals.

## 3. Group functions before interval enclosure

For each vertical physical interface retain its original function
psi_r(x,y)=k'_h(x−edge_r)*T_row_r(y), with the ACTUAL row width. The grouped
function for coordinate1 is

    Psi_g,11=sum_(vertical r in groupg) s_r*psi_r.

Similarly define Psi_g,22 from horizontal interfaces and the actual column
width. The exact derivative identities become

    (P_h f)_11=sum_g d_g*Psi_g,11,
    (P_h f)_22=sum_g d_g*Psi_g,22.                   (2)

Thus groups can be shared between the two coordinate directions while
keeping two different curvature functions. Exterior zero terms retain
their proper signs. A group with no contribution in one direction simply
has zero function in that direction.

On a target rectangle, suppose each physical function is enclosed by
[L_r,U_r]. An enclosing group interval is

    L_g=sum_(s_r=+1)L_r−sum_(s_r=−1)U_r,
    U_g=sum_(s_r=+1)U_r−sum_(s_r=−1)L_r.            (3)

These operations are done separately for the two directions. In a
rigorous implementation the sums must be outward-rounded. The proposal
uses ordinary MATLAB arithmetic and is not itself a directed certificate.

Split each shared d_g=dplus_g−dminus_g with nonnegative parts. The same
proved signed inequality gives

    d_g*Psi_g <= dplus_g*U_g−dminus_g*L_g.

Together with the origin curvature and nonnegative eta_x,eta_y this yields
valid sufficient upper-curvature rows when the intervals are valid. The
original full positive vertex convolution, target density, anisotropic
ell_x²/8 and ell_y²/8 factors, objective mass and D sampling are unchanged.
The exact source checker remains unchanged: it reconstructs all actual
physical density jumps from the original density vector, independently of
GD,SG or the solver's shared auxiliary variables.

This proves the exact function identity and safe interval interface.
It does not claim that floating group-summation order, ordinary interval
bounds or shared auxiliary constraints give precisely the same numerical
LP feasible set or optimum as V2. Any numerical optimum comparison needs
its actual run and is still separate from a general C2 bound.

## 4. Prepared grid counts and zero-group handling

The exact geometry ledger `computation/r065_group_map_exact.py` actually
ran with Python3.15.0, exit0, about0.0105s. Its output is
`computation/results/r065-tensor-occupation/exact-group-map.json` and binds
the exact grid SHA53cc33d7.... It finds:

- 45 axis intervals,2025 physical cells,276 density orbits;
- 4140 physical interfaces and529 canonical groups;
- group multiplicities4 through8;
- no identically zero physical jump rows for this particular grid;
- 1104 target subrectangles at subdiv2;
- variable count10765 before grouping,3543 after grouping.

The code nevertheless handles an identically zero physical row by omitting
its group ID and sign contribution. A nonzero symbolic group whose actual
optimized density difference happens to vanish remains a valid group;
it is not confused with a zero row in the symbolic map.

## 5. New code, solver controls and optional export

`computation/cauchy_occupation_tensor_jump_lp_v3.m` has signature

    (h,grid_input,subdiv,output_root,run_mode,export_sparse).

run_mode is solve(default),diagnose orlambda-only, as in V2. The optional
export_sparse flag defaults false. If true, the PRE-solver numeric matrices
are saved as `tensor-lp-matrices.mat` in the output directory using MATLAB
v7 format, so SciPy loadmat can read sparse matrices if root chooses a
backup solver. This is a temporary computational product, not a Git
payload. Its contents include Amat,bv,Aeq,beq,cost,lower_bounds, the exact-
integer-valued maps, float geometry and group pairs. No external solver is
launched by this file or this author.

The variable order is(base lambda,a),dplus_groups,dminus_groups,eta_x,eta_y.
The equality block is[GD,−I,+I,0]. Group interval arrays use signed physical
sparse maps: L=physicalL*Positive−physicalU*Negative, and U reversesL/U in
negative contributions. The pre-solver JSON records group counts,
zero-interface count,multiplicity ranges and complete matrix scale
summaries. The full status/failure JSON is retained from V2.

For the first root-owned comparison, h=.075,the exact grid,subdiv2,both
reported tolerances10^-9 and MaxIterations20000 remain unchanged. There is
no new mesh, changed source condition, automatic row scaling or larger
iteration cap in V3. Root owns actual MATLAB launch and prompt exit. At
preparation V3 has not been run by the author; no optimized candidate,
source/full-D certificate or new C2 bound is claimed.
