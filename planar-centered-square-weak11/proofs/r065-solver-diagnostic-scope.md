# R65 solver failure and deterministic feasibility diagnostic scope

2026-11-15. This note adds diagnostic scope without changing the frozen
nonuniform tensor theorem, original proposal or exact source checker.

Root actually ran the original proposal at h=.075 on the prepared grid,
subdiv2. MATLAB PID66280 exited1 and was verified absent after80.9316s.
The only saved stdout classification is the `assert(ef>0)` failure after
linprog returned. Original code did not persist ef or lpout; their actual
values are UNKNOWN and must not be reconstructed from runtime. In
particular, this record does not establish infeasibility, iteration limit,
numerical breakdown, or a lower bound on the mathematical kernel family.
The original failed source bytes are frozen SHA5f00c7f4..., with exact grid
SHA53cc33d7...; stdout and receipt are preserved byte-for-byte under
`computation/results/r065-tensor-occupation/v1-h0075-failed/`.

Root separately executed the exact nonuniform tiny self-test successfully;
the result is copied from its actual data-volume cache to
`computation/results/r065-tensor-occupation/exact-nonuniform-self-test.json`.
This verifies the small geometry/source checks, not an optimized candidate.

## 1. A mathematically legal lambda-only kernel always exists

For f=0 and lambda>0, the source is lambda*p_h>0 on all of R2 independently
of the cell partition. K=lambda*Q_h is nonnegative and integrable. For a
compact diamond D, positivity and continuity of p_s on
D times[h/2,h] give

    Q_h(x) >= 2*pi*integral_[h/2,h] p_s(x)ds
            >= pi*h*min_(x in D,s in[h/2,h])p_s(x) > 0.

Hence inf_D Q_h>0 and some finite lambda makes K>=1 on D. This establishes mathematical feasibility without claiming a
computed sharp lambda or an improved comparison cost.

The finite LP uses additional sufficient curvature rows and sampled D.
Those sufficient rows need not admit a legal lambda-only kernel: an overly
large curvature remainder can exceed a positive vertex value. Failure of
such a row is a conservative-enclosure diagnostic, not source negativity.
Conversely, a positive floating row witness is evidence about the saved
finite formulation; ordinary doubles alone are not a rigorous global
feasibility certificate.

## 2. New independent diagnostic implementation

`computation/cauchy_occupation_tensor_jump_lp_v2.m` leaves the original code
unchanged. It has the same first four arguments plus optional run_mode:

- `lambda-only`: only exact-grid-to-double geometry, origin curvatures,
  four vertex coefficients and coarse D values are evaluated. It does not
  build the large cell/interface matrices and does not call linprog.
- `diagnose`: full matrices and scale summaries are built and saved, but
  linprog is not called.
- `solve` (default): diagnostics are saved BEFORE linprog, and the complete
  returned ef, lpout, objective, exception and returned vector are written
  to `tensor-lp-solver-status.json` before any failure is raised.

The lambda-only unit vector has lambda1, f0, all jumps0, and per-subrectangle
eta_x=max(0,ox), eta_y=max(0,oy). Its unit source-row residual is

    p_h(vertex)−source_margin*max(1,|p_h(vertex)|)
      −ell_x²*eta_x/8−ell_y²*eta_y/8.

The diagnostic saves its worst rectangle, vertex, signed origin curvature
bounds, anisotropic errors and minimum residual. If the coarse D minimum
qmin is positive, lambda=(1+10^-6)/qmin is a concrete sampled-D normalized
candidate. Multiplying the unit source row by this lambda gives the source
row margin; eta scales with lambda and jump equalities remain zero. This
allows root to separate a conservative-row problem from a solver issue
before deciding another costly solve.

The full matrix diagnostic records dimensions, nnz, signed/nonzero ranges,
row/column infinity-norm ranges and zero-row/column counts, plus cost and
right-side ranges. Solver controls remain MaxIterations20000 and both
reported tolerances10^-9; no automatic scaling, extra mesh, relaxed source
constraint or higher iteration limit is applied in this diagnostic version.

If no solver result returns before an external time bound, the pre-solver
files still exist; the wrapper's timeout remains separate evidence. If
linprog throws, the caught report/stack is saved. If ef<=0 or empty, the
code raises an error only AFTER saving its status file. A returned partial
iterate is saved for diagnosis and is not labeled a certified kernel.

At preparation the new diagnostic version is unrun by the author. Root
chooses the next mode and owns all MATLAB processes. No new grid, source
certificate, full-D certificate or improved C2 bound is claimed here.
