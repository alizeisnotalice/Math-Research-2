# Numerical Root obstacle for the narrow Gaussian hole

Take \(\eta(dx)=(1-\epsilon/2)g(x)\mathbf1_K(x)dx\), \(K=[2,2+\epsilon^{3/2}/100]\), with \(\mu=\gamma-P_\epsilon\eta\) and \(\lambda=\gamma-\eta\). Use their unnormalized potentials \(u_\rho=-\int |x-y|\rho(dy)\). The measures have equal mass, and normalizing both rescales \(w\) and its sink by the same constant, leaving contact times unchanged.

The starting gap was computed from the positive integral
\[
 w(0,x)=\int_K\int_0^\epsilon p_t(x-y)\,dt\,\eta(dy),\qquad
 \int_0^\epsilon p_t(r)\,dt
    =2\sqrt\epsilon\int_0^\infty v\,g(v+|r|/\sqrt\epsilon)\,dv .
\]
The last identity follows from the Gaussian formula for
\(\mathbb E|r+\sqrt\epsilon G|-|r|\), using symmetry. The run integrated the positive kernel by Gauss–Legendre at \(v\le14\), and split the \(y\)-integral at \(x\) when \(x\in K\). The audit independently compares selected initial values with 55-digit quadrature of the Gaussian closed formula. The exact continuum identity \(\frac12w_{xx}(0,x)-\lambda(x)=-\mu(x)\) away from the layer endpoints provides an initial monotonicity check.

The spatial grid is aligned to the two layer endpoints, is uniformly fine immediately around \(K\), and grows in spacing toward the exterior. The sink uses dual-cell Gaussian mass divided by its width, split at both layer endpoints. Continuation obeys \(w_t=\frac12w_{xx}-\lambda\) with \(w=0\) at contact. Implicit Euler solves the discrete complementarity system by active-set tridiagonal solves. The explicit comparison instead advances one CFL-controlled finite-volume heat step followed by positive-part projection. The initial gap and the sink are shared between the two solvers. The implicit run compared different half-widths, zero boundary values and boundaries held at the initial gap. At the observed stage the central contact values were insensitive to the tested wide-box boundary condition within a time step.

## Contact times

The intervals in this table bracket the last positive and first zero numerical samples of their particular implicit runs. They represent the time step, and spatial and time-discretization sensitivity was assessed separately. Times are in microseconds. The first row uses 256 subdivisions inside \(K\), 3,747 total nodes and time step \(0.67611\) nanoseconds. The second uses 192 subdivisions in the larger \(K\), 3,899 nodes and time step \(4.54168\) nanoseconds.

| \(\epsilon\) | \(x_L-2|K|\) | Left endpoint | Center | Right endpoint | \(x_L+3|K|\) |
|---:|---:|---:|---:|---:|---:|
| \(1/8\) | [116.414695,116.415372] | [117.027251,117.027928] | [117.039421,117.040098] | [117.027251,117.027928] | [116.740580,116.741257] |
| \(1/4\) | [432.350008,432.354551] | [437.364025,437.368568] | [437.454859,437.459402] | [437.368567,437.373110] | [435.937937,435.942480] |

The fine explicit comparison at \(\epsilon=1/8\) gave center [117.037556,117.037822] microseconds, and the coarser wider explicit mesh gave [437.452073,437.452221] microseconds. The spatial and time convergence diagnostics support a small interior contact delay relative to both edges. The center time in the \(1/8\) case is about \(0.000936\epsilon\). No tested contact reopened. The largest absolute implicit complementarity defect was \(1.09\times10^{-18}\). These are numerical findings on the specified obstacle problem, without a directed PDE discretization certificate.

## Reproduction

Use **python run.py --mode main** for the twelve base cases, **python run.py --mode fine** for the finer cases, **python explicit.py** for the explicit comparison, and **python audit.py** to check source initialization and exits. The JSON case files contain the requested offsets, coordinate values, crossing samples, initial gaps, target densities, and runtime diagnostics. All calculation workers have exited.
