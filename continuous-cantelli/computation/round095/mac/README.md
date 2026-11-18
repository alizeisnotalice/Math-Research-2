# R088 moving-tail and endpoint budget

The computation uses the reflected halfline survival bound with \(T,B>0\) and adjoining low level \(H_{\rm low}=1+T\). Write \(p_T\) for the centered Gaussian density of variance \(T\), \(\phi=p_1\), and \(K=2Bp_T(B)/T\). The cap after multiplication by \(\lambda^2e^{-B\lambda}\) is \(KJ_{\rm cap}\), with

\[
J_{\rm cap}=\int_0^\infty v e^{-v}\exp[-v^2/(2T\lambda^2)]\,\mathrm{sinhc}(Bv/(T\lambda))\,dv
=1+\frac{B^2/T^2-3/T}{\lambda^2}+O(\lambda^{-4}).
\]

For the far profile \(H_{\rm low}-f(r)=2B/r-(\beta\log r+C)/r^2\), set \(L=\log\lambda\), \(r=\lambda+z\), \(z=y-B\), \(w=z/\lambda\), and \(Q=\beta L+C\). The recentered far residual relative to the leading Gaussian is

\[
R_\lambda(y)=\frac{Bz^2}{\lambda+z}-\frac{\beta\log(1+w)-Qw(2+w)}{2(1+w)^2}.
\]

The far contribution at leading order is \(e^{(B^2-C)/2}\lambda^{2-\beta/2}\), with first relative correction \(B(B^2+1+\beta/2-Q)/\lambda\). The first grid uses \(C=B^2-2\log K+\delta\), so the GPU log margin is \((2-\beta/2)L-\delta/2+\log\mathbb E e^{R_\lambda}-\log J_{\rm cap}\). The expectation uses the conditional normal on \([-8,8]\) for finite comparison. The cap's correction integral uses \(0\le v\le64\).

## Joint endpoint and far-tail condition

At an edge with expansion \(f(B+u)-T=cu+du^{3/2}+O(u^2)\), \(c>0\), the scaled endpoint contribution is \(E_0+E_1/\lambda+O(\lambda^{-2})\), where

\[
E_0=2\phi(B)/c,\qquad E_1=\phi(B)\left[4/c^2-(d/2)\Gamma(5/2)(2/c)^{5/2}\right].
\]

Integrate the positive endpoint kernel on \(u\ge0\) within a shrinking window using \(v=c\lambda^2u/2\). Its integrand relative to \(e^{-v}\) has first correction \((2v/c-(d/2)(2v/c)^{3/2})/\lambda\). A logarithmic cutoff in \(v\) gives negligible remaining exponential mass and an integrated \(O(\lambda^{-2})\) Taylor error.

The endpoint and moving far windows are disjoint and contribute to the same cap. Thus \(\beta<4\) is excluded. At \(\beta=4\), set \(F_0=e^{(B^2-C)/2}\). The leading condition is \(E_0+F_0\le K\). If this is saturated and the far deficit also has \(+(\eta\log r+D)/r^3\), then \(\eta\le8B\). The first constant at \(\eta=8B\) must satisfy

\[
E_1+F_0[B(B^2+3-C)+D/2]\le0.
\]

If this first constant is also saturated and the far deficit has the additional term \(+\theta(\log r)^2/r^4+o((\log r)^2/r^4)\), then \(\theta\le-16\). Under \(z\sim N(-B,1)\) the exponent's \(L/\lambda\) component is \(4(z+B)L/\lambda\), whose mean is zero. Half its second moment contributes \(8L^2/\lambda^2\) to the Gaussian average. The only \(L^2/\lambda^2\) term from the exponent itself has coefficient \(\theta/2\). The cap and endpoint remainders are \(O(\lambda^{-2})\). A far cutoff \(|y|\le\sqrt{16\log\lambda}\) keeps both Taylor errors and discarded reference normal moments below the \(L^2/\lambda^2\) scale. Consequently the excess above the cap has leading coefficient \(F_0(8+\theta/2)\) at that scale. These necessary conditions use the two stated profile expansions and the reflected cap.

## Native runs and CPU checks

- The main scan uses geometric \(B,T\in[1/4,4]\) and \(\lambda\in[32,65536]\), plus \(\beta\in[3.25,4.75]\) and \(\delta\in[-1,1]\). All 51,212,403 samples are finite. At \(B=1/4,T=4,\lambda=32,\beta=4,\delta=0\), the Metal log margin is 0.01236509 and the independent CPU value is 0.012365060487008743.
- The joint scan allocates a far budget share \(F_0/K\in[0.1,0.9]\), samples \(\eta-8B\in[-1,1]\), sets \(d=\theta=0\), and chooses \(c,D\) to saturate the leading and first constant budgets at \(\eta=8B\). The edge cutoff \(u=128/(c\lambda^2)\) is accepted only for \(c\lambda\ge16\), yielding \(u\le8/\lambda\le1/4\) and keeping it disjoint from \(r=\lambda-B+y\) with \(|y|\le8\). Of 5,419,375 submissions, 5,330,775 were accepted and finite. At exact \(\eta=8B\), 213,214 of 213,231 accepted combinations gave positive finite-grid margins in this model. The preliminary joint run lacked the edge-window guard and is archived as a diagnostic outside these counts.
- SciPy adaptive quadrature used an independent scaled complementary-error-function formula for the cap. The 68-digit mpmath check used the unexpanded tail exponent and high-precision cap, parameterized by the same finite-window conventions. Ten replay points match exactly between local and remote numeric reports. The largest GPU versus CPU absolute log margin error is \(1.441\times10^{-6}\). The double versus 68-digit discrepancy is at most \(1.061\times10^{-11}\).
- High-precision replay at \(B=0.75,T=0.5,F_0/K=0.5\) extended \(\lambda\) to \(10^8\). At \(\theta=-20,-16,0\) the scaled log margin \(2\lambda^2\log(\mathrm{ratio})/(\log\lambda)^2\) gave \(-2.358762822,-0.358762795,7.641237310\), approaching the respective limits \(-2,0,8\). The analytic expansion supplies the limit.

The Apple M5 Pro native timings were 4.326583125 s and 0.558558333 s, from Metal command buffers. OS device utilization readings were 100% during the first run and 99–100% during the joint run. The OS readings describe system-wide load. Both runs and the two diagnostics exited code zero, and an independent post-run process lookup found all four PIDs absent. See `summary.json`, `exit-and-hardware.log`, and `remote.sha256`.

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc88-temp"
xcrun swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o laplace_scan
sh run.sh sync 41 37 33 31 33 768 1024
(
  cd joint_admissible
  xcrun swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o laplace_scan
  sh run.sh joint_admissible 25 23 13 25 29 768 1024
)
python check_cpu.py . cpu.replay.json
```
