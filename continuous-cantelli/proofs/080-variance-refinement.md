# Bounded exact Gaussian variance graphs impose a low-frequency Hermite constraint

Let \(X,Z\) be independent standard Gaussians, \(a\) be a bounded nonnegative measurable function, and \(X+\sqrt{a(X)}Z\) be centered Gaussian. Its variance is \(1+\mathbb E a(X)\). Write \(h=a-\mathbb Ea\), \(\epsilon=\|h\|_{L^\infty(\gamma)}\), \(\sigma=\|h\|_{L^2(\gamma)}>0\), and \(p_k=\mathrm{He}_k/\sqrt{k!}\). Put \(c_k=\langle h,p_k\rangle_\gamma\).

Conditioning on \(X\) gives the identity
\[
  \mathbb E\!\left[e^{tX-t^2/2}e^{t^2h(X)/2}\right]=1.
\]
Both factors and their product can be expanded on compact complex \(t\)-sets under the expectation, dominated by a fixed multiple of \(e^{R|X|}\) when \(|t|\le R\). The coefficient of \(t^{k+2}\) yields
\[
 c_k=-\sum_{\ell=1}^{\lfloor k/2\rfloor}
 \frac{2^{-\ell}}{(\ell+1)!}\sqrt{\frac{k!}{(k-2\ell)!}}\,
 \langle h^{\ell+1},p_{k-2\ell}\rangle_\gamma.
\]
In particular \(c_0=c_1=0\) and \(c_2=-\sigma^2/(2\sqrt2)\). On the finite vector indexed by \(0\le k\le K\), apply Minkowski to the sum over \(\ell\) and Bessel to each shifted Hermite coefficient vector. Its diagonal weights are bounded by \(\sqrt{(K)_{2\ell}}\) and \(\|h^{\ell+1}\|_2\le\epsilon^\ell\sigma\). Therefore
\[
 \sum_{k=0}^K |c_k|^2/\sigma^2
 \le B_K:=\min\left\{1,S_K(\epsilon)^2\right\},\qquad
 S_K(\epsilon)=\sum_{\ell=1}^{\lfloor K/2\rfloor}
 \sqrt{(K)_{2\ell}}\frac{(\epsilon/2)^\ell}{(\ell+1)!}.
\]
Here \((K)_j=K!/(K-j)!\).

## A variance-ratio refinement

Set \(r=\sigma/\epsilon\). For \(K\ge3\), estimate the vector of modes \(3\le k\le K\) by the same argument. If \(\ell\) is odd, Jensen gives \(|\langle h^{\ell+1},p_0\rangle|^2\ge\sigma^{2\ell+2}\). This constant component of the input then enters the output with squared weight \(d_\ell\), where \(d_1=0\) because mode 2 has been removed and \(d_\ell=(2\ell)!\) for odd \(\ell\ge3\). For even \(\ell\) retain the unrefined term. Consequently a refined cap is
\[
 \widetilde B_K =\min\left\{B_K,(\epsilon r)^2/8+\widetilde S_K^2\right\},
\]
where \(\widetilde S_K\) is \(S_K\) with each odd-\(\ell\) weight replaced by
\(\sqrt{(K)_{2\ell}(1-r^{2\ell})+d_\ell r^{2\ell}}\). At \(K=0,1\) take \(\widetilde B_K=0\), and at \(K=2\) take \(\widetilde B_K=(\epsilon r)^2/8\). Use this equality only for feasible \(\epsilon,r\). For a constant graph \(\sigma=0\), normalized energy is undefined.

## Finite and asymptotic energy bounds

With \(w_k=|c_k|^2/\sigma^2\), mass is one. For every nondecreasing nonnegative \(f\) on the nonnegative integers, summing the nonnegative tail contributions gives
\[
 \sum_{k\ge0} f(k)w_k\ge f(0)+\sum_{K\ge0}
 (f(K+1)-f(K))(1-B_K).
\]
The corresponding refined inequality uses \(\widetilde B_K\). For \(f(k)=k\), the left side equals \(\|h'\|_2^2/\sigma^2\) for Gaussian Sobolev \(h\), with the spectral sum interpreted in \([0,\infty]\) in general. For \(f(k)=k^2\) it is the squared norm of the Ornstein–Uhlenbeck generator divided by \(\sigma^2\) on its domain. A single cutoff gives only \((K+1)(1-B_K)\) for \(f(k)=k\). Since \(B_K\) is nondecreasing, bounded by one and tends to one, the distribution with CDF \(B_K\) attains the all-prefix lower bound in the relaxation where only the prefix inequalities are imposed. It need not arise from a graph.

Fix \(r\) and let \(K\epsilon/2\to x\) as \(\epsilon\to0\). Each individual summand tends to its counterpart in
\[
 s_r(x)=\sum_{\ell=1}^\infty
       \frac{x^\ell}{(\ell+1)!}
       \begin{cases}
       \sqrt{1-r^{2\ell}},&\ell\ \text{odd},\\
       1,&\ell\ \text{even}.
       \end{cases}
\]
The summands are dominated uniformly for \(K\epsilon/2\) in a compact set by the corresponding exponential series. At \(x=0\), both sums tend to 0. Thus the rescaled caps converge on compact sets by the same majorant and coefficient convergence. In the base case \(s_0(x)=(e^x-1-x)/x\) and the unique positive solution of \(s_0(x)=1\) solves \(e^x=1+2x\). At \(r=1\), \(s_1(x)=\sinh(x)/x-1\) and the positive cutoff solves \(\sinh x=2x\). The refined cap uses a minimum with the base cap. In the limit this minimum coincides with \(\min(1,s_r(x)^2)\) since \(0\le s_r(x)\le s_0(x)\).

If \(x_r\) denotes the unique positive cutoff for \(s_r\), then for \(p\ge1\) and a sequence of exact graphs with \(\epsilon\to0\) and ratio \(r\) fixed,
\[
 \liminf_{\epsilon\to0}\epsilon^p\sum_{k\ge0} k^p w_k
 \ \ge\ 2^p p\int_0^{x_r} x^{p-1}(1-s_r(x)^2)\,dx .
\]
For \(p=1\), the individual cutoff method has largest limiting lower coefficient \(2\max_{0\le x\le x_r}x(1-s_r(x)^2)\). For the integral inequality one may sum only to \(2x_r/\epsilon\) and apply Riemann convergence. The same lower bound holds along profiles with \(\sigma/\epsilon\to r\) by continuity and the uniform majorant. These implications all start from the exact graph assumption.

## Computational verification

Run the deterministic exact-fraction finite scan and a separate mpmath check with:
\[
 \texttt{python experiment.py --duration 55 --profiles 240}.
\]
The duration gate determines the number of computed rows. Each shard ends on a complete row. A rational upper bound on each square root uses 80 fractional bits and an integer square root. The omitted positive tail after term 22 is bounded by a geometric series with ratio \((K\epsilon/2)/24\). Every emitted finite lower bound is rounded downward to 16 decimal places. Finite evaluations stop at first refined saturation or a preset finite cutoff. An unsaturated cutoff still gives a lower bound, since only nonnegative contributions were omitted. The sampled independent 49-digit computation directly evaluates factorial terms through order 32. The separately computed limiting coefficients in the summary are mpmath numerical calibrations.

## Quantitative exclusion of nearly two-valued exact graphs

For any bounds \(m\le h\le M\) with \(m<0<M\), let \(d=M-m\) and choose \(t=(1+i)\sqrt{2\pi/d}\), so that \(t^2h/2=2\pi i h/d\). Put \(z=e^{2\pi i m/d}=e^{2\pi i M/d}\). Complex MGF exactness and \(\mathbb E e^{tX-t^2/2}=1\) give
\[
 1-z=\mathbb E\!\left[e^{tX-t^2/2}
    \left(e^{2\pi i h/d}-z\right)\right].
\]
The phase difference has absolute value at most \((2\pi/d)\min(h-m,M-h)\), and \(\min(h-m,M-h)^2\le(h-m)(M-h)\). Apply Cauchy–Schwarz and \(\mathbb E|e^{tX-t^2/2}|^2=e^{4\pi/d}\). This proves the strictly positive Bhatia–Davis variance gap
\[
 -mM-\sigma^2\ \ge\
   \frac{d^2}{\pi^2}\sin^2\!\left(\frac{\pi m}{d}\right)e^{-4\pi/d}.
\]
In particular, a nonconstant exact graph cannot have only two values. With \(m=-\epsilon\), \(M=\epsilon\) this yields \(1-r^2\ge (4/\pi^2)e^{-2\pi/\epsilon}\). The coefficient listed at \(r=1\) is the limit of the energy bound as \(r\to1\).
