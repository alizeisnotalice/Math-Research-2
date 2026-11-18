# Sharp accumulation of the global Hermite tail bound

Let \(h=f-1\) be bounded and satisfy the full Cantelli transform identity under standard Gaussian measure \(\gamma\). Put \(\delta=\|h\|_\infty>0\). Write \(a_j=\langle h,\mathrm{He}_j/\sqrt{j!}\rangle_\gamma\). Expanding
\[
 \mathbb E_\gamma[e^{tX-t^2/2}(e^{-t^2h/2}-1)]=0
\]
at zero and comparing degree \(j+2\) gives the exact relation
\[
 a_j=\sum_{s=1}^{\lfloor j/2\rfloor}
  \frac{(-1)^{s+1}}{2^s(s+1)!}
  \sqrt{\frac{j!}{(j-2s)!}}\,
  \left\langle h^{s+1},\frac{\mathrm{He}_{j-2s}}{\sqrt{(j-2s)!}}\right\rangle_\gamma .
\]
The expansion is justified by bounded \(h\) and Gaussian exponential moments. For each \(s\) the vector on \(j\le K\) maps a truncated Hermite coefficient vector of \(h^{s+1}\) into distinct coordinates by raising index by \(2s\). Its \(\ell^2\) norm is at most \(\sqrt{K!/(K-2s)!}\|h^{s+1}\|_2\). The triangle inequality and \(\|h^{s+1}\|_2\le\delta^s\|h\|_2\) give
\[
 \frac{\|P_{\le K}h\|_2}{\|h\|_2}\le B_K(\delta),
 \qquad B_K(\delta)=\sum_{s=1}^{\lfloor K/2\rfloor}
     \frac{(\delta/2)^s}{(s+1)!}\sqrt{\frac{K!}{(K-2s)!}} .
\]

## Optimality from the cumulative constraints alone

Set \(u_K=\min(1,B_K^2)\) for \(K\ge0\), \(u_{-1}=0\). The normalized spectral law \(p_j=|a_j|^2/\|h\|_2^2\) has cumulative probabilities at most \(u_K\). Summing each nonnegative tail gives
\[
 \frac{\int |h'|^2\,d\gamma}{\int h^2\,d\gamma}
    =\sum_{K\ge0}\sum_{j>K}p_j
    \ge S(\delta):=\sum_{K\ge0}(1-B_K(\delta)^2)_+ .
\]
When the Sobolev energy is infinite the lower bound still holds. Since \(B_K\) is increasing in \(K\) and eventually exceeds one, \(q_j=u_j-u_{j-1}\) is a probability law, with finite support. It simultaneously saturates all cumulative caps and realizes \(S(\delta)\) as its first moment. This proves optimality over abstract spectral laws satisfying these bounds. It does not provide an exact graph or a bounded realization of that law.

## Small-amplitude asymptotics

For \(\delta\downarrow0\), set
\[
 \phi(x)=\sum_{s\ge1}\frac{(x/2)^s}{(s+1)!}
       =\frac{e^{x/2}-1-x/2}{x/2},
 \qquad x_*>0,\quad \phi(x_*)=1.
\]
Then
\[
 S(\delta)=C_*/\delta+D_*+o(1),
 \quad C_*=\int_0^{x_*}(1-\phi(x)^2)\,dx,
 \quad D_*=2(x_*-1)-2\int_0^{x_*}x\phi'(x)^2\,dx.
\]
The constants are approximately \(C_*=1.84746413904031\), \(D_*=1.66332055722075\), with \(x_*=2.51286241725234\).

For \(K\ge1\) the normalized falling factorial, set to zero for \(2s>K\), satisfies
\[
 K^{-s}\sqrt{K!/(K-2s)!}
 =1-\frac{s(2s-1)}{2K}+O(s^4/K^2).
\]
If \(s^2/K\) is bounded below the estimate follows by enlarging the constant, as the left side lies between zero and one. Otherwise expand \(\frac12\sum_{j<2s}\log(1-j/K)\) and exponentiate; its error and the square of its first-order term are \(O(s^4/K^2)\). For \(x=\delta K\) in a fixed bounded interval, summing the resulting series gives
\[
 B_K=\phi(x)-\frac{\delta}{2}\big(\phi'(x)+2x\phi''(x)\big)+O(\delta/K).
\]
The last error comes from \(\sum_{s\ge1}s^4(x/2)^s/(s+1)!=O(x)\) on such an interval. The cap is already one once \(K\) is larger than a constant times \(\delta^{-1}\), by just its first summand. Thus the sum over relevant \(K\) of the squared-cap expansion error is \(O(\delta\log(1/\delta))\). At \(K=0\) one treats the exact term one separately.

The positive-part operation changes the linearized expansion only within \(O(\delta)\) of the simple root \(x_*\), contributing \(O(\delta)\) to the sum there. The trapezoidal rule for the compactly supported continuous \((1-\phi^2)_+\), piecewise smooth, gives endpoint term \(1/2\). The first-order correction is
\[
 \int_0^{x_*}\phi(x)(\phi'(x)+2x\phi''(x))\,dx .
\]
After integration by parts the endpoint \(1/2\) and the integral's internal \(-1/2\) cancel. The relation \(e^{x_*/2}=1+x_*\) gives \(\phi'(x_*)=1-1/x_*\), yielding \(D_*\).

The worker code brackets each selected finite \(S(\delta)\) using exact integer square roots of scaled falling-factorial integers, bounded positive Taylor tails, and rational endpoint sums. Floating grids are a separate diagnostic.

## A simple bound at every amplitude

For \(K\ge1\), each successive pair in the falling factorial square root is bounded by \(K-1/2\), since \(\sqrt{x(x-1)}\le x-1/2\) at all integers \(x\ge1\). Thus \(B_K\le\phi(\delta(K-1/2))\). For the decreasing function \(q=(1-\phi^2)_+\), the deficit of the midpoint rule on any cell \([l,r]\) of length \(\delta\) is bounded above by \(\tfrac{\delta}{2}[q(l)-q(r)]\). Summing over cells gives a total deficit at most \(\delta/2\). Including the exact \(K=0\) term one yields the finite-amplitude bound
\[
 S(\delta)\ \ge\ C_*/\delta+1/2\qquad(\delta>0).
\]
