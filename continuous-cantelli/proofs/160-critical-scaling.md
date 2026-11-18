# Critical smoothing radius of the normalized fat Cantor set

Start with \(I=[-1/2,1/2]\). At stage \(n\ge1\) remove the open middle fraction \(\theta_n=2^{-n}/256\) from every surviving interval. Let \(K_n\) be the stage set and \(K=\bigcap_n K_n\). Write \(u_v(x)=(p_v*1_K)(x)\) and \(C=p_{7/8}(2)/p_1(2)=\sqrt{8/7}e^{-2/7}\). Here \(v\) is Gaussian variance.

## Unique maximizing location near the critical noise

The removed set has length at most \(\sum_n\theta_n=1/256\). On \(0\le x\le1/2\) the full interval convolution satisfies
\[
(p_v*1_I)''(x)=-[(x+1/2)p_v(x+1/2)+(1/2-x)p_v(1/2-x)]/v .
\]
The first term alone has magnitude at least \(e^{-1/(2v)}/[2\sqrt{2\pi}v^{3/2}]\). On the other hand \(\|p_v''\|_\infty=1/(\sqrt{2\pi}v^{3/2})\), so the hole integral changes curvature by at most \(1/[256\sqrt{2\pi}v^{3/2}]\). For \(v\ge1/8\) the strict inequality \(e^{-4}/2>1/256\) gives \(u_v''<0\) throughout \(I\). Evenness and the sign of \(p_v'(x-y)\) for \(x>1/2\), \(y\in I\), imply the unique global maximum is at \(x=0\). The proof applies equally to every stage \(K_n\).

For \(v>0\) the global maximum is attained, since the Gaussian convolution vanishes at infinity. For \(h>0\), \(u_{v+h}=p_h*u_v\) is strictly below \(\max u_v\) at every point: \(u_v\) is continuous and its submaximum set contains a nonempty open interval, and \(p_h\) is strictly positive everywhere. Hence \(\max u_v\) decreases strictly with \(v\). The certified signs at \(v=1/8\) and \(v=1/6\) give a unique critical variance in that interval.

## Controlled high-precision computation

Let \(L_n=|K_n|/2^n\) denote the leaf length. The uniform probability measure on \(K\) is the law of the series \(\sum_{k\ge1}a_k\xi_k\) with independent symmetric Bernoulli signs and \(a_k=L_{k-1}(1+\theta_k)/4\). Every leaf has the same length and future retained fraction, so the stage choices have equal conditional probabilities. The series remainder after stage \(N\) is centered conditionally on that leaf and lies in \([-L_N/2,L_N/2]\).

The interval computation obtains the moments of the truncated series by multiplying the even exponential generating series of its Bernoulli terms. It applies the alternating expansion of \(e^{-X^2/(2v)}\) at the center, where \(X^2/(2v)\le1\), forms the curvature from the same moments, and encloses both analytic remainders. Bounds \(\|p_v\|_\infty\le2,\|p_v''\|_\infty\le16,\|p_v^{(4)}\|_\infty\le4000\) hold for \(v\ge1/8\). The last follows from the fourth Gaussian derivative and \(z^4e^{-z^2/2}\le16\), \(z^2e^{-z^2/2}\le2\). With \(L_N\le2^{-N}\), the conditional Taylor errors in the value and curvature are at most \(2L_N^2\) and \(500L_N^2\). The total remaining mass difference is at most \(2^{-N}/256\). The script uses \(N=100\), expansion degree 42, and rational bisection endpoints.

The enclosure for \(v_*\) has width below \(3.37\times10^{-29}\). The critical curvature is about \(-2.99738371237800904523\), so the peak's variance slope is \(u_{v_*}''(0)/2\).

## Stage convergence and mass correction

Let \(A_n=|K_n|\), \(A=|K|\), \(q_n=A/A_n\), and let \(m_n(v)=(p_v*1_{K_n})(0)\), \(m(v)=u_v(0)\). The remaining gaps inside every leaf have a common relative mass \(1-q_n\), and their conditional law is symmetric about the center of that leaf. The whole leaf law is also symmetric. Taylor expanding the fixed Gaussian kernel at each leaf center therefore gives
\[
|m(v)-q_nm_n(v)|\le A_n(1-q_n)L_n^2\|p_v''\|_\infty/4=O(8^{-n}/256)
\]
uniformly for \(v\in[1/8,1/6]\). As \(1-q_n=2^{-n}/256+O(4^{-n}/256^2)\), if \(m_n(v_n)=C\) then the mean value theorem gives
\[
2^n(v_n-v_*)\longrightarrow \frac{C}{256\,|m'(v_*)|}
=\frac{2C}{256\,|u_{v_*}''(0)|}.
\]
Solving the mass-rescaled equation \(q_n m_n(v)=C\) changes the root approximation error to \(O(8^{-n}/256)\). The analytic rate explains why finite-stage Float32 errors eventually dominate the acceleration.
