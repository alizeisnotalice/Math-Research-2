# Cantor cylinder coupling for the complex-averaged heat kernel

Let \(\mu\) be the standard Cantor probability on \([2,3]\), let \(g=p_1\), and set \(\sigma=g\mu\). Let \(K_b(z)=p_b(z+iD)\), \(\omega=D/b\), and \(r=\sqrt b\). All numerical values are divided by the common positive factor
\[
C=g(2)^2p_b(0)e^{D^2/(2b)}.
\]
In these units the underlying kernel on \((a,x)\in[2,3]^2\) is
\[
W(a,x)=\frac{g(a)g(x)}{g(2)^2}
\exp[-(a-x)^2/(2b)]\exp[-i\omega(a-x)].
\]
The energy \(E=M_b/C=\iint W\,d\mu\,d\mu\) is positive, by the positive Fourier representation in round 185. Strict positivity follows from the continuity of the Fourier transform of the nonzero positive measure and the strictly positive frequency weight.

A depth-\(j\) parent cylinder has length \(\ell=3^{-j}\). Let \(G_j\) be the complex integral where the source measure is replaced by a uniformly chosen parent cylinder, then uniform on its middle gap of length \(\ell/3\). Let \(H_j\) instead use its midpoint. The removed gaps here are those between the two children of depth-\(j\) parents. All three computations use identical \(b,D\) and the same Gaussian weights.

For a shared parent prefix, a normalized Cantor suffix \(Y\) is in either the first or last third with equal probability and conditional means \(1/6\) and \(5/6\). For every fixed \(u\) in the central gap it follows that \(E|Y-u|=1/3\). This gives the conservative normalized Lipschitz estimate
\[
|G_j-E|,\ |H_j-E|\le \frac{\ell}{3}
  \left[2+e^{-1/2}/r+|\omega|\right].
\]
Indeed \(\sup_{[2,3]}g/g(2)=1\), \(\sup|g'|/g(2)=2\), and \(\sup|p_b'|/p_b(0)=e^{-1/2}/r\). The same bound holds for each fixed source position in the central gap, before averaging the gap.

Define
\[
U=3+4e^{-1/2}/r+1/r^2
   +2|\omega|(2+e^{-1/2}/r)+\omega^2.
\]
This bounds the normalized second \(x\) derivative of \(W\); it uses \(\sup |g''|/g(2)=3\) on \([2,3]\) and \(\sup |p_b''|/p_b(0)=1/b\). The Cantor suffix and gap variables have the same mean \(1/2\) and variances \(1/8\) and \(1/108\), respectively. Taylor expansion around the parent midpoint yields
\[
|G_j-E|\le \frac{29}{432}\ell^2 U,\qquad
|H_j-E|\le \frac{\ell^2}{16}U.
\]
The second inequality also bounds replacing one Cantor argument by its depth-\(N\) midpoint measure when \(j=N\). Thus approximating both arguments for the energy costs at most \(2\cdot3^{-2N}U/16\).

For the radially varying age \(b(u)\), the reflected gap positions \(u\) and \(1-u\) have equal age. At each such pair the average source position is the parent midpoint. The same Taylor argument for the paired law gives normalized second-order error at most \(11\ell^2 U(b(u),\omega(u))/144\), because \((u-1/2)^2\le1/36\). When reverting to the original units the prefactor \(C\) is the same for both members of this pair.

## Rigorous sign checks for coarse and fine equalization

At \(b=49/4,\ D=49\pi/2\), with resonance index \(m=0\), interval arithmetic over complete depth-six cylinders directly encloses the true infinite-Cantor energy and the true uniform gap at \(j=0\). The normalized real values satisfy
\[
0.06075<E<0.06343,\qquad
-0.04894<\operatorname{Re}G_0<-0.04573.
\]
The interval quadrature treats the full mass of each cylinder by a range enclosure, including the remaining Cantor suffix. Thus even actual gap averaging and Gaussian weighting can reverse the energy sign at the same \(b,D\). At the same \(b,D\) the depth-four Lipschitz estimate is less than \(0.03481\), yielding \(\operatorname{Re}G_4>0.02595\). These checks concern this coupling step, with the rest of the heat identity handled by the analytic comparisons in the project proof.
