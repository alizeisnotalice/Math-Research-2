# Sign of the fractional heat test on a biased Cantor cluster

For an integer \(m\ge9\), put \(\rho=m^{-2}\). The similarities have ratio \(\rho\) and left endpoints \(a_0=0\) and \(m-1\) equally spaced points \(a_1,\ldots,a_{m-1}\) on \([9/10,1-\rho]\), including the two endpoints. They are strongly separated, since \((m-1)\rho<1/10\). Their uniform Bernoulli probability is \(\mu\), the law of \(Y\). The dimension is \(1/2\).

At variance \(b=\rho^{2n}/D^2\) evaluate
\[
T_\kappa=\int_0^\infty[(p_b*\mu)(0)-(p_{b+t}*\mu)(0)]t^{-1-\kappa}\,dt,\qquad \kappa\in\{1/4,1/2\}.
\]
The positive normalization \(S_n=\sqrt{2\pi}m^n b^{\kappa+1/2}T_\kappa\) obeys
\[
S_n=E J_\kappa(DY)+\sum_{r=1}^n m^{r-1}\sum_{i=1}^{m-1}E J_\kappa(D\rho^{-r}(a_i+\rho Y)),
\]
by conditioning on the first nonzero address digit before depth \(n\). Integration by parts in \(t/b\), then \(q=(1+t/b)^{-1}\), gives
\[
J_\kappa(z)=\frac1{2\kappa}\int_0^1q^{a-1}(1-q)^{-\kappa}(1-z^2q)e^{-z^2q/2}\,dq,\qquad a=\kappa+1/2 .
\]
The heat difference is integrable at both ends. For the centre term, the interval certificate below encloses the full support of each cylinder group, including the Cantor suffix.

## Negative distant branches

For \(z\ge9/2\) use monotonicity of \(w(q)=(1-q)^{-\kappa}\). On either side of \(q_0=1/z^2\) the factor \(1-z^2q\) changes sign, so the integral with \(w\) is at most \(w(q_0)\) times the unweighted integral. The latter equals
\[
2e^{-z^2/2}-2\kappa\int_0^1q^{a-1}e^{-z^2q/2}dq .
\]
Here integration by parts uses \(q^ae^{-z^2q/2}\). The integral is at least \(e^{-1/2}z^{-2a}/a\), by restriction to \(0\le q\le z^{-2}\). The bracket is strictly negative for \(z\ge9/2\), since \(z^{2a}e^{-(z^2-1)/2}<\kappa/a\) there. For example it is at most \(20e^{-19/2}<1/3\) by \(2a\le2\) and monotonicity for \(z^2\ge20\). Because \(w(q_0)\ge1\), it follows that
\[
J_\kappa(z)\le e^{-z^2/2}/\kappa-\int_0^1q^{a-1}e^{-z^2q/2}dq<0 .
\]

For an outer absolute bound split the integral for \(J_\kappa\) at \(q=1/2\). On the lower half use \((1-q)^{-\kappa}\le3/2\), \(|1-z^2q|\le1+z^2q\), and extend the resulting gamma integrals to infinity. With \(\Gamma(a)<2\) this part is bounded by \(36z^{-2a}\). On the upper half the absolute contribution is at most \(8(1+z^2)e^{-z^2/4}\), by \(B(a,1-\kappa)<4\). For \(z^2\ge20\), \(z^2(1+z^2)e^{-z^2/4}\) decreases and is at most \(420e^{-5}\). In total, \(|J_\kappa(z)|\le100z^{-2a}\) on this range. Thus all outer terms are negative, and the magnitude of their infinite sum is at most
\[
100(0.9D)^{-2a}\frac{m^{-(4\kappa+1)}}{1-m^{-(4\kappa+1)}}.
\]
This bound also controls every finite \(n\).

## Uniform negativity for \(m\ge256\)

For \(D\in[5,15]\) and a right first cylinder, \(z=DY\in[9/2,15]\). Its negative contribution has a uniform explicit magnitude. Write \(\epsilon=(3/8)^{10}\), so \(e^{-z^2/2}<\epsilon\). When \(\kappa=1/2\), the integral in the preceding upper bound is at least \(2(1-\epsilon)/225\). When \(\kappa=1/4\), substitute \(u=z^2q/2\) and retain \([0,1]\) and \([1,2]\). Using \(e<11/4\), \(2^{-1/4}>4/5\), and \((2/225)^{3/4}>1/35\), the integral is larger than
\[
\frac1{35}\left(\frac43\frac4{11}+\frac45\frac{16}{121}\right).
\]
Subtract \(4\epsilon\) from each lower bound to cover the exponential term. The resulting positive magnitudes \(\delta_{1/2},\delta_{1/4}\) satisfy \(255\delta_{1/2}>2\) and \(255\delta_{1/4}>4\) by exact rational arithmetic.

On the left cylinder \(J_\kappa(z)\le J_\kappa(0)\). For \(\kappa=1/2\) this value is 2. For \(\kappa=1/4\) it is less than 4, by splitting the symmetric beta integral into halves and using \(2^{1/4}\int_0^{1/2}q^{-1/4}dq=(4/3)/\sqrt 2<1\). The left cylinder has mass \(1/m\), so the bound on the right dominates for \(m\ge256\). With all outer terms negative, \(T_\kappa<0\) for both orders, every \(D\in[5,15]\), and every \(n\ge0\).

## Certified mixed signs

The interval check integrates the kernel over rational subboxes using \(q=u^4/2\) and \(q=1-u^4/2\) on the two halves, which remove the singular endpoint powers. It encloses each group of complete first-level cylinders, then adds the conservative bound for all omitted outer shells. It proves, uniformly in \(n\ge0\), that at \(m=64,D=10\),
\[
S_n(1/4)\in[-0.023798,-0.0131709],\qquad
S_n(1/2)\in[0.0074225,0.0097355].
\]
At \(m=128,D=15\) the corresponding enclosures are \([-0.0148327,-0.0092823]\) and \([0.0051458,0.0061918]\). This computation concerns the prescribed cluster test.
