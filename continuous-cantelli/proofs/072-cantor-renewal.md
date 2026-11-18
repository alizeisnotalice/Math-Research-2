# Cantor minimum sets and resonant circular Gibbs coefficients

Fix an integer \(b\ge3\), let \(C_b\subset[0,1]\) be the self-similar Cantor set with base-\(b\) digits \(0,b-1\), and let \(\nu_b\) be its equal-weight Bernoulli law. Identify \(0\) and \(1\) on the circle. Put \(d(x)=\operatorname{dist}(x,C_b)\), \(D=\log_b 2\), and \(H=(b-2)/2\). Define
\[
 Z_\lambda=\int_0^1 e^{-\lambda d(x)^\alpha}\,dx,\qquad
 \mu_\lambda(dx)=Z_\lambda^{-1}e^{-\lambda d(x)^\alpha}\,dx .
\]
The circular Gaussian output with variance \(v_0+c\,d(x)^\alpha\) has coefficient
\(e^{-2\pi^2n^2v_0}Z_\lambda\widehat\mu_\lambda(n)\) at \(n\ne0\), where \(\lambda=2\pi^2 c n^2\) and \(\widehat\mu(n)=\int e^{2\pi i n x}\mu(dx)\).

## The subquadratic endpoint limit

Write \(\epsilon=\lambda^{-1/\alpha}\). The complement at level \(l\ge1\) has \(2^{l-1}\) equal gaps of half-length \(H b^{-l}\), each with distance to \(C_b\) equal to its distance to the nearer endpoint. The mass of that level is comparable to \(2^l\min(\epsilon,Hb^{-l})\). Summing the two geometric regimes gives \(Z_\lambda\asymp\epsilon^{1-D}\) as \(\epsilon\downarrow0\).

Condition on a uniformly sampled level-\(l\) gap under \(\mu_\lambda\). Couple its location to a point drawn from \(\nu_b\) with the same first \(l-1\) digits. The copied prefix is uniform, and this coupling has displacement at most \(b^{-(l-1)}\). Its mean displacement after averaging levels is at most
\[
 \frac{C\epsilon}{Z_\lambda}\sum_{l\ge1}(2/b)^l=O(\epsilon^D).
\]
Thus the fixed first harmonic agrees with \(\widehat\nu_b(1)\) to \(O(\epsilon^D)\).

Now take \(n=b^k\). The union of level-\(k\) cylinders gives, after applying the harmonic, exactly a rescaled copy of the first Gibbs harmonic at \(\lambda'=\lambda n^{-\alpha}\). All prefix phase factors are one. The total Gibbs probability of the remaining early gaps \(l\le k\) is at most \(C2^k\epsilon/Z_\lambda=O((n\epsilon)^D)\). Hence as \(n\epsilon\to0\),
\[
 \widehat\mu_\lambda(b^k)=\widehat\nu_b(1)+O((n\epsilon)^D).
\]
The limiting Bernoulli coefficient is
\[
 \widehat\nu_b(1)=-\prod_{j\ge1}\cos\big(\pi(b-1)b^{-j}\big)>0.
\]
The first factor is strictly negative, every later factor strictly positive, and the infinite product of the latter factors is nonzero by the summability of their deviations from one. For \(0<\alpha<2\) and \(\lambda=2\pi^2 c n^2\), this proves nonuniformity of the circular output.

The same nonuniformity applies to a continuous variance \(v\) with \(v-v_0\sim c\,d^\alpha\) uniformly as \(d\downarrow0\), and \(v>v_0\) off \(C_b\). On a small distance neighborhood the difference of the two unnormalized kernels is bounded by \(C\delta z e^{-z/2}\), where \(z=2\pi^2c n^2d^\alpha\) and \(\delta\) is the uniform relative deviation there. Its integral is \(O(\delta Z_\lambda)\) by the mass estimate and the inequality \(z e^{-z/2}\le C e^{-z/4}\). Away from this neighborhood, both integrals are exponentially small. First send \(n\to\infty\) and then \(\delta\to0\). The normalized laws differ by \(o(1)\) in total variation, uniformly for testing any harmonic.

## Exact critical renewal formula

Now let \(\alpha=2\), \(\lambda=a n^2\), with \(a>0\) fixed. Denote \(I_a(T)=\int_0^T e^{-ay^2}\cos(2\pi y)\,dy\) and \(J_a(T)=\int_0^T e^{-ay^2}\,dy\). Write \(G_a=\int_0^1 e^{2\pi i x}e^{-a d(x)^2}\,dx\). Inside the level-\(k\) cylinders the same scaling gives exactly \((2/b)^k G_a\) and \((2/b)^k Z_a\). In an earlier gap \(l\le k\), both endpoints have phase one at \(n=b^k\); the two edge integrals then add as a cosine. Dividing numerator and denominator by \((2/b)^k\) gives
\[
 \widehat\mu_{a n^2}(b^k)=
 \frac{G_a+\sum_{s=0}^{k-1}2^{-s} I_a(Hb^s)}
      {Z_a+\sum_{s=0}^{k-1}2^{-s} J_a(Hb^s)} .
\]
The infinite sums are absolutely convergent and the truncation error is \(O_a(2^{-k})\). This limit depends on the resolved gaps and does not generally coincide with a normalized Gaussian convolution of the Bernoulli measure.

Let \(\mathcal N_b(a)\) denote the infinite numerator. At small \(a>0\),
\[
 I_a(Hb^s)
 =-\frac{2a}{(2\pi)^2} H b^s e^{-aH^2b^{2s}}\cos(2\pi H b^s)
   +O(a^{3/2}),
\]
uniformly in \(s\). Integrate by parts four times, using that \(\sin(2\pi H b^s)=0\), the odd Gaussian derivatives at zero vanish, and \(\|f^{(3)}\|_\infty+\|f^{(4)}\|_1=O(a^{3/2})\) for \(f(y)=e^{-ay^2}\). The cosine has the same sign for every \(s\), positive for even \(b\) and negative for odd \(b\). Also \(G_a=O(a)\) after subtracting the constant Fourier integral, and
\[
 a\sum_{s\ge0}2^{-s}H b^s e^{-aH^2b^{2s}}\asymp a^{(1+D)/2}.
\]
For the last comparison, split the sum where \(H b^s\asymp a^{-1/2}\), retain one index there for the lower bound, and use geometric growth before that point and Gaussian decay after it for the upper bound. Since \(D<1\), \(\mathcal N_b(a)<0\) for even \(b\) and \(>0\) for odd \(b\) at small \(a\). At large \(a\) the endpoint estimate makes the normalized critical limit positive. Continuity in \(a>0\) follows from uniform absolute convergence locally there. Thus for even \(b\) the critical leading coefficient has at least one zero.

For the exact distance-square profile, any tuned zero still has nonzero finite resonant coefficients. If \(\mathcal N_b(a)=0\), its finite numerator is precisely the negative omitted tail, with
\[
 G_a+\sum_{s=0}^{k-1}2^{-s}I_a(Hb^s)
 \sim -2^{1-k}I_a(\infty)
 = -2^{-k}\frac{\sqrt\pi}{\sqrt a}e^{-\pi^2/a}.
\]
Consequently an exact distance-square field has nonuniform output for every \(a>0\), including the tuned values. Uniform \(o(d^2)\) errors preserve any nonzero leading critical coefficient by the normalized total-variation argument. At a tuned zero that argument provides no control at the smaller \(2^{-k}\) scale.

## Superquadratic powers in the exact distance family

The same finite renewal identity holds for arbitrary \(\alpha>0\) after replacing \(e^{-ay^2}\) by \(e^{-ay^\alpha}\) and allowing \(a=\lambda/n^\alpha\). For \(\alpha>2\) and the heat scale, \(a=a_k=(2\pi^2c)b^{-(\alpha-2)k}\to0\). Three integrations by parts, with \(f(y)=e^{-ay^\alpha}\) and \(\omega=2\pi\), give uniformly for \(T_s=Hb^s\),
\[
 I_{a,\alpha}(T_s)=
 -\frac{a\alpha}{\omega^2}T_s^{\alpha-1}
          e^{-aT_s^\alpha}\cos(\omega T_s)
 +O(a^{2/\alpha}).
\]
Here \(f'(0)=0\), the endpoint sine factors vanish, and \(f^{(3)}\) is integrable for \(\alpha>2\) with \(\|f^{(3)}\|_1=O(a^{2/\alpha})\) by scaling.
Now \(a\sum_{s=0}^{k-1}2^{-s}T_s^{\alpha-1}e^{-aT_s^\alpha}\asymp a^{(1+D)/\alpha}\). The same split at \(T_s\asymp a^{-1/\alpha}\) gives the comparison because the maximizing-scale index is \(s=(1-2/\alpha)k+O(1)\) and lies strictly within the retained sum for large \(k\). The endpoint terms all have the same parity sign. Their sum dominates \(G_a=O(a)\) and \(O(a^{2/\alpha})\) because \(D<1\). Finally the denominator is comparable to \(a^{-(1-D)/\alpha}\) by the corresponding split of
\(\sum_{s=0}^{k-1}2^{-s}\min(T_s,a^{-1/\alpha})\). Thus, for the exact distance-power profile,
\[
 |\widehat\mu_{a_kn^\alpha}(b^k)|\asymp a_k^{2/\alpha}
\]
with eventually negative sign for even \(b\) and positive sign for odd \(b\). This also proves nonuniformity for \(\alpha>2\). Once the normalized leading term goes to zero, the total-variation argument for arbitrary relative perturbations does not preserve its sign or nonvanishing.

## Quantitative transfer to profiles comparable at a specified rate

If the uniform relative variance error near \(C_b\) is \(O(d(x)^\beta)\), \(\beta>0\), the same comparison of positive kernels gives normalized total-variation error \(O(\epsilon^\beta)\). Indeed the pointwise error in the distance tube is bounded by \(C\,d(x)^\beta z e^{-z/2}\) with \(z=(d(x)/\epsilon)^\alpha\). Since \(d(x)^\beta=\epsilon^\beta z^{\beta/\alpha}\), its integral is \(O(\epsilon^\beta Z_\lambda)\) by absorption of the polynomial weights into \(e^{-z/4}\). Uniform \(o(d(x)^\beta)\) analogously gives \(o(\epsilon^\beta)\). Consequently, at a tuned quadratic critical zero the nonzero resonant residual transfers under relative error \(o(d(x)^D)\), as \(\epsilon^D\asymp2^{-k}\). For \(\alpha>2\) it transfers under relative error \(o(d(x)^{\alpha-2})\), as \(\epsilon^{\alpha-2}\asymp a_k^{2/\alpha}\).
