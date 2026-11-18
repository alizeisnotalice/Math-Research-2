# Shared sharp constants at the right edge

Write η_T(y) for the density in round 137, J=−∂_yη_T(b−)>0. For the canonical boundary at a positive time the edge flux result of rounds 88 and 90 gives J=2g(b)b'(T). Thus

    lim_{z→∞} z² ∫ exp(z(y−b))η_T(y)dy = J.                         (1)

This also follows directly from η_T(b−u)=Ju+o(u), the linear edge upper bound and u→u/z. The same limit holds for the corresponding outer graph integral in round 137. Every source there is nonnegative. Contributions from disjoint spatial regions therefore share J.

Consider finitely or countably many distinct right-hand points x_i>b with f(x_i)=T. At any specified one-sided flank suppose, as u↓0,

    v(x_i±u)=f(x_i±u)−T=K_i/[L−β_i log L+c_i+o(1)],
    L=log(1/u),         K_i=(x_i−b)²/2.

Then β_i≥5/2. For β_i=5/2 each such flank consumes at least 2√π g(x_i)exp(c_i)/K_i from (1). Each of the two flanks at a single point is counted separately when it satisfies the hypothesis.

Indeed put d=x_i−b and K=d²/2. In a sufficiently small fixed flank the distribution of v near zero has sublevel mass
M(s)∼A s^β exp(−K/s), A=g(x_i)exp(c_i)K^(−β). This follows by sandwiching the denominator by L−β log L+c±ε, inverting the two monotone comparison functions, and sending ε↓0. This argument does not require monotonicity of the flank itself. Its Laplace integral is asymptotic to

    A√π K^(β/2+1/4) q^(1/4−β/2) exp(−2√(Kq)).

For detail integrate by parts to q∫M(s)exp(−qs)ds, set s=√(K/q)w, and use w+1/w−2=(w−1)²/w. On [1/2,2] the rescaling (Kq)^(1/4)(w−1) yields the Gaussian integral √π by dominated convergence. Outside that interval the excess in the exponent bounds both integrals exponentially relative to the saddle (split w≥2 at 4, use linear decay beyond 4, and invert w for w≤1/2). Heights bounded away from zero are exponentially negligible. Restricting to u≤exp(−√z) does not affect this asymptotic: in the removed part the contact form gives v≥c'/√z, yielding O(exp(−c''z^(3/2))) before spatial tilting. On the retained part exp[z(x−b)]=exp(zd)(1+o(1)). With q=z²/2 the exponent cancels and multiplication by z² yields a constant times z^(5/2−β). This proves the divergent and critical cases, including the constant stated above.

Two optional contributions can be included in the same sum:

* If v(b+u)=k u+o(u) with 0<k<∞, the edge costs at least 2g(b)/k. Restrict to u=w/z² with bounded w, and then let the upper cutoff tend to infinity in ∫exp(−kw/2)dw.
* If f(x)=1+T−2b/x+[4 log x+c∞+o(1)]/x² at +∞, the far tail costs at least exp((b²−c∞)/2). Restrict to x=z+w with bounded w. Completing the Gaussian square, multiplying by z² and passing to the limit gives g(w)exp(−bw−c∞/2); then expand the w window to the full line.

The sum of all applicable costs is at most J. For countably many terms apply (1) to any finite set of separate truncated windows before increasing the set. In particular if the right slope at the core edge matches the canonical inverse slope k=1/b'(T), it consumes all of J. No critical β=5/2 flank at any further right point and no far tail of the displayed saturated type is then possible. Reflection gives a distinct left-edge budget. These results use the canonical core identity and do not presume the outer graph is realized by a stopping time.
