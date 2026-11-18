# Recovering a signed higher-order residual at the actual noise width

Set the contact at the origin, g(u)=p₁(a+u), and suppose the general exact Gaussian output identity for v. Let w be smooth, nonnegative, w(u)~c u^q for a positive c and even integer q≥4. Assume |v−w|≤C|u|^α locally, α>q, with both signs of u allowed. Write h=r^(q/2), r>0; extend w smoothly and boundedly outside this neighborhood. The functions

    B_r(t)=(v−w)(r+ht)/r^α

have a unique weak-star local limit as r↓0. It is constant. It equals zero for nonintegral α; for integral α its value is the coefficient of u^α in Q'(u), divided by g(0), with Q the round 52 primitive polynomial. On the negative side the same result applies after reflection. In particular the average over any fixed positive multiple of the actual noise width has the same leading value. This resolves the equality of scales in the earlier signed interval estimate for a higher-order perturbation of a positive leading term. It does not identify behavior on every interval of length o(h).

**Microscopic constraint.** Local L-infinity bounds and diagonal extraction give a weak-star subsequential limit B on every compact t-interval. Since for any fixed t the factor |r+ht|/r tends uniformly to 1 on compact intervals, this B is globally bounded by the same finite constant.

Test the exact density with added variance h² by its k-th spatial derivative at a+r+hs and multiply by h^k, for k≥1 chosen so (h/r)^k=o(r^(α−q)). The right side is O(h^k). Integrate in t=(u−r)/h. The scaled kernel is p_T^(k)(s−t) with T=1+v(r+ht)/h². On |t|>r/(2h) its absolute value is bounded by C_k |s−t|^(−k−1), uniformly in T, so this whole tail is O((h/r)^k). Constants here and below may depend on fixed s,k.

On |t|≤r/(2h) subtract the smooth baseline with T₀=1+w(r+ht)/h², including the exact factor g(r+ht). Its integral is O((h/r)^k). Indeed Taylor-expand the kernel as a function of the displacement x=ht in g(r+x)p_{1+w(r+x)/h²}^(k)(s−t) at x=0 through order k−1. All fixed-order derivatives of w at r+x satisfy |w^(j)(r+x)|≤C_j r^(q−j) on this range, including j>q. The relevant normalized variances stay in a fixed compact positive interval and ∂_T p_T=p_T''/2. Hence the integral remainder is O((h/r)^k) by polynomial Gaussian domination. Each integrated Taylor term is an integral of t^j p_{1+w(r)/h²}^(k+2l)(s−t), j<k, times a bounded factor times an appropriate power of h/r; it vanishes on the whole line by integration by parts. Restoring its omitted tails costs a Gaussian tail, uniformly in r.

Here T−T₀=r^(α−q) B_r(t). On this growing input window B_r is uniformly bounded. Taylor expansion in T produces an integrable O(r^(2(α−q))) remainder, by the same fixed positive variance bounds. On bounded t-sets T₀→1+c. Passing to the limit after division by r^(α−q), first truncating by a fixed t-window and then using the uniform Gaussian tails, gives

    p_{1+c}^(k+2) * B = 0.

The Fourier transform of a bounded B is tempered. The only real zero of the multiplier is frequency zero. Division against compactly supported frequency tests away from zero shows that the transform has point support; hence B is a polynomial and its boundedness makes it a constant D.

**Identification across subsequences.** Apply the round 52 bilateral primitive estimate to [r,r+Th], for any fixed T>0. After division by Th r^α, its remainder is bounded by C/T with C independent of T in the limit r↓0, since r+Th≤2r eventually. Along the extracted subsequence the weighted normalized residual average tends to g(0)D, by testing local weak-star convergence on [0,T]. The polynomial part of the average tends to the claimed coefficient (or zero): all degrees of Q' below α vanish and higher degrees give zero in the limit. Let T→∞. Thus D is uniquely determined, giving convergence of the full family. This argument combines the microscopic Gaussian smoothing identity with the macroscopic control of interval masses.
