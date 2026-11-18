# Universal logarithmic lower bound for every real tilt

Assume the exact bounded graph identity ∫exp(tx−t² f(x)/2)dγ=1, with f≥0, γ=N(0,1), and write P_t for the tilted probability. For either sign t→±∞ define

    L=log|t|, A=t²E_{P_t}f, Z=tX, T=t²f(X), D=Z−T, V=Z−2L.

Then, with absolute constants,

    A+log A ≥ 2L−1−O((log L)^(3/2)/sqrt(L)),
    A ≥ 2L−log L−1−log 2−o(1).                                     (1)

No tail bound, zero-set geometry, pointwise rate or monotonicity of the graph is assumed. The first and second derivatives of the normalization give

    ET=EZ=A, ED=0, E D²=A,  T≥0.

Equivalently, the first derivative gives E(X−tf)=0 and the second gives E(X−tf)²=Ef. The density of Z under P_t and its differential entropy h satisfy

    p(z)=g(z/t)exp[(V(z)+D(z))/2],
    A/2=L−h(Z)+log(2π)/2+E X²/2,                                  (2)

where g is the standard normal density, V(z)=z−2L and D(z)=z−t²f(z/t). The entropy formula also follows from D_KL(P_t||γ)=A/2. All entropies at a fixed nonzero t are finite because f is bounded and the tilt against γ has a linear-Gaussian tail.

Only A≤2L needs a lower estimate for large L. On V=−w<0 split at D=w/2. The contribution of the region D≥w/2 to E V⁻ is at most 2E|D|≤2sqrt(A); on the other region p(z)≤g(0)exp(−w/4). It follows that

    E V⁻ ≤ 2sqrt(A)+16g(0),   E|V|≤4sqrt(A)+32g(0).

The second inequality uses EV=A−2L≤0. Comparison with a Laplace density of scale 1+E|V| gives h(Z)=h(V)≤log[2(1+E|V|)]+1. Thus (2) first gives A≥2L−log(1+A)−O(1), in particular A≈2L as L→∞ with A≤2L.

Consider B=V+D=2log[p(Z)/g(Z/t)]. Write q for the density of the same Laplace comparison centered at 2L on the Z scale. Its maximum is ≤1/2. Since g≤g(0)<1, we have

    E B⁻ ≤ 2E log_+(1/p(Z))
           ≤ 2E log(1/q(Z)) + 2E log_+(q(Z)/p(Z))
           ≤ O(log A).

Indeed ∫_{p<q} p log(q/p)≤1/e because u log(1/u)≤1/e for 0<u<1 and ∫q=1. The first expectation is at most log[2(1+E|V|)]+1. Since EB=A−2L≤0, E|B|≤2EB⁻=O(log A).

Put W=V/sqrt(A) and S=−D/sqrt(A). Then ES²=1 and E|W−S|=O(log A/sqrt(A)). For K≥1 set ψ_K(w)=w² for |w|≤K and 2K|w|−K² outside. It is 2K-Lipschitz and ≤w², so Eψ_K(W)≤1+O(K log A/sqrt(A)). Its Gibbs normalizer N_K=∫exp(−ψ_K/2)dw is between sqrt(2π) and sqrt(2π)+2exp(−K²/2)/K. Nonnegativity of relative entropy against this density gives, at K=sqrt(log A),

    h(W) ≤ log N_K+Eψ_K(W)/2
         ≤ log(2πe)/2+O((log A)^(3/2)/sqrt(A)).                    (3)

Use h(Z)=h(W)+(log A)/2 in (2), discard the nonnegative EX² term and use A≈2L. This proves (1). If A>2L the same two lower bounds hold directly.

## Equality forces a microscopic parabola

Suppose along a sequence of real tilts

    A+log A−2L+1 → 0.                                               (4)

Then A≤2L eventually, and (2) gives h(W)=[log(2πe)+EX²]/2+o(1). By (3), EX²→0 and h(W)→log(2πe)/2. In the proof of (3) this makes the relative entropy of W against the normalized Huber density exp(−ψ_K/2)/N_K tend to zero. Pinsker's inequality yields total-variation convergence to that comparison, and the comparisons themselves converge in total variation to the standard normal since K→∞. Thus W→N(0,1) in total variation and S→N(0,1) weakly. Because f is a graph, its values can be recovered from this entropy saturation:

    t² f((2L+sqrt(A)w)/t)
         −[2L+log A+2sqrt(A)w+w²] → 0                              (5)

in Lebesgue measure on every bounded w-interval. Indeed the density of W is

    p_W(w)=sqrt(A)g((2L+sqrt(A)w)/t)
             exp[L+sqrt(A)w−t² f((2L+sqrt(A)w)/t)/2].

It tends to the strictly positive standard normal density in L¹. Its logarithm therefore converges locally in Lebesgue measure to −log(2π)/2−w²/2. Also (2L+sqrt(A)w)²/t²→0 uniformly on compact w sets. Rearranging gives (5). No existence of a nonconstant exact graph saturating (1) is asserted.

The denominator form at the origin follows without assuming a logarithmic ansatz. With x=(2L+sqrt(A)w)/t and l=log(1/|x|), (5) gives, locally in Lebesgue measure in w,

    x²/(2f(x))−l−(1/2)log l → (1/2)log 2.                         (6)

Indeed z=2L+sqrt(A)w satisfies z²/(2L)=2L+2sqrt(A)w+[A/(2L)]w², where A/(2L)→1. Thus t² f(x)=z²/(2L)+log A+o(1) in the same sense. Since z/(2L)=1+O(1/sqrt(L)) on bounded windows, z²/[2t²f(x)]=L−(log A)/2+o(1). Finally l=L−log(2L)+o(1) and log A=log(2L)+o(1). A possible zero of f may be assigned the extended denominator +∞; such points are negligible in this convergence.

The companion r079_tilt_check.py checks the entropy constant at a single chosen tilt with a Gaussian calibration, separately from the all-tilt theorem. If the assigned distribution of Z is N(A,B), the density identity in (2) defines a quadratic T(z). In the resulting standard Gaussian coordinate w, matching ED=0 and ED²=A requires

    A=2L−log B−1+(A²+B)e^(−2L),
    A=B(1−2A e^(−2L))²+2(1−B e^(−2L))².

The constructed T has minimum A−(1−B e^(−2L))−B(1−A e^(−2L))²/(1−B e^(−2L)). The 80-decimal numerical solutions in the accompanying JSON approach the lower bound constant for L=10,30,100,1000. This benchmark varies with the chosen tilt and is not proposed as a solution of the original all-tilt graph problem.
