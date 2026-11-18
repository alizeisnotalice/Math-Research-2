# Exponential inverse moments and a squared-gap bound for isolated highest points

Let H=sup f be attained, v=H−f, and suppose the nonempty maximum set A={v=0} is locally finite. Write the points in their spatial order. At each a∈A define
$$
 L_a=\limsup_{x\to a}\frac{|x-a|}{v(x)},\quad
 K_a=\frac{\pi^2}{8 L_a^2},\qquad
 C_a=\sup\{c>0:\ \int_{0<|x-a|<r}v(x)^{-1/2}e^{c/v(x)}dx<\infty
 \text{ for some }r>0\}.
$$
Use C=0 for an empty positive set, K=∞ for L=0, K=0 for L=∞. On an isolated neighborhood changing r to any sufficiently small positive radius does not affect C.

**Local inequality.** For every a in A,
$$
 C_a\ \ge\ \min\left(K_a,\ \inf_{b\in A,\ b\ne a}
        [C_b+|a-b|^2/2]\right).
$$
For ε>0 the exact transform implies ∫p_{v(x)+ε}(y−x)g(x)dx=p_{H+ε}(y), also at complex y by holomorphic domination for fixed ε. Evaluate at y=a+id where d>0 and d²/2 is strictly below the right side. Near a the real integrand is the positive Gaussian prefactor times
$$
 \exp\big((d^2-(x-a)^2)/(2(v(x)+\epsilon))\big)
 \cos\big(d(x-a)/(v(x)+\epsilon)\big).
$$
The cosine is bounded below by a strictly positive number by dL_a<π/2. Outside a small neighborhood of a its integral is uniformly bounded in absolute value for 0<ε≤1. Indeed for |x−a|≥d+δ the nonpositive exponential with a fixed strictly negative numerator controls the inverse square-root prefactor. Within that compact range only finitely many other maxima b require attention. When |b−a|>d use the same negative-numerator bound near b. When |b−a|≤d the strict inequality assumed above allows c_b>max(0,(d²−|a−b|²)/2) below C_b. Shrink the neighborhood so its exponent numerator divided by two is ≤c_b. Off the single point b the kernel is dominated by a constant times v^(-1/2)exp(c_b/v), locally integrable. The remaining compact set has variance bounded below. The target p_{H+ε}(a+id) stays bounded. Hence the positive near-a integral stays bounded, and Fatou on a sufficiently small punctured neighborhood gives integrability for every c<d²/2. Letting d increase proves the inequality, with the zero right side automatic.

**Ordered global bound.** For a,b∈A let D(a,b) be half the sum of the squared successive spacings from a to b, with D(a,a)=0. In each direction where there are infinitely many points define D(a,±∞) by summing the successive squared spacings in that direction, divided by two. In a direction with finitely many points put this escape cost equal to ∞. Then
$$
 C_a\ge\min\left(\inf_{b\in A}[D(a,b)+K_b],\
                    D(a,-\infty),\ D(a,+\infty)\right).
$$
To see this start with finite C_a. If C_a≥K_a stop. Otherwise the local inequality supplies a distinct point b with C_b+|a−b|²/2≤C_a. The infimum is attained when finite because the relevant points lie in a bounded interval and A is locally finite. Repeat while the local cap has not been reached. Capacities C strictly decrease along the resulting chain, and the sum of its half-squared steps is at most C_a. A finite chain ending at b costs at least D(a,b) since every gap between the endpoints must be crossed, and a step across several gaps costs at least the sum of the individual squared-gap costs. An infinite chain leaves every compact interval and must escape along at least one unbounded direction; its cost dominates the sum of the consecutive squared-gap costs in that direction. This proves the assertion.

In particular if L_a=0 at every maximum and both escape costs are infinite, then each C_a=∞. This includes finite maximum sets, uniformly separated sets, and sequences with shrinking gaps but divergent squared-gap sum on every infinite end. At any of these points even a one-sided upper bound v(a+u)≤B/log(1/u) is impossible for fixed B>0: an inverse exponential with c>B then dominates a nonintegrable power. No positive variance gap at spatial infinity is needed. This criterion is conditional on the stated local geometry; accumulations of maxima in compact sets and vanishing variance with other relative rates are not covered.
