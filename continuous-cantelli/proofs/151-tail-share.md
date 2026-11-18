# A monotone far tail and its actual share of the edge flux

Assume an exact graph with canonical core [−b,b] at time T>0. Write v=f−T on the exterior, and suppose v(x) is eventually nondecreasing to 1 on the right. Arbitrarily many other exterior wells and arbitrary continuous edge slopes are allowed. Let A(t)=∫_{|x|>b}exp(t(x−b)−t²v/2)dγ∼J/t² and let π_t be its normalized probability law. Put p(t)=π_t{X>t/2}, and

    D(x)=x²(v(x)−1)+2bx−4 log x.

Then exp(−D(t)/2)=J exp(−b²/2)p(t)+o(1), and sup_{|s|≤K}|p(t+s)−p(t)|→0 for every fixed K. Thus a positive limiting fraction determines the constant term in the critical tail expansion; a vanishing fraction forces D→+∞. The fraction itself need not be assumed to converge.

**Macroscopic sources.** Put U=X−b−t v. The ratio A(t+s)/A(t)=E_{π_t}exp(sU−s²v/2) tends to one for fixed s. Since v is bounded, this gives uniform exponential moments of U of every fixed order by taking both signs of s. Also X/t−v→0 in probability. The left exterior has exponentially small probability since v≥0. Inside any fixed bounded right interval, v→0 in probability on that interval because v=(X−b−U)/t. Outside arbitrarily large right intervals v is uniformly close to one. Thus any weak limit of (U,v) has v in {0,1} and its mass on v=1 is the subsequential limit p of p(t). All shifted likelihood identities pass to the limit by uniform integrability. The part at v=0 is supported on U=0: its moment generating integral is ≤1 for every real s, excluding positive or negative U by sending s to the corresponding infinity. The remaining U measure at v=1 then has moment generating function p exp(s²/2), so it is p times the standard normal. Reweighting by exp(sU−s²v/2) therefore preserves the mass p on this component. The distinction between thresholds (t+s)/2 and t/2 vanishes in π_t probability since X/t is close to {0,1}; uniform integrability transfers the claim to shifted tilts. This argument along subsequences also works for any bounded s tending to a limit, proving the stated uniform shift estimate.

**Reading the far graph from its flux.** Complete the Gaussian square and restrict t²A(t) to x∈[t−1,t]. Monotonicity v(x)≤v(t) lower bounds this portion by a fixed positive constant times exp(−D(t)/2), so W=exp(−D/2) is bounded on the far right. Extend it to a bounded function on the line. On x=t+w>t/2 the normalized true tail integral t²A(t)p(t) has integrand

    g(w) exp[2log t−bt+bt²/x−2(t/x)²log x] W(x)^((t/x)²).

For w in a fixed compact interval the exponential factor tends uniformly to exp(−bw); powers of W differ from W by o(1) uniformly there since W is bounded and the powers tend uniformly to one. On x>t/2 the entire integrand is bounded by Cg(w)exp(C|w|): the term bt²/x−bt is bounded above by 2b|w|, and the logarithmic difference vanishes at w=0 and has bounded derivative uniformly for x>t/2, t large. The power of W is uniformly bounded. Consequently

    J p(t) = ∫g(w)exp(−bw) W(t+w)dw+o(1).

On any sequence t_j where p(t_j)→p, extract a local weak-star limit of the bounded translates W(t_j+·). The convolution identity and shift estimate imply that its convolution with the shifted Gaussian is the constant Jp. Injectivity as in round 143 implies the limit is the constant Jp exp(−b²/2), uniquely.

Finally for 0≤u≤1, monotonicity gives D(t+u)≥D(t)+(2tu+u²)(v(t)−1)+2bu−4log(1+u/t)≥D(t)−Cu for large t, using the lower bound for D(t). Hence W(t+u)≤exp(Cu/2)W(t). Any fixed excess of W(t_j) above the weak-star limit persists to the left on a small interval; any fixed deficit persists to the right. Thus W(t_j) itself tends to that constant. Applying this to every subsequence proves the assertion. No rates of regularity of v or assumed higher-order expansion of its tail are used.
