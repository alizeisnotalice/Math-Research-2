# A one-wing logarithmic improvement of the remote-window bound

Let w>1, v=1−1/w and suppose that on the favored wing (say +∞) f≥1+η eventually, η>0. No assumptions are made on the other wing, on intermediate zero points or on the escape of the tilted mass. Then

    liminf_{x→+∞} [b_w(x)−x²/(2w)−v log x+(v/2)log(2log x)] ≥ 0.

Write t=x/√w, ρ=1/√w and b₀=t²/2. The fractional envelope gives Iρ(t)=∫exp(ρty−t²f/2)dγ≤(ρ^(−ρ²)+o(1))(√(2log t)/t)^v, using base s=t/ρ=x and k=ρ². For a fixed ε>0 test b=b₀+d, d=v log x−(v/2)log(2log x)−ε. In |y|≤δ the contribution to the window is bounded by exp(d+vδ²/2) Iρ(t)/√w, since d is positive eventually and exp(−df)≤1. Its limsup is at most exp(−ε+vδ²/2) by cancellation of the constants.

For δ≤y≤R the contribution is exponentially small in t after division by exp(d), by comparing Iρ with R_t: the extra factor is exp(−(1−ρ)ty) times the bounded exp(vy²/2). For y≤−δ use the original Lebesgue integral with factor exp(ρty−y²/(2w)−(t²/2+d)f), giving an exponentially small bound. On y≥R the original normalized window is ≤exp(−ηb). Taking a small δ makes the total strictly below one and hence b_w(x)>b. Let ε↓0. This lower bound is saturated under the stronger early canonical hypotheses and the separated right tail: local saddle asymptotics supply the matching lower contribution in the window and hence the matching upper bound on b_w. Reflection gives the left-wing version.

The separation condition can be weakened: it suffices that for some R>0 the integral over y>R of exp(ρty−t² f(y)/2−y²/(2w))dy is exponentially small on the t² scale, since exp(d)=exp(O(log t)) and d≥0. This formulation allows sufficiently thin separated dips in the favored tail.
