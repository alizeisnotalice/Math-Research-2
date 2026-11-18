# Initial heat-ball comparison and the remote critical coupling

Let g=p₁, and let the bounded continuous normalized Cantelli graph f be exact. Assume the early positive first-exit band of R084/R088 on some (0,T), and f≥1+η on both sufficiently far wings for some η>0. The early band includes a unique zero at the origin and strictly expanding continuous symmetric radii b(s). The argument below only requires the early band for the local expansion; the gap on the far wings enters for the final integral and remote critical coupling.

## Initial profile with absolute error in the exponent

As s decreases to zero and y tends to zero,

    b(s)²/s = log(1/s) + O(1/log(1/s)),
    y²/f(y) = 2L_y + log(2L_y) + O(log L_y/L_y),    L_y=log(1/|y|).   (1)

Take a positive constant c near g(0) and real time u such that
v_c(u)=log(1/(2πc²u))>1. The exact heat ball has survivor
q_c(u,x)=(p_u(x)−c)_+ and radius b_c(u)=√(u v_c(u)). Its initial law is δ₀: the removed density min(p_u,c) has mass tending to zero. The outgoing flux is c b_c/(2u) per wing per unit time, by the distributional heat equation for (p_u−c)_+ at a simple root. Since b_c'=(v_c−1)/(2b_c), the associated stopping density per unit distance inside the swept collar is

    k_c(x)=c v_c(u)/(v_c(u)−1),    |x|=b_c(u).

Put U(u,x)=∫₀ᵘ q(r,x)dr for the actual early survival density, and similarly U_c. The R088 distributional equation, integrated once in time, gives

    (∂_u−½∂²_x)U = δ₀ − g 1_{|x|<b(u)}.

The same equation for U_c has k_c, b_c on the right. Apart from the two boundary curves the indicator is 1_{U>0}: outside the nested interval U=0, and for each interior point at positive time there is an earlier interval strictly containing it with strictly positive Brownian survival density there. Positivity follows, for example, from the positive surviving mass at a sufficiently early time, the positive killed Brownian transition density on a compact interior interval during a later time segment, and continuity and strict expansion; an interior starting subinterval of positive survivor mass can be chosen by survival up to that earlier time. The heat ball has the same property directly.

If two swept densities k₁≥k₂ on the small spatial region swept by either profile, then U₁≤U₂. Indeed on {U₁−U₂>0} the source of their difference is −k₁+k₂1_{U₂>0}≤0. The difference starts at zero and solves the heat equation with bounded source and compact support. Testing by its positive part proves the comparison (or use the parabolic maximum principle by mollification). Hence their radii are ordered b₁≤b₂. This argument is local in the final time; the common δ₀ cancels.

For a fixed small final time s set m=log(1/s). Take c_high=g(0), so k_high≥g(0)≥g, and c_low=g(0)(1−4/m). For all 0<u≤s along the low heat-ball boundary, v_low(u)≥m, and k_low≤g(0)(1−4/m)(1+1/(m−1)). On the ball |x|≤b_low(s)=O(√(sm)), g(x)≥g(0)(1−O(sm))>k_low for small s. To apply the lower-density comparison when a point lies beyond this ball, the comparison only needs g≥k_low wherever both profiles are active (on a point where only the actual profile is active the source difference is −g≤0). Equivalently extend the low swept density by zero outside its final ball. For the high-density comparison extend k_high by g(0) if needed. Therefore b_high(s)≤b(s)≤b_low(s). Their squared radii divided by s are m and m−2log(1−4/m). This proves the first assertion in (1), including the constant zero specific to g=p₁. Substitution s=f(y) gives v=2L_y+log v+O(1/v), v=y²/f(y), proving the second from the first.

## Incomplete tilt

For t→∞ and fixed 0<ρ≤1, define

    I_ρ(t)=∫ exp(ρty−t² f(y)/2)g(y)dy.

The asymptotic is

    I_ρ(t) = (1+o(1)) ρ^(−ρ²) (2log t)^((1−ρ²)/2) t^(ρ²−1).          (2)

In particular the constant at ρ=1 agrees with exactness. The same computation applies with the favored wing reversed. All estimates are uniform for ρ in a compact subinterval of (0,1].

Here are the localization details. Write L=log t, u=ty, v(y)=y²/f(y). Outside a fixed small origin interval the compact part is exponentially small in t² by strict positivity of f there, and the remote integral is bounded by exp(−(1+η−ρ²)t²/2). The negative origin half and 0≤u≤1 contribute O(1/t). For 1≤u and u/t in the fixed small positive interval, (1) implies v(u/t)≤2L+C log L (choose the initial interval so the asymptotic bound increases with L_y). Completing the square with this upper bound shows that outside |u−2ρL|≤√L log L the Gaussian integral in u is bounded by exp(ρ²L−c(log L)²) times a polynomial in L. This is negligible. In the window, (1) gives uniformly

    v(u/t)=V+o(1),   V=2L−log(2L)−2log ρ.

Since u/v=O(1), replacing v by V in ρu−u²/(2v) gives an o(1) absolute error. On the same window g(u/t)=(1+o(1))/√(2π). Integrating the resulting Gaussian over u gives (1+o(1))t⁻¹√V exp(ρ² V/2), proving (2). The Gaussian saddle lies at u=ρV in the equivalent quadratic exponent; a unique stationary point of the actual integrand is unnecessary.

## Remote-center root

Fix w>1 and let b_w(x)>0 denote the nonzero root P_w exp[−b_w(x)(f−1)](x)=1 of R115. Then

    b_w(x) = x²/(2w) + (1−1/w)[log|x|−½log(2log|x|)] + o(1)          (3)

at either end. In particular the constant following the separate log and log-log terms is −(1−1/w)log(2)/2. To see this take x>0, t=x/√w, ρ=1/√w, b=t²/2+d. Then the test integral is

    exp(d) w^(−1/2) ∫ exp(ρty−t²f/2) exp[(1−1/w)y²/2−d f(y)]dγ(y).

For |d|≤C log t the second exponential is 1+o(1) uniformly on the saddle window. The near-origin complement still has the stated negligible bound: for d≥0 the multiplier involving f is at most one (this is all that is needed to bracket the root); if d<0 one can instead replace t² by t²−2C log t in the localization Gaussian. On fixed compact sets with f bounded below the penalty remains of order t². On the far tail f≥1+η, and the Gaussian with variance w now has complete-square gain wρ²t²/2=t²/2; hence this part remains O_w(exp(−ηt²/2+C H log t)). Thus the logarithm of the test equals

    d −½log w + (ρ²−1)log t + ½(1−ρ²)log(2log t) −ρ²log ρ +o(1).

Set d to the cancelling expression plus or minus any fixed ε. The two roots are bracketed in the sense that the test values tend to exp(±ε) respectively, and continuity and uniqueness give the stated o(1) asymptotics of the positive root. On substituting t=|x|/√w all fixed log w constants cancel. This result is conditional on an early first-exit band. Tail gaps alone do not establish it for arbitrary continuous graphs.
