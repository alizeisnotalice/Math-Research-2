# A level-by-level distance budget from a single early cut

Use the notation and reverse graph hitting representation of [round 199](199-graph-trace.md). Let y be any point with f(y)<1. The component I of {E>f} containing y has constant roof σ>1 and W(σ,y)=0. The conclusions also hold with σ=H and I the whole line. No finite maximum, derivatives of f, finiteness of branches or tail conditions are required.

Take 0<δ<1 with f≤1−δ on J=[y−L,y+L], L>0, and put ε=σ−1. Write m= min_{z∈J}∫_{1−δ}¹p_s(z)ds>0. For 1<s<σ let A_s^± be the portions of {f≥s}∩I to the right and left of y. Put d_±(s)=dist(y,A_s^±), c_±(s)=sup_{z∈A_s^±} p_s(z), taking distance infinity and supremum zero for an empty set. Let S(L,T) be the probability that a standard Brownian motion from 0 stays in (−L,L) for time T. Then exactness requires
```
m S(L,ε+δ) ≤ ∫₁^σ Σ_± c_±(s) 2 Φ̄(d_±(s)/√(σ−s)) ds.                 (1)
```
In particular, writing c=(2π)^(-1/2), one obtains
```
m S(L,ε+δ) ≤ 4c ε Φ̄(L/√ε).                                           (2)
```

Indeed start reverse Brownian motion at (σ,y), run to the graph or the deterministic time height 1−δ, whichever comes first. The stopped expectation of W is zero. This follows by bounded caloric optional stopping as in round 199; if σ=H the solution continues as zero above H, and the initial point lies strictly above its graph. Staying inside J throughout ε+δ yields terminal W≥m from formula (3) there, even if the path touches the graph at the final instant. Every negative terminal value must occur at a graph hit Z with h=f(Z)>1. On this event w^−(Z)≤∫₁^h p_s(Z)ds, by dropping the nonnegative earlier-source term. Hits stay in I: each finite endpoint of I has height σ and reaching that endpoint at positive reverse time first forces a graph hit. On {Z∈A_s^±,h≥s} the one-sided Brownian maximum before time σ−s must reach d_±(s), so the reflection principle, Tonelli and the two payoffs prove (1). Using d_±≥L, c_±≤c yields (2).

**Fixed-window small-amplitude obstruction.** For a sequence of exact graphs with a common bounded spatial observation interval, suppose the preceding construction holds with L bounded below by a positive constant and ε→0. If δ→0 as well, then S(L,ε+δ)→1 and m≥c_J δ with c_J>0 uniform (choose a common bounded interval J of fixed radius). Mills' inequality now gives
```
δ ≤ C_J ε^(3/2) L^(−1) exp(−L²/(2ε)).
```
In particular a fixed-power relative valley depth there is impossible. With no small-amplitude limit, (1) weights the distance to *every* high level by the time remaining before that level can be hit. This condition constrains both tails even if high levels escape to infinity. It allows exponentially shallow valleys and does not by itself prove constancy.
