# Spatial alternation of the heat density on every delayed copy of a bounded roof valley

Suppose I=(a,b) is a bounded component of {E>f} for a continuous exact graph. Its roof is constant σ, with f(a)=f(b)=σ. Extend Q by zero beyond H. For every
```
0<τ<σ−inf_I f
```
the continuous function
```
u_τ(z)=q(f(z)+τ,z),    z∈I,
```
has arbitrarily many strict spatial sign alternations. Its support is contained in the compact set {z∈[a,b]: f(z)≤σ−τ} inside I. In particular a sequence of distinct zeros and alternating nonzero values accumulates inside I. For τ increasing to σ−min f, all the oscillations are confined to arbitrarily small neighborhoods of the minimum set of f on [a,b].

Indeed q is smooth on the open region {(s,z): s>f(z), s>0}, by the heat equation with graph source; the distribution extended past H is still caloric on this open region. It vanishes for z∈I, s≥σ by local extinction and continuity when s>f(z). Define U(t,z)=q(t+τ,z) for z∈[a,b] and f(z)≤t≤σ. This whole compact set, after the positive time shift, is inside the open smooth region. U is therefore bounded and continuous to its boundary, and caloric in the interior. For z∈I and f(z)<t≤σ the stopped reverse Brownian representation gives U(t,z)=∫u_τ dK_{t,z}. No lateral exit is possible before first graph hit because f(a)=f(b)=σ. At t=σ this annihilates u_τ for every start in I. The proof of ordered strict positivity in [round 202](202-graph-trace-alternation.md) applies to this same hitting kernel.

If u_τ vanished identically, U would vanish throughout the region just described. Choose a small nonempty open subinterval J compactly contained in I with sup_J f<σ−τ. We would have q(s,z)=0 on J for s in an open neighborhood of σ: below σ the representation applies whenever s−τ>sup_J f, and at and above σ there is extinction. Thus the spacetime distribution Q has no support in that neighborhood of σ over J. It also vanishes over I for all larger times, contradicting E=σ on J by the definition as the upper support time. Hence u_τ is nonzero, and the sign-block determinant argument forces all the asserted alternations. Compact containment follows from continuity at the two spatial endpoints. All accumulation statements allow shrinking amplitudes; they refer to the heat density along a delayed graph, not to the sign of f−1.
