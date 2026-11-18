# An endpoint obstruction strictly after the mean time

Let a bounded deterministic Cantelli graph have exact exponential transform, terminal spatial marginal g dx=N(0,1), 0≤f≤H, and mean time 1. Use Q_t for the signed remainder including stops at t. Let I=(a,b), L=b−a, with a neighborhood of [a,b] contained in {f<τ} for some 0<τ<1. For continuous f this holds on every compact subinterval of {f<1} with a suitable τ. Let φ(x)=sin(π(x−a)/L) on I, λ=π²/(2L²), and κ=π/(2L).

On a neighborhood of [τ,H]×[a,b] the signed remainder has a smooth caloric density q. Indeed sources from a spatial neighborhood there occur with a uniform time gap (take a slightly smaller neighborhood), and for all other sources the Gaussian heat kernel and all derivatives extend by zero to nonpositive ages because of the spatial gap. Differentiation is locally dominated at each order. At the terminal time the transform implies Q_H=0 by Fourier uniqueness. The same argument applies to measurable f under the essential neighborhood bound. We assume H>1, since H≤1 and the mean condition already imply f=1 almost everywhere.

Write W(t,x)=∫_t^H q(r,x)dr on the cylinder. Then W(H)=0, W_t=−q=ΔW/2 there. The global occupation identity ∫₀ᴴ Q_r dr=∫₀¹ p_r dr and Q_r≤p_r as measures (see Round 3) imply for τ≤t≤1

$$
W(t,x)\ge \int_t^1 p_r(x)\,dr.
$$

The measure inequality is pointwise on the cylinder by continuity. For M(t)=∫_I W(t,x)φ(x)dx, twice integrating by parts gives M'+λM=κ[W(t,a)+W(t,b)]. Integrate from τ to H, using M(H)=0. The endpoint terms up to time 1 are nonnegative. With A=∫_τ¹∫_I p_r(x)φ(x)dx dr>0 this yields

$$
\int_1^H e^{\lambda(t-\tau)}[W^-(t,a)+W^-(t,b)]\,dt\ge A/\kappa.
$$

As W^-(t,x)≤∫_t^H q^-(r,x)dr, changing integration order gives the explicit late-time obstruction

$$
\int_1^H [e^{\lambda(r-\tau)}-e^{\lambda(1-\tau)}](q^-(r,a)+q^-(r,b))\,dr\ge (\pi/L)A.
$$

Here u^-=max(0,−u). Thus in each connected component of {f<1}, at most one spatial point can satisfy q(t,x)≥0 for all 1<t<H. The same at-most-one conclusion holds for W, a stronger tail-in-time test. This makes no assumption on the number or shape of other components, the behavior at infinity, or the sign of Q before time 1.

For a spatially integrated version translate both endpoints by h∈[0,d], assuming the whole translated family stays in the same early neighborhood at the same τ. Write A_h for the resulting initial lower bound. Adding the two endpoint strips, with multiplicity if overlapping, their spacetime integral of q^- after 1 is at least (π/L)∫₀ᵈ A_h dh/[exp(λ(H−τ))−exp(λ(1−τ))]. Division by H−1 gives a lower bound for their supremal single-time negative mass. If d≤L the strips have disjoint interiors. This is a necessary condition on any putative continuous counterexample. The inequalities permit a signed remainder with sufficiently large late negative part.
