# An output-distance certificate from finitely many spatial sign windows

The input f is any bounded continuous time mark 0≤f≤H; no transform identity is assumed. Fix T>max(H,1), let g=p₁ and denote by Y=X+sqrt(H−f(X))G the output, with independent standard normal X,G. Use the computable function
```
U(t,y)=∫₀¹p_s(y)ds−∫₀ᵗp_s(y)ds+∫g(x)∫₀^((t−f(x))₊)p_s(y−x)ds dx,
w(y)=U(f(y),y).
```
It is bounded and continuous up to time T by the integrated Gaussian estimates in round 199. It is caloric above the graph, including for non-exact inputs. If K_y is the first graph hit from (T,y) under reverse Brownian motion, optional stopping gives U(T,y)=∫w dK_y.

**Certificate.** Partition the line into m ordered intervals B_j and fix signs e_j in {−1,+1}. Suppose that in B_j a closed interval centered at y_j with radius r_j>0 has e_j w≥δ_j>0. Write a_j for the infimum of f on this closed window. Let M>0 bound |w| and let η≥0 bound (e_j w)_− throughout each B_j. Then
```
W₁(law(Y),N(0,H)) ≥ [ Π_j δ_j S(r_j,T−a_j) ] / [m M^(m−1)] − η,
```
where S(r,t)=P(sup_{u≤t}|B_u|<r). The bound holds also when the right side is nonpositive. The constants and sign assumptions can be checked on the time graph itself and do not presuppose a lower bound for the output defect.

**Proof.** Put A_ij=∫_{B_j}(e_jw)_+ dK_{y_i}. The noncrossing determinant formula of round 202 applies. A sufficient event for a matching without collisions before all hits is that each path remains in its own window until it stops. Requiring its underlying Brownian path to stay there until T−a_j suffices, since by then its time height has reached the infimum on that window. After graph hits take the spatial paths to be frozen. By independence and positivity of the endpoint weights,
```
det A ≥ Π_j δ_j S(r_j,T−a_j).
```
The row sums of A are at most M, so Hadamard bounds each cofactor by M^(m−1) (the empty minor equals one). The error in replacing ∫w dK_{y_i} by (Ae)_i is at most η. Applying the adjugate to Ae and using |e_j|=1 proves
```
sup_y |U(T,y)| ≥ det A/(m M^(m−1))−η.
```
Finally, in the sense of spatial distributions, U(T,·) has second derivative 2(law(Y_T)−N(0,T)), where Y_T=Y+sqrt(T−H)G' with an independent Gaussian G'. This follows directly from ∂²p_s=2∂_s p_s and integration at zero age. The difference E|Y_T−y|−E|N(0,T)−y| has exactly the same second derivative. Both functions tend to zero at the two infinities: for U use dominated convergence on the finite age integrals, and for the potentials use finite first moments and equal means zero. Thus the functions agree. Since the absolute distance to y is 1-Lipschitz and convolution with a common independent Gaussian cannot increase W₁, this completes the proof.

**Checking a family from one graph estimate.** If ||f−f₀||∞≤b then ||w_f−w_{f₀}||∞≤sqrt(2/π)(1+sqrt(2))sqrt(b): the single upper age in the negative free term moves by at most b, and the upper age in the source term by at most 2b. Use the bound for integration of a Gaussian over any age interval. Accordingly positive margins and wrong-sign tolerances on a reference graph extend with explicit loss to a whole sup-norm neighborhood, with a_j reduced by b if needed. The sign assumptions over the two unbounded blocks may be verified with tail bounds.

**Independent finite calibration.** The accompanying standard-library script computes an ordered first-hit matrix on a nonflat time barrier for three nearest-neighbor walks with common parity before stopping. It separately propagates the killed three-path process. The two resulting integers agree exactly at 652288, giving the common probability 637/16384. This verifies the determinant orientation in a finite model; the Brownian certificate follows from the proof above.
