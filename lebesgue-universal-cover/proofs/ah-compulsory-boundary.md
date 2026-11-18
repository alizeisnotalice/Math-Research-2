# Compulsory pieces of the AH boundary

Fix the Gibbs angle $1.5497^\circ$, write $D=T_a-Z$, and retain the original normalization into the doubly truncated regular hexagon $P$. Every closed universal cover $C_{\rm cov}\subset P\setminus A_H^\circ$ contains the straight segment
$$
 \{Z+wD:0.429\le w\le1\}
$$
and the whole boundary chain from $q(0.429)$ to $P_0$ through $T'=q(1)$. Here $q$ is the AH pairing of the straight segment and $P_0$ is the endpoint of the retained $G$-centered arc. The subscript distinguishes this point from the polygon $P$.

For a closed universal cover contained in the original Gibbs body, a further disjunction holds on the $C$-centered arc $T_aY$. At each point $x$ of this arc the cover contains $x$ or its reflection in the vertical axis. In particular it contains $Y$.

## Width-one witnesses and normalization

Use the chart $X(\ell,t)$ and pairing $(\psi(\ell),-t)$ of the AH proof. At a fixed ideal parameter pair, put $x=X(\ell,t)$, $x'=X(\psi(\ell),-t)$, and use the same five cyclic anchors as in the single-point obstruction,
$$
 L(\ell),\quad Q_A L(\psi(\ell)),\quad Q_Ax',\quad N(t),\quad x.
$$
Consecutive distances are one by the exact DC pairing and the two unit-circle equations. On both closed parameter ranges checked below every other distance is strictly less than one. These properties persist under small perturbations of $\ell$ and $t$. Completing the anchors to a width-one body determines its support on five vertex fans and their antipodal fans. On a fan it has the formula $u\cdot v+c$, with $c=0$ for a vertex and $c=1$ for an antipodal arc center.

The three original outward normals $n_3,n_4,n_5$ support, respectively, the EF, DE and DC anchors. The fan signs below establish this throughout both families. The anchor positions on the original faces and width one therefore give $h(n_j)=1/2$ for all six normals. Thus the compatibility function
$$
 F(\phi)=h(u(\phi))-h(u(\phi+60^\circ))
                    +h(u(\phi+120^\circ))-1/2
$$
vanishes exactly at each multiple of $60^\circ$, provided its value is strictly negative on the open interval $(0,60^\circ)$. The width relation guarantees $F(\phi+60^\circ)=-F(\phi)$.

On each parameter rectangle the certificate tests strict negativity of the angular derivative near $0$, strict positivity near $60^\circ$, and strict negativity of $F$ in the middle. The middle test uses the signed derivatives in the parameter and the angle together with the function at the center of the rectangle. Across a bounding ray between two fan formulas, their marked center difference is parallel to the ray, so their angular derivatives agree. The derivatives in the parameter agree as well, by differentiating equality of the two formulas along the moving ray. The enclosure therefore applies also to rectangles including fan edges. All strict sign properties persist for small perturbations of $\ell,t$ at any fixed interior ideal parameter. The exact equalities at angle zero remain true because the three supporting anchors remain on the same faces.

Any isometric placement into $P$ must fill the six original unit strips. Eliminating translation gives a zero of $F$, with a reversal of the three terms for orientation-reversing maps. Only the twelve dihedral orientations are possible, and all translations are zero since the unshifted original support values are $1/2$. The certificate tests the genuine slanted supports $m_3,m_5$ for the orientations that require numerical exclusion.

## The straight and paired tails

For $0.429\le w\le1$, take $(\ell,t)$ from $x=Z+wD$, with $t=g(\ell)$. The uniform angle and vertex-sector tests are strict throughout this closed $w$ interval. Nine of the twelve orientations violate a true slanted support with uniform positive margin. The analytic exclusion of the vertical reflection uses $L(\ell)>C$ along the DC direction $e_L$ when $w<1$. Indeed $T_a-Z$ is tangent to the unit $C$-circle, so
$$
 |Z+wD-C|^2=1+(1-w)^2|D|^2 .
$$
The projection $(Z+wD-C)\cdot e_L$ is positive by its two endpoint values. Solving for the smaller unit-circle root on DC therefore puts $L(\ell)$ strictly beyond $C$. The reflected $m_3$ normal is $m'_5$, with $m'_5\cdot C=1/2$ and $m'_5\cdot e_L>0$. Thus the vertical reflection violates $m_3$, also for small perturbations of $t$ at fixed $w<1$.

Only the identity and $Q_A$ remain. To show that the identity anchor $x$ must be retained, use $t=g(\ell)-\epsilon$. Its paired point then has $t'=g(\psi(\ell))+\epsilon$, and $Q_A$ sends the copy into old AH. The identity anchor tends to $x$, so omission of $x$ by a closed cover excludes that placement for sufficiently small $\epsilon$. To force retention of $q(w)$, instead use $t=g(\ell)+\epsilon$. The identity now enters old AH, while the paired anchor in the $Q_A$ image tends to $q(w)$. This applies to every $w<1$ in the range. Closedness gives the two endpoints at $w=1$.

## The retained $G$ arc and the point $Y$

Set $\ell=\ell_C$ with $L(\ell_C)=C$. The paired DC center is $L(\psi(\ell_C))=G$. As $t$ runs from $t_T$ to $\tau$, defined by $N(\tau)=E_3$, $x(t)=X(\ell_C,t)$ traces $T_aY$ and $x'(t)=X(\psi(\ell_C),-t)$ traces the retained $G$-arc $T'P_0$. These are the two paired arcs in the original AH chart.

The certificate repeats the compatibility calculation on the closed $t$-range. The pinning for $n_3$ uses one exact endpoint inequality. At the EF vertex the last fan ray is $-Q_A(x'-G)$. Its cross product with $n_3$, in that order with $n_3$ first, is $(x'-G)\cdot e_L\ge0$, by the smaller-root branch on DC. Equality holds at $P_0=G+n_2$. The other bounding ray has strict sign. This supplies the exact value $h(n_3)=1/2$ including the endpoint, while $n_4,n_5$ are strictly pinned.

The unchanged eight other motions fail by the numerical slanted support test. The rotation $R_{300^\circ}$ is tested on the initial portion of the $t$-range. Near its endpoint the certificate pins $h(m_3)=m_3\cdot N(t)$. For $t<\tau$ this is strictly below $1/2=m_3\cdot E_3$, since $m_3\cdot e_N>0$. Consequently $h(m_0)>1/2$, and the rotated copy violates $m_5$, whose inverse image is $m_0$. These strict inequalities persist for small perturbations of $\ell$ at fixed interior $t$.

For $t_T<t<\tau$, perturb to $\ell>\ell_C$. The vertical reflection fails by its strict DC cap violation. The identity now enters old AH since $t>g(\ell_C)$. The paired anchor in the $Q_A$ copy tends to $x'(t)$. Thus every point of the retained $G$-arc is compulsory, by the limits at both endpoints.

For the two-sided disjunction, instead perturb to $\ell<\ell_C$ at the same $t$. The identity and vertical images contain points tending to $x(t)$ and its vertical reflection. The $Q_A$ image now violates the old true $G$-arc support at $x'(t)$. To see this, put $u=x'(t)-G$, the outward unit radial. At the ideal parameter,
$$
 \frac{d}{d\ell}\big[u\cdot X(\psi(\ell),-t)\big]_{ell=\ell_C}
 =\psi'(\ell_C)u\cdot e_L<0.
$$
This follows by differentiating the unit-circle equation. The pairing derivative is strictly negative, and the smaller-root branch gives $u\cdot e_L>0$ on the interior arc. Decreasing $\ell$ strictly raises the projected point past that old support. If a closed cover omits both points of the disjunction, choose the perturbation small enough to fail the remaining two placements as well. At $t=\tau$, $x(\tau)=Y$ is fixed by vertical reflection. Closedness and the disjunction for interior $t$ force retention of $Y$.

## Global necessary parameter walls

Set $b=-D/|D|$, $n=Jb$, $w_0=0.429$, and
$$
 \lambda_0=\frac{b\cdot(q(w_0)-Z)}{|D|},\qquad
 \eta_0=\frac{n\cdot(q(w_0)-Z)}{|D|}.
$$
For any positive-angle cut parametrized by $(a,s)=(p\tan\delta,\tan\delta)$ with $p\ge0$ and $0<\delta<\pi/2$, preserving the two compulsory tails and the $G$-arc is equivalent to
$$
 a\le w_0s,\qquad a+\lambda_0s\le\eta_0 .
$$
Necessity follows by testing the first point of each tail. The straight tail is handled by $p\le w_0$. Along the paired chain the projection $b\cdot(q-Z)$ is strictly increasing. The program checks its derivative on the original paired arc and the tangent projections at both ends of the short $G$-arc. Convexity gives a lower boundary graph in this projection, with height zero at $Z$. For a convex nonnegative height function $f(\lambda)$ with $f(0)=0$, the quotient $f(\lambda)/(p+\lambda)$ increases with $\lambda>0$ for any $p\ge0$. This follows directly from $f(\lambda_1)\le(\lambda_1/\lambda_2)f(\lambda_2)$. Its minimum on the required chain is therefore at the first point. The numerical values are enclosed at $lambda_0=0.501057141818249ldots$, $eta_0=0.017440956357112ldots$.

In $(a,s)$ coordinates the two necessary walls and the two nonnegativity constraints define a triangle. The largest possible removal profile among all its affine lines is the maximum of zero and the two affine functions from the nonzero vertices. This is the exact global domination bound from the compulsory chain. The more restrictive sufficient profile of the certified double-line cover is stated separately.

## Computation

Run `python -B certify_ah_compulsory_tail.py` in `computation` with Python 3.13.18 and python-flint 0.9.0. It emits `ah-compulsory-tail.json`, with separate receipts for the (w)-tail and the true (G)-arc tests.
