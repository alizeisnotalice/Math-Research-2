# Global domination of positive AH hinge cuts

Fix the Gibbs angle $1.5497^\circ$, the old AH boundary, and the EH cut of depth $0.0005521$. Let $B$ be the resulting fixed convex body. Consider all further cuts by halfplanes of the form
$$
 \eta/|D|\ge a-sw,qquad
 w=-b\cdot(x-Z)/|D|,qquad
 \eta=n\cdot(x-Z),
 \quad (a,s)=(p\tan\delta,\tan\delta),
$$
where $D=T_a-Z$, $b=-D/|D|$, $n=Jb$, $p\ge0$ and $0<\delta<\pi/2$. If the intersection of $B$ with any family of these halfplanes is a universal cover, then its removed region is contained in one fixed chord cap. Its area is greater than $0.8440860344$.

## The dominating chord

Let $w_0=0.429$. The compulsory-boundary theorem forces $x_0=Z+w_0D$ and $q_0=q(w_0)$ in every such universal cover. Write the coordinates of $q_0$ as $w=-\lambda_0$, $\eta/|D|=\eta_0$. Every used line therefore satisfies
$$
 a\le w_0s,\qquad a+\lambda_0s\le\eta_0.
$$
The chord $x_0q_0$ has coefficients
$$
 s_* = \frac{\eta_0}{w_0+\lambda_0},\qquad a_*=w_0s_* .
$$

For $-\lambda_0\le w\le w_0$, each allowed affine height is at most the chord height by comparing at both endpoints. For $w\ge w_0$, each line has nonpositive height, while all of $B$ has $\eta\ge0$. For $w\le-\lambda_0$, the old body's lower height graph in the coordinate $lambda=-w$ is convex, starts at height zero at $Z$, and passes through $(\lambda_0,\eta_0)$. Hence its value is at least $(\lambda/\lambda_0)\eta_0$ there. The allowed line has height
$$
 a+s\lambda\le\eta_0+s(\lambda-\lambda_0)
 \le(\lambda/\lambda_0)\eta_0,
$$
where the last step uses $a\ge0$ and the wall at $q_0$. Thus no allowed line removes anything outside the middle coordinate interval. The chord $x_0q_0$ dominates every such cap on $B$, even for an arbitrary family of lines. It realizes the maximum cut allowed by the necessary condition of retaining the compulsory chain.

## Sharp area ceiling from the compulsory chain

The chord cap has oriented boundary $x_0\to Z\$, $q([0,w_0])$, $q_0\to x_0$. Translation by $x_0$ makes the straight contributions vanish, giving its exact area
$$
 \frac12\int_0^{w_0}(q(w)-x_0)\times q'(w)\,dw .
$$
A 64-cell outward rational calculation with tenth-order Taylor remainders encloses this value in
$$
 (0.0000074358439921451580997320,\;
   0.0000074358439921451580997321).
$$
Subtracting its upper endpoint from the lower area endpoint for $B$ gives the stated lower bound for this family of covers. The EH and AH chord caps are disjoint because the entire AH chord cap lies in the upper part of $B$.

Run `python -B certify_ah_compulsory_chord_area.py` in `computation`. The output is `ah-compulsory-chord-area.json`.

## Allowing the EH depth to vary

Hold the same original Gibbs angle, the same lower EH endpoint and the same positive-turn AH family, but vary the EH upper endpoint along its boundary arc on the nested branch of Round 117. The non-universality witness at depth 0.00055222 excludes that depth and all deeper cuts in this branch, even without the extra AH cuts. For any shallower depth the EH cap area is at most its value at 0.00055222. The compulsory AH boundary argument only requires the original containment in the old AH cover and so applies at every one of these depths. Subtract the two upper bounds for removed areas from the original Gibbs area lower endpoint. This yields the uniform strict bound 0.844086034386379 for this entire family, including arbitrary intersections of positive-turn hinge halfplanes. The two separate area ceilings need not be attained together. Run `python -B bound_fixed_hinge_family.py` to evaluate the bounds with outward endpoints.
