# Gibbs 2018 Definition-3 base polygon audit

Source reviewed: Philip Gibbs, *An Upper Bound for Lebesgue's Covering Problem*,
arXiv:1810.10089v1, pp. 4--5 (Definition 3, Figure 3, Proposition 2).
PDF SHA-256:
`df420d8754fd888a2adc8a1392d92f4d2f9ef4789adf263437c0afffad48712f`.

## Coordinate reconstruction

Use the regular hexagon of inradius `1/2` with centre zero and clockwise labels
starting at the top:

\[
\begin{aligned}
A&=(0,1/\sqrt3),&B&=(1/2,1/(2\sqrt3)),&C&=(1/2,-1/(2\sqrt3)),\\
D&=(0,-1/\sqrt3),&E&=(-1/2,-1/(2\sqrt3)),&F&=(-1/2,1/(2\sqrt3)).
\end{aligned}
\]

Figure 3 shows that the two rotated-hexagon supporting lines cutting the corners
labelled `C` and `E` have outward unit normals

\[
u_C=(\cos(\sigma-\pi/6),\sin(\sigma-\pi/6)),\qquad
u_E=(\cos(\sigma-5\pi/6),\sin(\sigma-5\pi/6)).
\]

Thus the only unambiguous base-polygon reconstruction supported jointly by the
text and Figure 3 is

\[
P(\sigma)=H\cap\{x:u_C\cdot x\le1/2\}\cap\{x:u_E\cdot x\le1/2\},
\]

where \(H=\bigcap_{j=0}^{5}\{x:(\cos(j\pi/3),\sin(j\pi/3))\cdot x\le1/2\}\).
The two normals belong to the hexagon rotated by \(\pi/6+\sigma\), exactly as
Definition 3 says.  They are related by the original hexagon's `-120` degree
symmetry, so the two removed triangles have equal area.

## Area formula for the base polygon

Put \(\beta=\sigma-\pi/6\),
\(c=C=(1/2,-1/(2\sqrt3))\), and
\(v=(-\sqrt3/2,-1/2)\), the unit direction from `C` to `D`.  The two side
offsets of the `C` cut are

\[
d(\sigma)=\frac{(1-\cos\beta)/2}{\sin\beta}+\frac1{2\sqrt3},
\qquad
\ell(\sigma)=\frac{1/2-u_C\cdot c}{u_C\cdot v}.
\]

The cut triangle has area \(\sqrt3d\ell/4\), giving

\[
|P(\sigma)|=\frac{\sqrt3}{2}-\frac{\sqrt3}{2}d(\sigma)\ell(\sigma).
\]

This is only the area before the later arc cuts `CS`, `ES`, `AS`, `AH`, and
`EH`; it is not Gibbs's final upper-bound area.  At \(\sigma=0\), it gives
\(2-2/\sqrt3=0.845299461620748\ldots\), the Pál base cover, as required.

## Actual local verification

`reconstruct_p_sigma.py` clips the polygon independently and compares its
shoelace area with the displayed formula.  Run:

```sh
python3 reconstruct_p_sigma.py
```

The checked angles `0`, `1`, and `1.5494` degrees all give eight vertices and
formula residual at most `2.22e-16`.  At `1.5494` degrees the base-polygon area
is `0.845351972027020`; this is larger than the final quoted value because the
five later regions have not been subtracted.

## Scope and remaining specification work

The reconstruction fixes the two Definition-3 half-planes only.  It does not
interpret later labelled intersections from Figures 4--15.  In particular,
Definition 8's shortest path in \(\mathcal R\) and Definition 9's curved
convex-hull path must be converted to branch-stable line/circle formulas before
the final area of \(\mathcal H(\sigma)\) can be reproduced.
