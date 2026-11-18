# Supplement: every Sprague fan is contained in `P`

This is a post-freeze supplement to `certify_sprague_fixed.py`; it does not
alter that certificate or its hash.

At \(\sigma=1.5494^\circ\), write
\[
P=\bigcap_{j=0}^5\{n_j\cdot x\le1/2\}\cap
  \{m_3\cdot x\le1/2\}\cap\{m_5\cdot x\le1/2\}.
\]
The supplementary certificate checks every one of these eight halfplanes for
each of the five fans used in the main proof:
\[
\operatorname{Fan}(C_2;\widehat{KG}_{F_3}),\quad
\operatorname{Fan}(E_2;\widehat{IJ}_{B_3}),\quad
\operatorname{Fan}(E_2;\widehat{JH}_{C_3}),\quad
\operatorname{Fan}(A_1;\widehat{QX}_{E_3}),\quad
\operatorname{Fan}(A_1;\widehat{XP}_{G}).
\]

For an ordinary strict support inequality on a minor circular arc, the script
uses 256 rational radial-chord cells.  For a boundary equality it uses the
analytic tangent identity, rather than accepting a rounded zero:
\[
K-F_3=m_5,\quad I-B_3=m_3,\quad H-C_3=n_3,\quad
Q-E_3=n_1,\quad P-G=n_2.
\]
Thus the appropriate support functional is bounded on the *entire* unit
circle by its tangent value.  The one non-tangent boundary incidence,
\(G\in\{n_5\cdot x=1/2\}\), is handled by a strictly signed derivative of
its slack along the complete selected arc; the slack decreases to zero at
\(G\).  Apex coordinates are checked against all eight halfplanes (strictly
except for their stated incident edges).  Since a halfplane is convex, an arc
in `P` and its apex in `P` imply the full fan is in `P`.
