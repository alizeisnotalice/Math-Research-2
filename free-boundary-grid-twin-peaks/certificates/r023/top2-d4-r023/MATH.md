# Exact D4 inclusion–exclusion for tau=2 roots

For an independent set (I) of grid vertices, the recurrence computes the
integer (C_G(I)=N!\,\mathbb P(I\text{ are all peaks})), once for each D4
orbit representative.  This is the same recurrence used for the inherited
K=2 peak histogram.  Let

\[
B(I)=V\setminus N[I].
\]

For (d\ge2), let (T_d(I)) be the number of **ordered** pairs
((x,y)\in B(I)^2) with original-grid Manhattan distance (d).  Define
(T=\sum_dT_d); define (D,F,G) by summing respectively squared Euclidean
distance, doubled-centre norm sum, and doubled-centre dot product over the
same ordered pairs.  The exact numerator is accumulated one time per D4
representative using its orbit multiplier:

\[
 H_d=\frac1{(N)_2}\sum_{I\text{ independent}}(-1)^{|I|}C_G(I)T_d(I).
\]

The (I=\varnothing) term is included explicitly.  It is essential because
(B(\varnothing)=V).  The total (B_{\rm top2}=\sum_dH_d) counts the
complete tau=2 ordered-root words; both root directions are already in
(T_d), hence the division by ((N)_2=N(N-1)) is exactly once.

For each representative the code independently checks

\[
T=q(q-1)-2e_B,qquad
D=2\bigl(q\sum_{v\in B}|v|^2-|\sum_{v\in B}v|^2\bigr)-2e_B,qquad
F-2G=4D,
\]

where (q=|B|) and (e_B) is the number of unoriented grid edges inside
(B).  The final signed identity is checked again after all D4 layers.

The exact relation
\(\mathbb E_{\tau=2}L=A_2/B_{\rm top2}\) follows because completing the
second-root neighbourhood after (L) intervening ranks produces the
exactly-two-peak count (A_2).  This does **not** equate the tau=2 birth-rank
law to the law conditioned on exactly two peaks; `summarize_exact.py` reports
their rational (Q,V,C) quantities side by side.  Here
(Q=D/(n^2B)); the unnormalized root-distance moment (D/B) is separately named
`E_distance_squared`.  The checked identity is (Q=2(V-C)/n^2), and likewise
for K=2.
