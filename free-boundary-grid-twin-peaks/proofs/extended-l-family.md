# Extended-L family: exact two-component reduction

For \(n\geq4\), let
\[
 S_n=\{(2,0),(0,1),(1,1),(2,1)\}
 \cup\{(x,0):3\leq x\leq n-1\}
 \cup\{(0,y):2\leq y\leq n-1\}.
\]
It is an induced Manhattan geodesic from \((n-1,0)\) to \((0,n-1)\):
the path has \(2n-1\) vertices and \(2n-2\) edges, exactly the endpoint
Manhattan distance.  Direct coordinate counting gives
\[
 |S_n|=2n-1,\quad m_n=|V\setminus S_n|=(n-1)^2,\quad
 |\partial S_n|=2n-2.
\]
The complement has exactly two components:
\[
 A=\{(0,0),(1,0)\},\qquad
 B_n=(V\setminus S_n)\setminus A,
\]
with \(|B_n|=n^2-2n-1\).  Both vertices of \(A\) are already frontier
vertices.  The boundary vertices in \(B_n\) number \(2n-4\); there are
also exactly \(2n-4\) ordered admissible new two-step pairs, all inside
\(B_n\).  The small program verify_extended_l_geometry.py checks these
coordinate identities for \(4\leq n\leq20\).

Write \(H_n=h(S_n)\), \(D_n=\sum_{x,y}h(S_nxy)\).  Let \(H_n^B\) count
legal relative completion orders of \(B_n\), and let \(D_n^B\) be the
sum over its admissible ordered two-step prefixes, with the same fixed
initial \(S_n\).  An interleaving of an arbitrary legal \(B_n\)-history
with the two freely insertable vertices of \(A\) is a bijection with a
legal full completion.  Therefore
\[
 H_n=\binom{m_n}{2}\,2!\,H_n^B=m_n(m_n-1)H_n^B. \tag{1}
\]
No valid new two-step pair begins in \(A\).  After a specified \(B_n\)
pair, the same interleaving bijection leaves \(m_n-2\) positions, hence
\[
 D_n=\binom{m_n-2}{2}\,2!\,D_n^B
 =(m_n-2)(m_n-3)D_n^B. \tag{2}
\]
Consequently the exact pointwise kernel has the reduction
\[
 |\partial S_n|p(S_n)=
 \lambda_n (2n-4)\frac{D_n^B}{H_n^B},\qquad
 \lambda_n=
 \frac{n-1}{n-2}
 \frac{(m_n-2)(m_n-3)}{m_n(m_n-1)}\longrightarrow1. \tag{3}
\]

This makes the genuine missing estimate explicit.  A proof that
\(\liminf_n D_n^B/H_n^B>0\) would give
\(p(S_n)\geq c>0\) and thus rule out **every** bound
\(p(S)\leq C/|\partial S|\) for fixed \(C\).  The interleaving identity
alone does not supply that lower bound; it only converts such a
large-component estimate to the original grid prefix.

## Exact finite instances

For the saved \(n=7\) truncation,
\[
 |\partial S_7|p(S_7)\approx3.4361735253823023,\quad
 (2n-4)D_7^B/H_7^B\approx3.215670411454026,\quad
 \lambda_7=187/175.
\]
For the certified \(n=8\) instance,
\[
 |\partial S_8|p(S_8)\approx4.0552255354842766,\quad
 (2n-4)D_8^B/H_8^B\approx3.781375892477475,\quad
 \lambda_8=1081/1008.
\]
These are exact-rational evaluations from saved integers, but two values
are not an asymptotic lower-bound proof.  In particular, they do not
justify a claim that this family refutes all \(C/|\partial S|\) bounds.
