# r022: projection skeletons behind the binomial-square root comparison

Let \(\mathcal W_{i,j}\) be the connected orderings of
\(P_n\square P_n\) whose first vertex is \((i,j)\), and write
\(a_{i,j}=|\mathcal W_{i,j}|\).  This note gives a rigorous decomposition
which explains the candidate

\[
 {a_{i+1,j}\over a_{i,j}}
 \leq\left({n-1-i\over i+1}\right)^2, \tag{1}
\]

but does not prove it.  It identifies a specific missing compression
inequality for *full squares*, while recording why the familiar local and
deletion inductions do not supply it.

## Two one-dimensional skeletons of every connected ordering

For \(w=(v_1,\ldots,v_N)\in\mathcal W_{i,j}\), let \(S_t\) be its first
\(t\) vertices and let

\[
 X_t=\{x:\ (x,y)\in S_t\text{ for some }y\}.
\]

Because \(S_t\) is connected, \(X_t\) is an interval of integers.  Every
time a new column first appears it is therefore the immediate left or
right neighbour of the previous interval.  The resulting first-arrival
word \(\alpha_x(w)\) has length \(n-1\), with exactly \(i\) letters L and
\(n-1-i\) letters R.

The same prefix fact gives a second, independent-in-definition skeleton.
Read the ordering backwards as a deletion order: after removing
\(v_N,\ldots,v_{t+1}\), the remaining set is precisely \(S_t\), hence is
connected.  Record L or R whenever an extreme column becomes empty.  This
last-departure word \(\beta_x(w)\) also has \(i\) L letters and
\(n-1-i\) R letters.  This is valid even though arbitrary suffixes of the
original ordering need not be connected.

The horizontal projection hence partitions

\[
 \mathcal W_{i,j}=
 \bigsqcup_{\alpha,\beta\in\{L,R\}^{n-1}\atop
             |\alpha|_L=|\beta|_L=i}
 \mathcal W_{i,j}^{\alpha,\beta}. \tag{2}
\]

Applying the identical construction to the vertical projection yields two
more skeletons, each with \(j\) D letters and \(n-1-j\) U letters.
Consequently there are at most

\[
 {n-1\choose i}^{\!2}{n-1\choose j}^{\!2} \tag{3}
\]

possible quadruples of projection skeletons for a root \((i,j)\).  Some
formally possible quadruples can be empty: no assertion of simultaneous
realizability is being made.

All assertions through (3) are strict combinatorial facts, requiring no
probability approximation or asymptotic passage.

## A sufficient compression statement, and its sharp failure

If for every admissible four-skeleton \(\sigma\) there were an injection

\[
 \mathcal W_{i,j}^{\sigma}\hookrightarrow\mathcal W_{0,0}, \tag{4}
\]

then (2)--(3) would give

\[
 a_{i,j}\leq {n-1\choose i}^{\!2}{n-1\choose j}^{\!2}a_{0,0}. \tag{5}
\]

Vandermonde's identity would then give

\[
 {A_1\over a_{0,0}}
 \leq\left(\sum_{r=0}^{n-1}{n-1\choose r}^2\right)^2
 ={2n-2\choose n-1}^2\leq16^{n-1}. \tag{6}
\]

Thus (4), or even the weaker uniform bound
\(|\mathcal W_{i,j}^{\sigma}|\leq C^n a_{0,0}\), is sufficient for the
desired \(\exp[-O(n)]\) corner-root lower bound.  The latter weaker
statement already tolerates exponential interface information; it does not
need a bounded fibre for every local root move.

The pointwise injection (4) is false already on the full 3-by-3 square.
An exhaustive connected-ordering enumeration gives

\[
\begin{array}{c|c|c|c}
\text{root}&\text{four-skeleton}&
|\mathcal W_{i,j}^{\sigma}|&a_{0,0}\\ \hline
(1,1)&(\mathrm{LR},\mathrm{RL},\mathrm{LR},\mathrm{RL})&
2352&1344 .
\end{array}
\]

Thus this one fixed class is \(7/4\) times the entire corner-root family,
so no injection in (4) can exist.  The calculation also finds that the
centre has four nonempty four-skeleton classes, each of size \(2352\);
there were only *at most* sixteen formally possible classes.  This is a
finite exact counterexample to the stronger classwise compression, not to
(1) or to the original twin-peaks problem.

The remaining potentially useful statement is the genuinely weaker
\(C^n\)-weighted compression, or an inequality that averages class sizes
before comparing roots.  It is not proved here.  The exhaustive n=3 audit
is in audit_n3_four_skeletons.py and n3-four-skeleton-audit.json.  It was
actually run with Python 3.15.0 using

    python3 -B audit_n3_four_skeletons.py

and enumerates exactly the 9! label orders through the connected-prefix
recursion; it does not use the independent-set dynamic program.

## Why the most direct compression maps fail

The tempting relabelling which maps a column to its rank of first
appearance is not a graph homomorphism.  For example on the 3-by-3 square,
start at \((1,0)\) and begin

\[
(1,0),(2,0),(1,1),(0,1),\ldots.
\]

The first-column order is \(1,2,0\), so rank relabelling sends these first
four vertices to

\[
(0,0),(1,0),(0,1),(2,1).
\]

Its last vertex has no earlier grid neighbour.  Thus this natural attempt
does not map connected words to connected words, even before injectivity is
considered.  A global reflection only sends the root \(i+1\) to
\(n-2-i\), not to \(i\), and supplies symmetry rather than (1).

Nor can one prove (4) by applying a deletion induction to arbitrary
connected induced subgraphs.  The already verified induced path
\[
(1,0)-(2,0)-(2,1)-(3,1)\subset P_5\square P_5
\]
has root counts \(1,3,3,1\).  On its first edge the ratio is \(3\), larger
than the full-square right side \(9/4\) corresponding to \(n=5,i=1\).
This does not refute (1), but it excludes the hereditary induction needed
by that route.

## Current status

The strongest proved content of this note is the four-skeleton
decomposition (2)--(3), the conditional implication (4)--(6), and the
explicit n=3 refutation of (4).  The first unproved viable step is a
comparable averaged inequality with \(C^n\) loss.  Finite exact root-vector
evidence through \(n=7\) supports (1), but is not used as proof here.
