# Exact degree-two calibration

Use ranks \(1,\ldots,n\), where rank 1 is the largest label.  Thus a
peak is a local rank minimum.  Let \(a_{n,k}\) count rank permutations of
the path \(P_n\) having exactly \(k\) peaks, and let \(b_n\) count those
with exactly two peaks for which ranks 1 and 2 are both peaks.

## Path theorem

For \(n\geq3\),

\[
 a_{n,1}=2^{n-1},\qquad
 a_{n,2}=2^{2n-3}-n2^{n-2},\qquad
 b_n=2^{2n-4}-2^{n-2}.
\]

Consequently

\[
 \frac{a_{n,2}}{b_n}=2-\frac{n-2}{2^{n-2}-1}<2. \tag{P}
\]

### Proof

Delete the unique rank-1 vertex and join its former path neighbours.  The
inverse inserts a new rank-1 vertex in one of the \(n\) spatial slots of a
path of length \(n-1\): its two exterior slots and its \(n-2\) edges.
For a word with \(k\) peaks, precisely \(2k\) slots preserve the number
of peaks: the slots incident to a peak.  An endpoint peak has its exterior
slot and its unique incident edge, while an internal peak has its two
incident edges.  Insertion in one of those slots replaces that peak by the
new rank-1 peak.  All remaining slots make the new vertex a new peak
without destroying an old one.  (The incident-slot sets of distinct local
rank minima are disjoint because peaks are non-adjacent.)

Thus, with the zero conventions outside the possible range,

\[
 a_{n,k}=2k\,a_{n-1,k}+(n-2k+2)a_{n-1,k-1}. \tag{1}
\]

The initial values are \(a_{2,1}=2\) and \(a_{2,k}=0\) for \(k\ne1\).
Putting \(k=1,2\) in (1) and solving the two first-order recurrences gives
the displayed formulas for \(a_{n,1}\) and \(a_{n,2}\).

For \(b_n\), again delete rank 1.  An old word with one peak has exactly
\(n-2\) valid slots: every slot except the two incident to its sole peak;
the insertion then produces the two top-ranked peaks.  An old word with
two peaks has exactly two valid slots, namely those incident to the peak
other than its rank-1 peak.  Words with three or more peaks cannot become
two-peak words, since an inserted rank-1 vertex destroys at most one old
peak.  Therefore

\[
 b_n=(n-2)a_{n-1,1}+2a_{n-1,2}.
\]

Substitution of the two formulas above gives the asserted expression for
\(b_n\), hence (P).  This also verifies the equivalent recurrence
\(b_n=4b_{n-1}+a_{n-1,1}\), with \(b_3=2\).

## Cycle theorem

For a cycle, divide labelled counts by \(n\), equivalently mark and then
forget a cyclic origin.  Write \(c_{n,k}=A_k(C_n)/n\) and
\(d_n=B(C_n)/n\).  Deleting rank 1 leaves \(n-1\) cyclic insertion slots.
The same incident-edge argument gives

\[
 c_{n,k}=2k c_{n-1,k}+(n+1-2k)c_{n-1,k-1}.
\]

Starting with \(c_{3,1}=2,c_{3,2}=0\), one obtains for \(n\ge4\)

\[
 c_{n,1}=2^{n-2},\qquad
 c_{n,2}=2^{2n-5}-(n-1)2^{n-3}.
\]

The rank-1 deletion classification is identical except that a one-peak
cycle has \(n-3\) valid slots.  Hence

\[
 d_n=(n-3)c_{n-1,1}+2c_{n-1,2}
     =2^{2n-6}-2^{n-3},
\]

and therefore

\[
 \frac{A_2(C_n)}{B(C_n)}
 =2-\frac{n-3}{2^{n-3}-1}<2 \quad(n\ge4).
\]

The cases before two non-adjacent vertices exist are excluded.  These
degree-two identities are a calibration theorem only; they do not establish
the desired bound for all bounded-degree graphs or for square grids.
