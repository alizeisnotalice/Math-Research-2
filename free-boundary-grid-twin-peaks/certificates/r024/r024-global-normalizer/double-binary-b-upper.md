# A uniform one-peak normalization bound for the double-binary bridge

Let \(D_h\) be two height-\(h\) complete rooted binary trees joined at
their roots.  Set \(S=2^{h+1}-1\) and \(N=2S\).  For a tree \(T\), let
\(a_r(T)\) be its one-peak count with prescribed root \(r\), so
\(A_1(T)=\sum_r a_r(T)\).

We use the exact tree cut-edge formula established separately:

\[
B(T)=\sum_{e=uv,\ T-e=A\sqcup B}
 {N-2\choose |A|-1}\bigl(A_1(A)A_1(B)-a_u(A)a_v(B)\bigr). \tag{E}
\]

For every non-bridge edge of \(D_h\), call the smaller, complete descendant
binary component \(B\), with \(s=|B|\le(S-1)/2\); for the bridge, choose
either full block as \(B\).  The following proves the useful analytic
upper bound

\[
 \boxed{B(D_h)\le720(h+1)A_1(D_h).} \tag{U}
\]

## Root ratios

The hook formula has the adjacent-root consequence

\[
 {a_{\rm child}(T)\over a_{\rm parent}(T)}={q\over M-q}, \tag{R}
\]

where the child-side component has size \(q\) in an \(M\)-vertex tree.
For a complete binary tree rooted at its top, a depth-\(d\) vertex has
root ratio at most
\[
 \prod_{j=1}^d{1\over2^j-1}.
\]
After multiplication by the \(2^d\) vertices at that depth, the first two
layers contribute \(2\) and \(4/3\), and the later layer contributions
have successive ratio at most \(2/7\).  Hence

\[
 A_1(B)/a_v(B)<5. \tag{B1}
\]

Now delete a descendant binary subtree of size \(s\) from one block and
write \(A=D_h\setminus B\), of size \(M=N-s>3N/4\).  Choose as \(r\) the
root of the other, intact block.  At graph distance one from \(r\), the
outward component has at most \(M/2\) vertices.  At level \(\ell\ge2\),
it has at most \(N/2^\ell\) vertices, so (R) gives a factor at most
\[
 {4\over3\,2^\ell-4};
\]
this is at most \(1/2\) for \(\ell=2\) and at most \(1/5\) thereafter.
There are at most \(3\,2^{\ell-1}\) vertices at level \(\ell\).  Summing
the resulting geometric majorant gives

\[
 A_1(A)/a_r(A)<1+3+3+{3(2/5)\over1-2/5}=9. \tag{A1}
\]

For a bridge cut both components are complete binary trees, so (B1) also
supplies the corresponding bound for its large component.

## The ancestral hook correction

One must not omit the effect of ancestors of the cut root.  If the tree is
rooted at \(r\in A\), hook cancellation yields

\[
a_r(D_h)={N-1\choose s}a_r(A)a_v(B)G,
\quad
G=\prod_{w\text{ nonroot ancestor of }v}{t_w-s\over t_w}. \tag{G}
\]

Here \(t_w\) is the original rooted subtree size at \(w\).  At the
\(j\)-th such ancestor, \(t_w\ge2^js+2^j-1\); consequently

\[
G^{-1}\le\prod_{j\ge1}{2^j\over2^j-1}<8. \tag{G1}
\]

The last numerical bound follows, for example, from
\(\log(1/(1-2^{-j}))\le2^{1-j}\) and summation.  For the bridge there are
no such ancestors.

Divide the positive part of an edge summand in (E) by \(a_r(D_h)\), use
(G), discard the negative term in (E), and cancel factorials.  Equations
(B1), (A1), and (G1) give

\[
 {\text{edge summand}\over a_r(D_h)}
 \le360\,{s\over N-1}. \tag{H}
\]

Finally, summing the selected smaller descendant-subtree sizes over all
edges counts every vertex at most once per binary depth, plus the bridge:

\[
 \sum_e s_e\le N(h+1).
\]

Since \(a_r(D_h)\le A_1(D_h)\), summing (H) proves (U).

## The general promotion lower bound closes the comparison

The already proved weighted-promotion theorem applies to every connected
maximum-degree-\(\Delta\) graph.  Its \(r=0\) specialization is

\[
 {A_2(T)\over A_1(T)}\ge
 \sum_{q=1}^{N-2}{(N-1-\Delta q)_+\over N-q-1}. \tag{P}
\]

For \(1\le q\le\lfloor(N-1)/(2\Delta)\rfloor\), each summand in (P)
is at least \(1/4\).  Thus, for \(N\ge4\Delta+2\),

\[
 {A_2(T)\over A_1(T)}\ge {N\over16\Delta}. \tag{L}
\]

Apply (L) with \(\Delta=3\) to \(D_h\), whose order is
\(N_h=2^{h+2}-2\).  Combining it with (U) gives, for every \(h\ge2\),

\[
 \boxed{\ {A_2(D_h)\over B(D_h)}
 \ge {N_h\over34560(h+1)}\ \longrightarrow\infty.\ } \tag{D}
\]

Consequently there is no constant depending only on maximum degree for
which \(A_2(T)\le C(\Delta)B(T)\) holds for every finite connected graph:
the degree-three trees \(D_h\) already contradict it.  This is an
auxiliary graph-class counterexample, not a counterexample to the free
square-grid twin-peaks question.

The source for (P) is the independently archived proof
`free-boundary-grid-twin-peaks/proofs/weighted-promotion.md`, SHA-256
`26c13189b3d3a9957bf176cccfff755c8c226badbb765cf0f6a9ed06e4b1bc83`.
