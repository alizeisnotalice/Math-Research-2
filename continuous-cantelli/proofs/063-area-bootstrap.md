# Log-flat contact and the area bootstrap

Write $S(y)=\int g(u)K_{v(u)}(y-u)\,du$ with $K_v(z)=\int_0^v p_s(z)\,ds$, $g(0)>0$ and the exact Gaussian output density $P$. Then $S''=2(P-g)$ in distributions. Boundedness of $v$ and continuity of $K_v$ give a continuous $S$, so its analytic representative agrees with the pointwise value. Assume
$$
0\le v(u)\le C u^2/(\log(1/|u|))^k,\qquad k>0,
$$
on both sides of 0. Write $A(r)=\int_{-r}^r v(u)\,du$.

## Kernel and first boundary step

For $z\ne0$ and $v>0$,
$$
K_v(z)=2v p_v(z)-2|z|\Phi(-|z|/\sqrt v),\quad
K'_v(z)=-2\operatorname{sgn}(z)\Phi(-|z|/\sqrt v),\quad
K''_v(z)=2p_v(z),\quad K'''_v(z)=2p_v'(z).
$$
The derivative jump at zero is $-2$, giving the distributional Dirac term in $K''_v$. For $v=0$ use $K_0=0$. More generally all fixed real derivatives away from the mean obey
$$
|\partial_y^jK_{v(u)}(y-u)|
\le C_j v(u)|u|^{-j-1}e^{-c(\log(1/|u|))^k}\qquad (|y|\le|u|/4).
\tag{1}
$$
To see this, for $0<s\le v(u)$ put $a=|y-u|/\sqrt s$. The Gaussian derivative is bounded by $C_j |y-u|^{-j-1}$ times a polynomial in $a$ times $e^{-a^2/2}$. Extract $e^{-a^2/4}$ and use $a^2\gtrsim(\log(1/|u|))^k$ before integrating $s$. One may replace $v(u)$ by $C u^2$ in (1) for the first boundary estimate.

Let $a_j=\int g(u)\partial_y^jK_{v(u)}(-u)\,du/j!$ whenever this integral exists and set $b_j=[y^j]S-a_j$. By (1) the $a_j$ exist through degree 2. For $|y|\le r$, expand the integral over $|u|\ge8r$ to second order. The third derivative remainder from $8r\le |u|\le\sqrt r$ is $O(r^2 e^{-c'(\log(1/r))^k})$, whereas the region $|u|\ge\sqrt r$ gives $O(r^{5/2})$. Restoring the missing low coefficient integrals costs the same order up to adjusting $c'$. For $j=2$ the missing term is bounded by $Cr^2\int_{\log(1/r)-\log8}^\infty e^{-ct^k}\,dt=O(r^2e^{-c'(\log(1/r))^k})$. The cases $j=0,1$ also have an integrable power of $u$.

The inner potential is bounded by $\|g\|_\infty\sup_{|u|\le8r}v(u)=O(r^2(\log(1/r))^{-k})$ by monotonicity in $v$ and $\int K_v=v$. Comparing at $y=0,\pm r$ shows $b_0=b_1=b_2=0$. Integrating over $|y|\le2r$ captures a fixed fraction of the Gaussian mass for each $|u|\le r$ and time $s\le v(u)$. Thus the first boundary estimate gives $A(r)=o(r^3(\log(1/r))^{-N})$ for every fixed $N$.

## Dyadic improvement to a full power and iteration

Assume for some $m\ge2$ that $A(r)=O(r^{m+1})$ and $b_0,\ldots,b_m$ vanish. Let $r_j=2^{-j}$, $F_j=A(r_j)/r_j^{m+1}$, and $h_j=e^{-c j^k}$. On $E_j=\{r_{j+1}<|u|\le r_j\}$, (1) bounds the integrated order-$d$ derivative by
$$
C_dh_jr_j^{m-d}F_j.
$$
In particular coefficient integrals of order $d\le m$ are absolutely convergent because $\sum h_j<\infty$ and the sequence $F$ is bounded.

For $n\ge J+3$, split at $8r_n$ and consider $|y|\le2r_n$. In the exterior Taylor remainder of order $m+1$, the fixed region $|u|\ge r_J$ contributes $O_J(r_n^{m+1})$, and the other exterior annuli lie in $J\le j\le n-4$. Analyticity of $S$ bounds its own Taylor remainder by $O(r_n^{m+1})$. Restoring the omitted coefficients of degrees $0,\ldots,m$ uses $j\ge n-3$, and $r_j/r_n\le8$ in that range. Thus the inner potential divided by $r_n^m$ is bounded by
$$
C_{m,J}r_n
+C_m\sum_{j\ge n-3}h_jF_j
+C_m\sum_{J\le j\le n-4}2^{j-n}h_jF_j.
$$
The same bound holds for $F_n$ up to fixed constants: integrate the positive inner potential on $[-2r_n,2r_n]$ and use the Gaussian capture estimate from $|u|\le r_n$. For the three initial indices $n=J,J+1,J+2$, absorb $F_n$ into $C_{m,J}r_n$. This gives a positive vector inequality $F\le G+TF$, with $G_n=C_{m,J}r_n$ and the two sums defining $T$ on noninitial rows.

The operator norm is at most $C_m\sum_{j\ge J}h_j$ on ordinary $\ell^\infty$. It has the same bound up to a constant on the weighted space with norm $\sup_j |F_j|/r_j$. In the inward sum $r_j/r_n\le8$, and in the outward sum $2^{j-n}r_j/r_n=1$. Take $J$ large so both norms are below $1/2$. Iterate the vector inequality in ordinary $\ell^\infty$, where $F$ is initially bounded. The term $T^N F$ tends to zero, and the remaining Neumann series is bounded in the weighted norm. Consequently $F_n=O(r_n)$, so $A(r)=O(r^{m+2})$ for all small $r$ by monotonicity.

At order $m+1$ coefficient integrals exist by (1) and the improved area estimate. The same remainder bounds now show
$$
S_{\rm in}(r_n z)/r_n^{m+1}=b_{m+1}z^{m+1}+o(1)
$$
uniformly on bounded $z$ intervals, after splitting at any fixed larger multiple of $r_n$. The inward tail tends to zero by summability of $h_j$ and the outward sum tends to zero by bounded convergence. Positivity forces any odd $b_{m+1}$ to vanish, and gives $b_{m+1}\ge0$ when even. Starting at $m=2$ therefore gives $A(r)=O(r^5)$ for every log-flat contact under consideration.

If the first nonzero coefficient occurs at even degree $m\ge4$, the inner expansion gives the weak tangent $v(rz)/r^m\,dz\Rightarrow (b_m/g(0))z^m\,dz$ on bounded windows, for all $r\downarrow0$. Indeed test by a smooth compact function and move the Gaussian smoothing onto the test. Its added variance is at most $O((\log(1/r))^{-k})$ after rescaling, and the rescaled source areas are bounded. The dyadic expansion extends to arbitrary $r$ by comparable radii and the same fixed source window bounds. It follows that $A(r)\sim 2 b_m r^{m+1}/((m+1)g(0))$. If no nonzero coefficient occurs, the iteration proves $A(r)=O_N(r^N)$ for every $N$. This is an area result and does not require pointwise bounds on the rescaled variance.
