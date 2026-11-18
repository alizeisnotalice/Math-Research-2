# The Gaussian-capped moment problem on a fixed interval

Let \(I=[l,r]\), \(l<r\), let \(g=p_t\) be the centered Gaussian density of variance \(t>0\), and impose
\[
 0\le q\le g,\qquad \operatorname{supp}q\subset I,\qquad
 \int q=m,\qquad \int xq(x)\,dx=\mu .
\]
The conclusions for threshold existence and convex costs use only that \(g\) is continuous and strictly positive on \(I\).

Let \(M=\int_I g\). The cases \(m=0\) and \(m=M\) determine \(q\) uniquely. For \(0<m<M\), the feasible first-moment interval is precisely the one obtained by filling mass \(m\) entirely from the left and entirely from the right, respectively. The lower bound follows by subtracting a linear cost with its zero at the one-sided mass threshold. The upper bound follows by reflection. To see that all intermediate moments are attained, use the following two-tail construction.

For \(p\in[0,m]\), let \(a(p)\) and \(b(p)\) be the unique thresholds with
\[
 \int_l^{a(p)}g=p,\qquad \int_{b(p)}^r g=m-p.
\]
Since \(m<M\), \(a(p)<b(p)\) for all \(p\). The two-tail first moment is a continuous function \(\mu(p)\). Its derivative for \(0<p<m\) is \(a(p)-b(p)<0\), since \(a'(p)=1/g(a(p))\) and \(b'(p)=1/g(b(p))\). Its extreme values are the right-only and left-only first moments. Thus there is a unique allocation for every feasible \(\mu\). Write \(a,b\) for its thresholds.

The unique maximizing density for the second moment is
\[
 q_*(x)=g(x)1_{[l,a]\cup[b,r]}(x).
\]
Indeed the sign of \(P(x)=(x-a)(x-b)\) agrees with the cap-filled region, and
\[
 \int x^2(q_*-q)\,dx=\int P(x)(q_*-q)\,dx
 =\int |P(x)|\,|q_*-q|\,dx\ \ge0.
\]
Equality forces \(q=q_*\) almost everywhere. Equivalently the attained dual is
\[
 V(m,\mu)=-ab\,m+(a+b)\mu+\int_l^r[(x-a)(x-b)]_+g(x)\,dx.
\]
For the Gaussian one can instead evaluate
\[
 \mu=t[g(l)-g(a)+g(b)-g(r)],\qquad
 V=t[m+lg(l)-ag(a)+bg(b)-rg(r)].
\]

The same \(q_*\) simultaneously maximizes all convex costs. For a convex function replace \(P\) by its difference from the affine secant through \(a,b\), which is nonnegative outside \([a,b]\) and nonpositive inside. By complementing the cap, a single contiguous filled interval simultaneously minimizes all convex costs. These conclusions assume fixed mass, fixed first moment, and the density lower bound \(q\ge0\).

For an interior moment pair the map of thresholds is nonsingular, because
\[
 dm=g(a)\,da-g(b)\,db,\qquad
 d\mu=ag(a)\,da-bg(b)\,db.
\]
Consequently \(V\) is differentiable there with \(V_m=-ab,\ V_\mu=a+b\), and its Hessian quadratic form is
\[
 dm\,dV_m+d\mu\,dV_\mu
 =(a-b)\big[g(a)(da)^2+g(b)(db)^2\big]<0
\]
on nonzero perturbations. In particular \(V_{\mu\mu}=-(g(a)^{-1}+g(b)^{-1})/(b-a)\).

There is also a stability estimate. Put \(G=\|g\|_{\infty,I}\), \(D=V-\int x^2q\), and \(N=\|q-q_*\|_1\). If \(d(x)\) is the distance to the nearer threshold, then \(|P(x)|\ge(b-a)d(x)/2\). The mass of \(|q-q_*|\) in \(d(x)\le s\) is at most \(4Gs\). Layer-cake integration yields \(\int d(x)|q-q_*|(x)\,dx\ge N^2/(8G)\). Therefore
\[
 N\le4\sqrt{GD/(b-a)}.
\]

## A signed measure capped from above

If there is no nonnegativity assumption on \(q\), the two-tail conclusion need not hold. Write \(E=g1_I-q\), a nonnegative defect. Its mass is \(d=M-m\) and its first moment is \(G_1-\mu\), where \(G_j=\int_I x^j g(x)\,dx\). When \(d>0\) and \(x_c=(G_1-\mu)/d\) is feasible in \(I\), Jensen gives
\[
 \int x^2 q\le G_2-\frac{(G_1-\mu)^2}{M-m}.
\]
Equality as signed measures occurs at \(q=g1_I-d\delta_{x_c}\). For absolutely continuous signed \(q\), the bound is a supremum when \(x_c\) is in the interior, approached by defects of mass \(d\) supported on shrinking symmetric intervals about \(x_c\). It is strict for every such density when \(d>0\). A mean exactly at an endpoint cannot be realized by a nonzero absolutely continuous defect. If \(d=0\), the defect is zero.

## Exact numerical certificate

For \(I=[-2,3]\), \(t=1\), and unnormalized mass \(1/2\) and first moment \(-1/8\) relative to \(e^{-x^2/2}\), the exact rational certificate gives
\[
 a\simeq-1.02656167158426,\qquad
 b\simeq 1.46648784378548,\qquad
 V\simeq0.51961086193385
\]
after Gaussian normalization for \(V\). The certificate uses rational alternating bounds on the Gaussian density and its antiderivative, rational enclosures for both threshold quantiles at two mass splits, a strict bracket on the target mean, and Machin's arctangent formula to enclose \(\sqrt{2\pi}\). The final interval width is less than \(6.23\cdot10^{-14}\).
