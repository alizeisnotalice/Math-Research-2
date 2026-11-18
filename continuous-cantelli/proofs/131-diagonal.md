# Exact first negative Bernstein order in the three-layer model

For the positive non-graph source of round 5, the smallest total order \(j+k\) at which a nonnegative smooth compactly supported spatial test gives \(D_{j,k}(\psi)<0\) is **949**. The negative certificate for \((j,k)=(551,398)\) is in round 130. This round proves that every pair on \(j+k=948\) has positive spatial density on the entire real line. The identity
\[
 D_{j,k}=D_{j+1,k}+D_{j,k+1}
\]
then proves positivity of every lower order. The model has the exact Gaussian spatial marginal and complete exponential identities. It belongs to the non-graph relaxation.

## Real-space reduction

Use \(e=(3/5)e^{-16/3}\), and denote the centered Gaussian density of variance \(v\) by \(p_v\). The densities of \(Q_t\) between consecutive time layers are
\[
 p_t(y)\quad(0<t<1/2),\qquad
 p_t(y)-e p_{t-3/8}(y-2)\quad(1/2<t<1),\qquad
 e p_{t-7/8}(y-2)\quad(1<t<3/2).
\]
The endpoint values do not affect integration against \(W(t)=e^{-jt}(1-e^{-t})^k\). There is no remainder at times \(t\ge3/2\).

Only the second slab can be negative. The maximum in \(y\) of its source/free ratio is \((3/5)\sqrt{t/(t-3/8)}\), so negativity requires \(t<75/128\). The ratio at arbitrary \(y\) is that maximum times \(\exp[-K(t)(y-16t/3)^2]\), where \(K(t)=(3/8)/(2t(t-3/8))\ge1024/675\) on the potential negative range. The maximum is at most \(6/5\). Since \(\log(6/5)<1/5\) and \((1/5)/(1024/675)<(3/8)^2\), negativity is possible only when \(y\in[55/24,7/2]\). Outside that interval \(Q_t\) is nonnegative for all positive times and strictly positive during \((0,1/2)\).

For each pair, choose a cutoff \(c\le1/2\). Define positive and negative integrals
\[
 P_c(y)=\int_c^1 W(t)p_t(y)\,dt+e\int_1^{3/2}W(t)p_{t-7/8}(y-2)\,dt,\qquad
 N(y)=e\int_{1/2}^1W(t)p_{t-3/8}(y-2)\,dt.
\]
The full signed density is at least \(P_c(y)-N(y)\). All arithmetic for \(P_c,N\) drops the common positive factor \(1/\sqrt{2\pi}\) and scales by \(1/\sup W\).

## Two finite interval enclosures

**Rectangles.** With cutoff \(c=1/8\), split time at its layers and into dyadic cells. On a cell, \(W\) is unimodal and its global maximum is \((j/(j+k))^j(k/(j+k))^k\) for both indices positive, with maximum one for the single-zero index cases. The extrema of \(p_v(y-m)\) when \(y>m\) follow from its decrease in \(y\) and its unique variance maximum at \(v=(y-m)^2\). The outer bounds therefore give exact lower bounds on \(P_c\) and upper bounds on \(N\), uniformly on each spatial box. The rectangle pass certifies 6033 boxes, including the entire interval \([2.291,3.5]\) for \(j=0,\ldots,429\) and \(j=580,\ldots,948\).

**Certified spatial interpolation.** For the remaining boxes use cutoff \(c=3/8\) and \(F=\log(P_c/N)\). For a positive mixture of Gaussian terms, its logarithmic second derivative equals the weighted average of the log curvatures plus the variance of the log slopes. Here every positive-term slope in \(y\in[2.291,3.5]\) lies in \([-12,0]\). For the free term this follows from \(v\ge3/8, y\le3.5\), and for the late term from \(v\ge1/8, y-2\le1.5\). Every component log curvature is nonpositive, hence \((\log P_c)''\le36\) by the elementary variance bound on a bounded range. Those of \(N\) are at least \(-8\), so \((\log N)''\ge-8\) and \(F''\le44\). For spatial endpoints \(a,b\) with \(F\) lower bounds \(L_a,L_b\), this implies
\[
 F(y)\ge \min(L_a,L_b)-44(b-a)^2/8\qquad(a\le y\le b).
\]
The endpoint lower bounds use \(\log u\ge1-1/u\) and certified enclosures of the integrals. If the lower bound for \(P_c\) is \(P_l>0\) and the upper bound for \(N\) is \(N_u\), take \(L=1-N_u/P_l\).

The time integrals at each spatial endpoint are bounded separately with composite Simpson cells of width \(h=2^{-8}\) and error at most \(h^5\sup|\partial_t^4(Wg)|/2880\) in a single cell, where \(g=p_v(y-m)\) up to the common constant. If \(A,B,C,G\) are the first through fourth derivatives of \(\log(Wg)\), then the fourth derivative is \(Wg(A^4+6A^2B+3B^2+4AC+G)\). The derivatives of \(\log g\) on a term with variance \(v=t-b_0\) and squared distance \(d^2=(y-m)^2\) are
\[
 -\tfrac1{2v}+\tfrac{d^2}{2v^2},\quad \tfrac1{2v^2}-\tfrac{d^2}{v^3},\quad
 -\tfrac1{v^3}+\tfrac{3d^2}{v^4},\quad \tfrac3{v^4}-\tfrac{12d^2}{v^5}.
\]
With \(u=e^t\), the derivatives of \(\log W\) are
\[
 -j+\tfrac{k}{u-1},\quad -\tfrac{ku}{(u-1)^2},\quad
 \tfrac{ku(u+1)}{(u-1)^3},\quad -\tfrac{ku(u^2+4u+1)}{(u-1)^4}.
\]
Each derivative and each Gaussian is enclosed over its whole cell. Time-layer endpoints align with the Simpson subdivision. A separate ordinate mesh uses integer multiples of \(10^{-6}\), extending the input gap endpoints outward. The Simpson pass certifies 1052 spatial cells with smallest recorded log-ratio lower bound greater than \(3.28589\times10^{-6}\). All 1202 distinct spatial endpoints in this pass completed at depth 8.

An earlier independent calculation used midpoint error bounds and cutoff \(1/8\) in the 15 most dangerous neighborhoods. The same curvature argument gives \(F''\le204\) there since all positive logarithmic slopes are in \([-28,0]\). It certified 180 spatial cells, including \([2.926,2.938]\) for \((551,397)\), with lowest log-ratio lower bound greater than \(1.92230\times10^{-5}\). These boxes overlap the comprehensive proof.

All computations used 48-digit Decimal with outward rounding for basic arithmetic. Correctly rounded `exp` and `sqrt` were expanded by one adjacent representable value on both sides. Every spatial interpolation loss is a rational terminating decimal. `verify_coverage.py` merges the certified ranges for each integer \(j=0,\ldots,948\), checks every cell bound is strictly positive, and checks exact coverage of \([2.291,3.5]\). It found no gaps. Both the rectangle and the Simpson passes account for the entire proof.

## Compute and reproduction

The 15-worker MacBook CPU batches used Python 3.12.13 and standard-library Decimal. The initial 180-box pass had all-worker interval 07:20:00.018671–07:20:54.353383 UTC, totaling 678.522 CPU seconds. The rectangle pass had all-worker interval 07:23:05.372719–07:23:12.986958 UTC, totaling 113.680 CPU seconds. The final Simpson pass had all-worker interval 07:30:31.241797–07:30:40.462170 UTC, totaling 144.146 CPU seconds. All process exit receipts were zero, and separate process queries confirmed the PIDs absent. The 127 transferred file hashes were compared successfully. The 07:20 timed batch reported 0.0% idle in its three system snapshots. Its post-run clock bracket was +0.006586 to +0.115771 seconds for MacBook minus local.

In the remote working directory run `bash run.sh` for the initial neighborhood proof, `bash run_broad.sh` for rectangles, then `python make_remaining.py` and `bash run_simpson.sh` for the unresolved ranges. Fetch the output directories and run `python verify_coverage.py` to check each interval union. The archived MAC wrapper already contains its original UTC gate. For a run outside that gate, the workers start immediately when the scheduled UTC is in the past. The sources use no random seeds. The worker exit receipts and interval shards are archived with the computation.
