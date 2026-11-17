# r041：首次合并比例、原图根几何与 \(Q_n\) 的精确换测度

## 结论范围

本文件将 r037 的 top-two 比例律与原图根几何联合，但不声称已经证明任何一致可积性、
端点局部化或 \(Q_n\to0\)。严格得到的是：

1. 一个保留根位置的联合 \(B\)-law；
2. 从它到完整 \(K=2\) 条件律的精确非负核；
3. \(Q_n\to0\) 的精确必要充分几何判据；
4. 两类明确的、尚未证明的充分条件：常规 UI 路线和更弱的端点乘积几何路线。

## 1. 联合 top-two 律与根交换

令 \(G_n=P_n\square P_n\)，\(B(n)\) 为两个最高标号都是峰的双峰词数。在均匀
\(B\)-law 下，一个词的前两根为 \(r_1,r_2\)，其中 \(r_1\) 是全局高根。首次合并前，
它们所在两支的大小记为 \(a,b\ge1\)，并置

\[
q=a+b,\qquad \Theta=\frac aq,\qquad U=\frac1q,
\qquad D=\frac{\rho(r_1,r_2)}n,\qquad
E=\frac{\|r_1-r_2\|_2^2}{n^2}.
\tag{1}
\]

记这些量的律为 \(\lambda_n\)。它支撑于紧集
\[
[0,1]\times[0,2]\times[0,2]\times[0,1/2].
\]
更精确地，有限 \(n\) 支撑还满足 \(U\le\Theta\le1-U\)；以下涉及核的连续截断时，
均在这个闭子集上理解。
原图的几何关系是
\[
\frac{D^2}{2}\le E\le D^2.
\tag{2}
\]

交换词的前两个非邻根保持所有其后的已见顶点集合，故仍是一个 top-two 词；
它交换 \(a,b\)，但不改变 \(\rho,\|r_1-r_2\|_2\) 或 \(q\)。因此
\[
(\Theta,D,E,U)\overset{\rm law}=(1-\Theta,D,E,U).
\tag{3}
\]
这同时说明：固定两支内部词、首次合并点与合并后后缀后，所有余下交错都不改变两个根
的位置，因而也不改变 \(D,E\)。

r037 已严格证明 \(q\to\infty\) 依 \(B\)-law 概率收敛。因此每个
\(\lambda_n\) 的弱子序列极限 \(\lambda\) 满足 \(U=0\) 几乎处处，且仍有 (3) 的对称性。

## 2. 精确有限 \(n\) 核，包括 singleton 和 \(q=2\)

定义
\[
\kappa(\theta,u):=\frac{1-u}{2\theta(1-\theta)}.
\tag{4}
\]
在有限 \(n\) 的支撑上 \(0<\Theta<1\)，并且逐词严格有
\[
\kappa(\Theta,U)
=\frac{q(q-1)}{2ab}=:m(a,b).
\tag{5}
\]

这里不需要省略边界纤维：

- \(a=b=1,q=2\) 时，\(\Theta=1/2,U=1/2\)，(5) 给出 \(m=1\)，且 \(L=1\)；
- 一支 singleton、\(q>2\) 时，(5) 给出 \(m=q/2\)，而方向平均的 \(L\) 律在
  \(1,q-1\) 各有质量 \(1/2\)，故条件均值也是 \(q/2\)。

一般地，固定未记录根方向、两支内部词和后缀的纤维后，r022 的交错恒等式对任意
交换根不变的 \(H=H(r_1,r_2)\) 给出

\[
\mathbb E_B[LH\mid\mathcal D]=m(a,b)H.
\tag{6}
\]

特别地，令
\[
Z_n:=\int\kappa\,d\lambda_n=\mathbb E_{B,n}m=\frac{A_2(n)}{B(n)}
\tag{7}
\]
并定义非负 size-biased 概率测度
\[
d\widehat\lambda_n:=\frac{\kappa}{Z_n}\,d\lambda_n.
\tag{8}
\]
则其 \((D,E)\) 边际恰是原均匀标号条件律 \(\mathbb P_2\) 下的
\((\rho/n,\|r_1-r_2\|_2^2/n^2)\) 律。故以下全是精确有限 \(n\) 式：

\[
\boxed{\quad
Q_n=\mathbb E_2E
=\frac{\int E\kappa\,d\lambda_n}{\int\kappa\,d\lambda_n},
\qquad
\mathbb P_2(\rho>\varepsilon n)
=\frac{\int{\bf1}_{\{D>\varepsilon\}}\kappa\,d\lambda_n}{\int\kappa\,d\lambda_n}.
\quad}
\tag{9}
\]

这就是所需的首次合并加权空间几何；没有把 \(B\)-law 的几何误当作完整双峰律。

## 3. \(Q_n\) 的必要充分形式

由 (8)--(9) 和 \(0\le E\le2\)，以下等价：

\[
Q_n\longrightarrow0
\quad\Longleftrightarrow\quad
\widehat\lambda_n(E>\eta)\longrightarrow0
\quad\hbox{对每个 }\eta>0.
\tag{10}
\]

再由 (2)，这也等价于 \(\widehat\lambda_n(D>\varepsilon)\to0\) 对每个
\(\varepsilon>0\)，即原题的宏观 Manhattan 距离集中。量化地，
\[
\mathbb P_2(\rho>\varepsilon n)\le\frac{2Q_n}{\varepsilon^2}.
\tag{11}
\]

因此任何宣称 \(Q_n\to0\) 的路线必须证明的是 \(\kappa\)-size-biased 联合几何集中，
而不只是 \(\lambda_n\)（top-two）的 \(\Theta\) 或根几何弱收敛。

## 4. 弱收敛加 UI 时的充分且必要判据

取一个子列使 \(\lambda_n\Rightarrow\lambda\)。仅在可行闭集的两个端点
\((\theta,u)=(0,0)\) 或 \((1,0)\) 将 (4) 延拓为 \(+\infty\)；
对 \(0<\theta<1,u=0\)，仍取有限值 \([2\theta(1-\theta)]^{-1}\)。
对任意有限 K，截断核 \(\min(\kappa,K)\) 在该可行闭集上连续。考虑下列尚未证明的 UI 条件：

\[
\lim_{K\to\infty}\limsup_n
\int\kappa\,{\bf1}_{\{\kappa>K\}}\,d\lambda_n=0.
\tag{12}
\]

它强制 \(\lambda\) 不在 \(\theta=0,1\) 留正质量；在这个条件下，截断加弱收敛给
\[
\int\kappa\,d\lambda_n\to\int\kappa\,d\lambda,\qquad
\int E\kappa\,d\lambda_n\to\int E\kappa\,d\lambda.
\tag{13}
\]
所以
\[
Q_n\to
\frac{\int E[2\theta(1-\theta)]^{-1}\,d\lambda}
{\int[2\theta(1-\theta)]^{-1}\,d\lambda}.
\tag{14}
\]
在 (12) 下，这个子列的 \(Q_n\to0\) 当且仅当 \(E=0\) 在
\([2\theta(1-\theta)]^{-1}d\lambda\)-几乎处处成立。若每个弱子列都满足 (12)
且其极限具有这个性质，才可推出全序列 \(Q_n\to0\)。

式 (12) 是方便但强的充分工具；它不是本项目已证事实。

另一个已知形式、但同样未证的充分条件不要求 \(Z_n\) 收敛。令
\[
W_n:=\frac{\kappa}{Z_n},\qquad \int W_n\,d\lambda_n=1.
\]
若
\[
\lim_{K\to\infty}\limsup_n\int W_n{\bf1}_{\{W_n>K\}}\,d\lambda_n=0
\tag{15}
\]
且 \(\int E\,d\lambda_n\to0\)，则 \(Q_n\to0\)。事实上，将
\(\int EW_n\) 按 \(W_n\le K\) 与补集分开，前者至多
\(K\int E\,d\lambda_n\)，后者至多两倍 (15) 的尾项。这正是已有 r023 的
\(m/\mathbb E_Bm\) 一致可积路线以联合几何记号重写；它仍是一个待证输入。

## 5. 更弱的端点乘积几何充分条件

UI 并非必要，因为不平衡端点即使承载很大的 \(\kappa\)，也可能同时承载很小的根几何。
令
\[
I_\delta=[\delta,1-\delta],\qquad C_\delta=[0,1]\setminus I_\delta.
\]
以下两个条件足以推出 \(Q_n\to0\)：

\[
\text{对每个固定 }\delta>0,\qquad
\int_{\{\Theta\in I_\delta\}}E\,d\lambda_n\longrightarrow0,
\tag{16}
\]
\[
\lim_{\delta\downarrow0}\limsup_{n\to\infty}
\frac{\int_{\{\Theta\in C_\delta\}}E\kappa\,d\lambda_n}
{\int\kappa\,d\lambda_n}=0.
\tag{17}
\]

确实，在 \(I_\delta\) 上有 \(\kappa\le[2\delta(1-\delta)]^{-1}\)，而
\(Z_n=\int\kappa\,d\lambda_n\ge1\)。故 (16) 使内部部分的 (9) 分子除以 \(Z_n\)
趋零；(17) 消去端点部分。

相对于用 (12) 消去端点权重的路线，(17) 更弱：它允许端点的 \(\kappa\) 质量、
甚至 \(Z_n=A_2/B\)，发散，只要求带空间因子 \(E\) 的端点质量相对 \(Z_n\) 消失。近 \(\Theta=0\) 时
\[
\kappa E\sim\frac{E}{2\Theta}.
\tag{18}
\]
所以一种足够但尚未建立的局部化估计是端点区平均意义下 \(E=o(\Theta)\)，或更一般的
(17)。这正是可利用窄支本身使根距局部化的乘积几何版本。

现有资料没有给出 (17)：\(\rho\le q\) 是反向的一侧界，不能由 \(\Theta\) 控制
\(\rho/n\)；而平衡的顶行路径纤维已有 \(m<2\) 却有常数阶 \(E\) 的严格例子，
故单独控制不平衡或 UI 也不能替代内部几何条件 (16)。

## 6. 端点逃逸可造成非零极限的机制

弱极限 \(\lambda\) 单独不足以决定 (9)，因为它忘掉了 \(\Theta\to0,1\) 的速率以及
它和 \(E\) 的相关性。一个严格的诊断量是
\[
R_{\rm end}:=
\lim_{\delta\downarrow0}\liminf_{n\to\infty}
\frac{\int_{\{\Theta\in C_\delta\}}E\kappa\,d\lambda_n}{Z_n}.
\tag{19}
\]
若 \(R_{\rm end}>0\)，则 \(\liminf Q_n\ge R_{\rm end}>0\)，原题肯定结论失败。
这是条件性判据，不是已构造的方格反例。

例如，一个正 \(B\)-law 质量的 singleton 型端点序列有
\(\Theta=1/q\to0\)、\(\kappa=q/2\)。若其 \(E\) 不趋零且其 \(\kappa\)-质量在
\(Z_n\) 中不消失，它通过 (19) 保留正的 \(Q_n\) 贡献；反之 \(E\) 足够快地趋零时，
该同一端点 \(\kappa\) 发散仍可能不妨碍 \(Q_n\to0\)。所以端点质量本身既不能证明
也不能否证集中，必须保留乘积 \(E\kappa\)。

## 7. 首缺口

严格身份已将原题完全定位为 (9) 的 size-biased 联合几何。当前首个未证步骤是证明
(16)--(17)，或直接证明 (10) 的 \(\widehat\lambda_n\)-几何集中。没有已证的
首次合并大小、短间隙、粗边界或 \(B/A_1\) 估计给出这个乘积控制。

## 复核范围

解析输入为 r022 的首次合并 Rao--Blackwell 恒等式、r023 的 \(m\) 公式、r036 的
\(B/A_1\) 发散以及 r037 的 \(q\to\infty\)。audit_kernel.py 仅用 Fraction 核验
有限 \(a,b\) 纤维中 \(p_t\) 的归一化、\(\sum tp_t=m\) 与 (5)，包括 \(q=2\)；
其通过不承担 (9) 之后的任一渐近主张。
