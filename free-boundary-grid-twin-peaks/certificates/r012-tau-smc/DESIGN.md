# 原条件律的第二出生秩 τ 分层 SMC

对象是 (G_n=P_n\square P_n)，(N=n^2)。令

\[
X=(x_1,\ldots,x_N)
\]

是按标签从高到低的均匀顶点排列，(S_r=\{x_1,\ldots,x_r\})。在第
(r) 步出生指示量是

\[
b_r=\mathbf1\{x_r\notin N(S_{r-1})\}.
\]

一次出生恰是一座峰；令

\[
T_2=\min\{r:\sum_{j\le r}b_j=2\}.
\]

故 (K_n=2) 可分割为互斥事件

\[
E_\tau=\{T_2=\tau,K_n=2\},\qquad 2\le\tau\le N.
\]

这里的距离一直是原方格距离；但空间矩使用平方欧氏距离
\(d_E^2/n^2\)，因为 \(\rho^2\le2d_E^2\)。

## 固定 τ 的受限提议与无偏恒等式

对未访问集合中的顶点，写

\[
B(S)=N(S)\setminus S,\qquad C(S)=V\setminus N[S].
\]

固定 τ 时，第 (r) 步允许集是

\[
A_r^{(\tau)}(S)=
\begin{cases}
V,&r=1,\\
C(S),&r=\tau,\\
B(S),&r\ne1,\tau.
\end{cases}
\]

从 (A_r^{(\tau)}) 均匀取下一个顶点；若允许集空则该路径权重为零。
这恰好生成 (E_\tau) 中的排列。相对于均匀排列参考律，其第 (r)
步 Feynman--Kac 势为

\[
G_r=\frac{|A_r^{(\tau)}(S_{r-1})|}{N-r+1}\in[0,1].
\]

于是对任何末端函数 (f)，其中两根是 (x_1,x_\tau)，有严格恒等式

\[
\mathbb E_{q_\tau}\!\left[\prod_{r=1}^N G_r\,f\right]
=\mathbb E_{\rm unif}[\mathbf1_{E_\tau}f]. \tag{1}
\]

理由是每条允许路径上，提议密度是
\(\prod_r|A_r|^{-1}\)，势的乘积是
\(\prod_r |A_r|/(N-r+1)\)，两者相乘为 (1/N!\)。这不是固定前缀
completion kernel，而是原均匀全排列律的分层重要性恒等式。

特别令

\[
Z_\tau=\mathbb P(E_\tau),\qquad
U_\tau=\mathbb E_{\rm unif}\!\left[
\mathbf1_{E_\tau}\frac{d_E(x_1,x_\tau)^2}{n^2}\right].
\]

则

\[
\mathbb P(K_n=2)=\sum_\tau Z_\tau=\frac{A_2(n)}{N!},
\qquad
Q_n=\frac{\sum_\tau U_\tau}{\sum_\tau Z_\tau}. \tag{2}
\]

因而 (2) 的 (Q_n) 正是原问题所需的条件空间矩；再由
\(\rho^2\le2d_E^2\) 得到既有的尾界。

## Bootstrap SMC 及末步不重采样

固定 τ，粒子在第 (r) 个 mutation 前持有 (S_{r-1})。先计算全部
(G_r^i)，令 \(\bar G_r=P^{-1}\sum_iG_r^i\)。若 \(r<N\)，按

\[
\mathbb E[N_i\mid \mathcal F_r]
=P\frac{G_r^i}{\sum_jG_r^j}
\]

选择父粒子，再由相应 (A_r) 均匀 mutation。多项式重采样满足这个条件；
正确实现的系统重采样也有相同的边缘 offspring 期望，但子粒子不独立。

最后一步不再重采样。令最终 mutation 后的函数为
\(f=d_E^2/n^2\)，则使用

\[
\widehat Z_\tau=
\left(\prod_{r=1}^{N-1}\bar G_r\right)\bar G_N,
\]

\[
\widehat U_\tau=
\left(\prod_{r=1}^{N-1}\bar G_r\right)
\frac1P\sum_{i=1}^P G_N^i f_i. \tag{3}
\]

在理想随机数和上述 offspring 条件下，(3) 分别是 (Z_\tau,U_\tau)
的无偏估计。实际有限 PRNG 只提供数值实现，不能据此声称随机性证明。
本实现给每个 child 用 `(replicate, rank, child-index)` 的新键生成未来
RNG 状态，绝不复制父粒子的未来 RNG 状态。

无偏性的归纳不需要把 particle 当作独立，也不要求每个 child 有相同的
父粒子分布。设
\(\widehat\gamma_{r-1}\) 是上一步的未归一化粒子测度，且在所有测试函数
\(h\) 上 \(\mathbb E\widehat\gamma_{r-1}(h)=\gamma_{r-1}(h)\)。给定
当前粒子，若第 (i) 个势为 (G_i)，mutation 核为 (Q_i)，令 (N_i) 为该父
粒子的 offspring 数。所需的重采样条件是

\[
\mathbb E[N_i\mid\mathcal F_r]
=P\frac{G_i}{\sum_jG_j}.
\]

若给定 parent 后每个 mutation 的条件均值为 (Q_i h)，则取全体 child
平均有

\[
\mathbb E\left[\frac1P\sum_{\text{children }c}h(X_c)\middle|\mathcal F_r\right]
=\frac1P\sum_i\mathbb E[N_i\mid\mathcal F_r]Q_i h
=\frac{\sum_iG_iQ_i h}{\sum_jG_j}.
\]

乘回 \(\bar G=P^{-1}\sum_iG_i\) 后，条件均值正好是
\(P^{-1}\sum_iG_iQ_i h\)。再对前一层取期望，归纳给出
\(\mathbb E\widehat\gamma_r(h)=\gamma_r(h)\)。第 (r=N) 步的加权均值
有同一条件均值，故不重采样也保持该结论。这证明 (3) 的一阶无偏性；它
不声称 ratio 或有限 PRNG 的无偏性。

整数实现中令 (a_i=|A_r(X_i)|)、\(T=\sum_i a_i\)。每一个 child 从
整数区间 \(\{0,\ldots,T-1\}\) 均匀取目标；属于第 (i) 个连续区间的
目标恰有 (a_i) 个，故其父粒子的边缘概率严格为 (a_i/T)。GPU 先求每
256 粒子的整数和，CPU 只读取这些组和定位组，GPU 在组选中的区间扫描。
因此没有浮点 CDF 舍入。64 位拒绝采样给 child 的顶点选择避免模偏差。

若某一阶段所有势均为零，整次重复规定
\(\widehat Z_\tau=\widehat U_\tau=0\)，并保留在重复平均中；不能删去
这些死端重复后再归一化。

重复 (R) 次后，\(R^{-1}\sum\widehat Z\) 与
\(R^{-1}\sum\widehat U\) 分别是无偏的未归一化估计；两者之比只是数值
条件矩诊断，并非无偏估计。对单个 τ 有

\[
0\le\widehat Z_\tau\le1,\qquad0\le\widehat U_\tau\le2,
\]

所以在理想独立重复下，Hoeffding 给出
\(\Pr(|\bar Z-Z|\ge a)\le2e^{-2Ra^2}\) 与
\(\Pr(|\bar U-U|\ge a)\le2e^{-Ra^2/2}\)。这是粗的重复方差界；
对多 τ 相加需保留分层协方差或采用带已知 τ 概率的混合设计，不能把
各分层比值直接平均。

建议每次报告：每个 τ 的死端率、最小单步不同直接父索引数（不是跨代谱系祖先数）、跨重复的
τ-normalizer ESS、\(\widehat Z^2\) 与 \(\widehat U^2\) 样本矩、以及
不同 P/R 下的稳定性。它们是退化诊断，不是严格误差证书。

## 本地 CPU 锚点

`original_tau_smc_audit.py` 做两种独立有限核验：

1. 直接枚举全部 (N!) 均匀排列；
2. 枚举每个受限提议树，逐路径精确验证 (q_\tau W_\tau=1/N!\)。

得到 (A_2(2)=8,Q_2=1/2)，及

\[
A_2(3)=143808,\qquad Q_3=1100/2247.
\]

所有非零 τ 分层的精确 (Z_\tau) 与未除以 (n^2) 的平方欧氏距离矩
都逐项相等（后者再除以 (n^2) 即为 (U_\tau)）。128 个 256 粒子
CPU particle-filter 重复给出 (n=3) 的 ratio 诊断约 0.48729，接近但
不替代精确 0.48954；这个有意保留的误差表明小粒子比值不能被包装为证明。

## GPU 实现设计

每粒子存 (n^2\) bit revealed mask、两根、alive 标志和 nonce；每个 rank
GPU 计算 allowed count 与可选顶点，主机或 GPU prefix-sum 根据 (G_r) 做
满足边缘条件的 resampling。child nonce 必须由新的计数器键产生。每个
τ 分层可独立 dispatch；优先把 GPU 做成同一 τ 的大量粒子，避免在
不同 τ 上把单个粒子池混合而漏掉 (1/\pi_\tau) 修正。先以 n=2/3 的
exact strata 审计，再报告 n=8,12,16 等的未归一化 (Z_\tau,U_\tau) 和
退化诊断；未经这些步骤的频率不作为原条件概率。
