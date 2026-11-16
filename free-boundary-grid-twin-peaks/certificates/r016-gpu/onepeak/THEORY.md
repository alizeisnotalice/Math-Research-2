# 一峰 normalizer 的 FK 身份

令 \(x_1,\ldots,x_N\) 是方格 \(P_n\square P_n\) 的降序标号顶点次序，
\(S_r=\{x_1,\ldots,x_r\}\)，并令

\[
B(S)=N(S)\setminus S.
\]

严格局部峰正是出生点（与已有前缀无邻边的点）。因此恰有一个峰，当且仅当
\(x_1\) 任取而对每一个 \(r\ge2\)，\(x_r\in B(S_{r-1})\)；等价地每个
非空真前缀均连通地增长。把旧两峰记号形式化地取为 \(\tau=N+1\)，则 kernel
中的 `rank==tau ? !touch : touch` 在所有实际秩 \(r\le N\) 都选择 `touch`，
故 MSL 本体无需改变。

固定这种受限提议 \(q_1\)：第一步均匀选 \(V\)，其后在 \(B(S)\) 均匀选点。
令

\[
G_r(S_{r-1})=\frac{|A_r(S_{r-1})|}{N-r+1},\qquad
A_1=V,quad A_r=B(S_{r-1})\ (r>1).
\]

对任何完整次序 \(x\)，逐步相消给出

\[
q_1(x)\prod_{r=1}^{N}G_r
=\frac{1}{N!}\mathbf1_{\{K_n=1\}}(x).
\]

故

\[
Z_1:=\mathbb E_{q_1}\prod_rG_r
=\mathbb P(K_n=1)=\frac{A_1(n)}{N!}.
\]

理想随机数下，若每一步重采样的每个父粒子 offspring 条件期望等于其势
归一化权重乘 \(P\)，再配合独立正确 mutation 核，则标准 Feynman--Kac 条件
期望归纳给出未归一化 \(\widehat Z_1\) 的一阶无偏性。这里采用整数分层目标；
它满足该 offspring 条件期望。实际 PRNG、有限 \(P,R\) 仍只是数值实现。

本程序只输出 \(Z\)。它保留 kernel 中的 `root2`/终端 moment 路径以保持
Metal 核逐字节不变，但 host 不输出、也不解释这些无关的一峰量。rank-to-rank
parent 随机流使用 domain 1；选定 parent 后的顶点 mutation 使用独立 domain 2，
不会复制父粒子的未来 RNG 状态。

该身份只校准一峰 normalizer；它本身不估计远峰概率、也不证明任何 \(n\to\infty\)
结论。用于 \(R_n=A_2(n)/A_1(n)\) 时还需要独立、可比的 \(A_2\) 估计与误差控制。
