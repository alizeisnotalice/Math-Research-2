# r021 随机第二出生秩的重要抽样

令 (N=n^2)、(	au\in\{2,\ldots,N\})、(L=N-1)、
(K=\min(6,N-2))、(D_0=2^K(N-K-1))。令 (j=	au-1)，并令
`gnum` 为 (2^{K-j}(N-K-1))（若 (j\le K)）或 1（其余尾部）。提议分布为

\[
 q_\tau=\frac{3L\,\mathrm{gnum}_\tau+D_0}{4D_0L}.
\]

它具有全支撑。实现先以严格的 UInt64 拒绝法在
`0,...,qden-1` 均匀取整数，再按其 CDF 选层；不使用有偏的 `% qden`。

若固定层 SMC 的未归一化估计分别无偏估计 (z_\tau,d_\tau)，则单次随机层
输出 (X_Z=\widehat z_\tau/q_\tau)、(X_D=\widehat d_\tau/q_\tau) 无偏估计
\(\sum z_\tau,\sum d_\tau\)。条件根距的报告量只能是
\(\bar X_D/(N\bar X_Z)\)，不是每次 (D/Z) 的均值，也不是无偏估计。

对真值 (Q)，delta 方差的层混合分解为
\[
\sum_\tau\frac{\operatorname{Var}(\widehat d_\tau-NQ\widehat z_\tau)}{q_\tau}
+\sum_\tau\frac{(d_\tau-NQz_\tau)^2}{q_\tau}.
\]
归档汇总器跨 `globalrep` 的样本方差同时包含这两项，但仅为数值诊断而非认证置信区间。

层选择和轨迹使用独立 domain。冻结 GPU CLI 的种子是
`0xA022... + domain + n*10^6 + tau*10^4 + localrep`；调度器将每个独立的
SplitMix desired seed 反解为 CLI domain，并在汇总时重放检查。
