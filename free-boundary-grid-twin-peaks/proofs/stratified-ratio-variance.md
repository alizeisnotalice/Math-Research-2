# 第二出生秩分层的比值误差账本

令理想独立抽样的第 τ 层有 M_τ 个样本，U=W_τ，V=W_τ f，f=欧氏根距平方/n²。死路径 U=V=0。记 a_τ=E U，b_τ=E V，A=Σa_τ=A₂，B=Σb_τ=A₂Q。各层样本独立时，T=Σ平均U、S=Σ平均V 有 E T=A、E S=B。若 T>0，比值满足精确代数关系

\[
\widehat Q-Q=(S-QT)/T.
\]

对任意确定 q，严格有

\[
\operatorname{Var}(S-qT)=\sum_\tau\frac{E[W_\tau^2(f-q)^2]-(b_\tau-qa_\tau)^2}{M_\tau}.
\]

因此使用原始二阶矩构成的

\[
\widehat V(q)=\sum_\tau\frac{\sum_iW_i^2f_i^2-2q\sum_iW_i^2f_i+q^2\sum_iW_i^2-(\sum_iW_if_i-q\sum_iW_i)^2/M_\tau}{M_\tau(M_\tau-1)}
\]

在固定 q、每层 M_τ≥2、理想独立随机数条件下是上式方差的无偏估计。将 q 换为同批数据算出的 Qhat，再除 T²，只得到 delta 方法的经验方差诊断；不再声称有限样本严格无偏，也不是置信上界。实际有限 PRNG 也不自动满足独立性假设。

按每条样本实际贡献 U_i/M_τ 定义的全局权重 ESS 为

\[
\frac{T^2}{\sum_\tau M_\tau^{-2}\sum_i U_i^2}.
\]

该量与最小层内 ESS 不同。低质量晚秩层很差的层内 ESS 不必意味着整体 A₂ 估计退化；反之全局权重 ESS 也不能单独控制 f 加权分子的误差。

GPU 若保存 f_raw=欧氏根距平方，应在上式中用中心 n² Qhat，并在最终方差除以 n⁴；不能只归一化一阶矩而漏掉二阶矩。全死层对总和贡献为零，但 M_τ 仍为全部提议数。
