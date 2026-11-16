# 整数分层重采样：精确边缘与一步条件方差

状态：有限离散概率的严格推导；应用于原条件律SMC时，理想独立随机数是数学假设，有限PRNG实现仍是数值证据。

固定P个父状态X_i及非负整数权重a_i，T=Σa_i>0，累计整数C_i=Σ_{k≤i}a_k。对于整数u∈{0,…,T−1}，令I(u)是第一个满足u<C_i的父索引。零权重父不会被选中。

对每个child j=0,…,P−1，独立取t_j均匀于{0,…,T−1}，定义

\[
u_j=\left\lfloor\frac{jT+t_j}{P}\right\rfloor,
\qquad \text{parent}_j=I(u_j).
\]

## 1. 精确边缘总计

所有有序对(j,t)一一对应整数k=jT+t∈{0,…,PT−1}。映射floor(k/P)把每个u∈{0,…,T−1}恰取P次。因此对任意父i，其offspring数N_i满足

\[
\mathbb E[N_i\mid X,a]=\sum_j\Pr(I(u_j)=i)=P a_i/T.
\]

单个child的父分布不必等于a/T；要求的是以上总计。各child使用独立t_j时，它们条件独立，但不具相同分布。SMC未归一化无偏性只需这个offspring条件和正确mutation核，不要求子粒子同分布。

若所有a_i相同且为正，则T=P a，u_j=j a+floor(t_j/P)∈{j a,…,(j+1)a−1}，故parent_j=j恰成立。此时重采样本身不丢掉父状态；这比多项式重采样在同权情况下无谓复制和删除父状态更适合保留早期根信息。

## 2. 一步条件方差比较

允许选父后独立按mutation核M_i抽取Y_j，观测任意平方可积函数f(Y_j)。考虑K均匀于{0,…,PT−1}，令J=floor(K/T)，U=floor(K/P)。J均匀于P个层，U均匀于T个整数，故从I(U)再mutation得到的Y正是一个多项式重采样child的分布。

分层算法中各child的观测分布恰为f(Y)|J=j，且条件独立。于是

\[
\operatorname{Var}\left(\frac1P\sum_j f(Y_j)\mid X,a\right)
=\frac1P\mathbb E[\operatorname{Var}(f(Y)\mid J,X,a)\mid X,a]
\le \frac1P\operatorname{Var}(f(Y)\mid X,a).
\]

右侧恰为P个独立多项式child平均的条件方差。这是对固定当前父状态的一步比较；没有把它外推为所有多步SMC终端比值的无条件方差定理。

若改成所有child共用一个t（systematic重采样），offspring期望仍正确，但独立性步骤不成立，不能照搬此方差比较。

## 3. 对原条件律的应用和实现边界

在固定第二出生秩τ的SMC中，a_i是父前缀可选顶点数，势为a_i/(N−r+1)。组间CDF和组内CDF均可用精确整数；本变换只改变child的整数目标，不改变势、mutation或原条件事件。T=0时整次重复权重为零，必须保留。末步仍不重采样。

64位实现需保证PT≤2^64，因为最大中间整数为PT−1。常用P=32768,N≤144时T≤PN，故PT≤P²N<2^38，远低于界限。应显式检查参数，不能仅凭常用配置推断任意输入安全。拒绝法生成t_j不产生取模偏差，其结论仍以源64位字真正均匀为假设。

`audit_integer_stratified.py` 用Fraction穷举P≤6、a_i∈{0,1,2,3}的全部非零权重向量，核对整数边缘总计、同权无父丢失及一个确定父观测的一步方差式。计算是实现核验；一般结论由上述有限求和和全方差公式证明。

本轮r012冻结SMC使用多项式选父，尚未采用本引理的分层整数目标。一步方差比较不归因于该运行；后续实现还须独立审计，并实际检查整数范围和独立child随机化。
