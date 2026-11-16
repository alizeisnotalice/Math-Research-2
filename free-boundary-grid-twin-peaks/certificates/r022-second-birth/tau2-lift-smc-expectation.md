# r022：\(\tau=2\) 提升的理想 SMC 一阶无偏性

## 命题及范围

令 \(G=P_n\square P_n\), \(N=n^2\)，且 \(n\ge2\)。令
\(\mathcal X_k\) 为所有长度 \(k\) 的无重复有效前缀：它们的出生位置正是
\(\{1,2\}\cap\{1,\ldots,k\}\)。对定义在 \(\mathcal X_k\) 的函数 \(f\)，记

\[
 \mu_k(f)={1\over(N)_k}\sum_{x\in\mathcal X_k}f(x).
 \tag{S1}
\]

这是相对于均匀全排列参考律的**未归一化**前缀测度，而非 \(\mathcal X_k\)
上的概率律。

若 \(x\in\mathcal X_{k-1}\)，令 \(A_k(x)\) 是使延长前缀仍有效的允许点集，
\(a_k(x)=|A_k(x)|\)，并令 \(Q_k\) 均匀选取该允许点。采用概率归一化的势

\[
 G_k(x)={a_k(x)\over N-k+1}; \tag{S2}
\]

其中 \(a_k=0\) 的状态给出零势。每一个有效长度 \(k\) 词只有一个长度
\(k-1\) 父前缀，所以逐和式给出

\[
 \mu_{k-1}(G_kQ_kf)=\mu_k(f). \tag{S3}
\]

在每一个非末层，允许在 mutation 后施加一个核 \(M_k\)，只要它留在
\(\mathcal X_k\) 内且保持其均匀测度，即

\[
 \mu_kM_k=\mu_k. \tag{S4}
\]

文献目录中的 `prefix-resample-move.md` 所给交换/替换拒绝核是对称的，故满足
(S4)；它无需不可约或混合结论。

考虑有 \(P\) 个粒子的 bootstrap 算法。给定父粒子和势，设 offspring 数
\(N_i\) 满足唯一所需的一阶条件

\[
 \mathbb E[N_i\mid\mathcal F]=P{G_k(X_i)\over\sum_jG_k(X_j)}. \tag{S5}
\]

随后每个 child 按 \(Q_k\) mutation，再可按预定次数施加 \(M_k\)。例如以理想的
独立均匀随机整数驱动的整数分层重采样满足 (S5)，见
`integer-stratified-resampling.md`；多项式重采样也满足它。

令 \(\widehat\gamma_k=\widehat Z_k\eta_k^P\) 表示第 \(k\) 层的未归一化粒子
测度。则在**真独立随机字和精确算术**的理想算法中，对每个有界 \(f\)，

\[
 \mathbb E\widehat\gamma_k(f)=\mu_k(f),\qquad 1\le k<N. \tag{S6}
\]

证明是归纳。给定第 \(k-1\) 层，(S5) 和后续 conditionally-correct mutation 给出

\[
 \mathbb E[\widehat\gamma_k(f)\mid\mathcal F_{k-1}]
 =\widehat\gamma_{k-1}(G_kQ_kM_kf).
\]

取期望、使用归纳假设，再用 (S3)--(S4)，右侧为 \(\mu_k(f)\)。这只使用
offspring 的条件一阶矩，不要求 child 有相同边缘分布或粒子彼此独立。

## 末步和第二出生提升

长度 \(N-1\) 的有效方格前缀尚余唯一顶点；它有已揭示邻点，故
\(a_N=1\)、\(G_N=1\)，且 \(Q_N\) 是确定性补上末点。末步不应再 resample。
若 \(H\) 是完整词统计量，并把确定性补点后的函数仍记作 \(H\)，末端估计量应写为

\[
 \widehat T_H=\widehat Z_{N-1}\eta_{N-1}^P(G_NH)
 \quad\text{（等价地 }\widehat\gamma_{N-1}(G_NQ_NH)\text{）}. \tag{S7}
\]

由 (S6) 和 (S3)，

\[
 \mathbb E\widehat T_H=\mu_N(H)
 ={1\over N!}\sum_{\sigma\in\mathcal S_{\tau=2}}H(\sigma). \tag{S8}
\]

另一方面，稳定前移的严格插槽双射给出

\[
 \sum_{\omega:K(\omega)=2}h(\omega)
 =\sum_{\sigma\in\mathcal S_{\tau=2}}L(\sigma)h(\sigma). \tag{S9}
\]

故在 (S7) 中分别取

\[
 H=L,\quad L\|x-y\|_2^2,\quad Lf,\quad Lg
\]

时，\(N!\widehat T_H\) 分别是全 \(K=2\) 计数 \(A_2\)、欧氏平方距离总量
\(D\)、以及 r018 的根矩总量 \(F,G\) 的一阶无偏估计。这里
\(f=\|t_x\|^2+\|t_y\|^2\)、\(g=t_x\cdot t_y\)，\(t_v\) 是 doubled-centred
坐标。也可令 \(h\) 是任意固定曼哈顿距离层的指示函数，以无偏地估计该层的
未归一化质量。

这为以单个 \(\tau=2\) 目标层估计全原 \(K=2\) 未归一化总量提供严格的理想
算法命题；它消除了必须逐个运行 \(\tau=2,\ldots,N\) 层的必要性。

## 实现界限

上述结论不对有限样本 ratio 作无偏性声明：
\(\widehat T_{L\|x-y\|^2}/\widehat T_L\) 或由 \(F,G\) 形成的 ratio 通常有有限
样本偏差，即使基础未归一化估计各自无偏。

真实运行中的有限 PRNG 不是“真独立随机字”，而 `logFloat` / IEEE 浮点的 normalizer
和加总也不等于精确算术。因此它们只能作为此理想命题的数值近似，不能写成程序的
严格无偏性证明。预定 resample--move 的对称性只证明目标不变；本命题并不赋予它
混合率、不可约性或有效样本数的严格置信区间。
