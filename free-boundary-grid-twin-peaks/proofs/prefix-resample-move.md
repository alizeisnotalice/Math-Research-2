# 完整揭示词上的 resample--move 核：有限状态引理

固定 (G=P_n\square P_n)、(N=|V|) 与 (2\leq\tau\leq N+1)。长度为
(r) 的状态是无重复有序词 (x=(x_1,\ldots,x_r))。令 (b_t(x)=1)
当且仅当 (x_t) 与先前集合 (S_{t-1}=\{x_1,\ldots,x_{t-1}\})
没有边；称 (x) 有效，当且仅当

\[
 \{t\le r:b_t(x)=1\}=\{1,\tau\}\cap\{1,\ldots,r\}.
\]

以 ({\cal X}_r^\tau) 记有效词集合，并取未归一化 FK 目标
\(\gamma_r(x)=1_{{\cal X}_r^\tau}(x)\)。这逐位置检查 birth，绝未把
prefix 压缩成无序子集或 canonical word。

## 对称 move

给定 (x\in{\cal X}_r^\tau)，均匀抽取
((i,v)\in\{1,\ldots,r\}\times V)。若 (v=x_j) 已出现在词中，交换
第 (i,j) 位；否则用 (v) 替换第 (i) 位。若所得词无效就拒绝并留在
(x)，有效则接受。记核为 (K_r)。

**引理（对称性）。** 对不同有效词 (x,y)，
\(K_r(x,y)=K_r(y,x)\)。

**证明。** 顶点集合相同的非平凡转移只能交换 (i\ne j)。由
((i,x_j))、((j,x_i)) 两个提议均能从 (x) 到 (y)，总概率
(2/(rN))；逆向恰有相同计数。顶点集合不同的相邻词只能在同一位置
(i) 以删去点和一个新点不同；前向、逆向均只有一个替换提议，概率
(1/(rN))。两种情形不重叠。`swap(i,i)` 只给对角提议，无效候选也只
增加对角元。□

故均匀 (eta_r=gamma_r/|\mathcal X_r^\tau|) 对 (K_r) 可逆，对每个
预定有限 (L\ge0)，(K_r^L) 不变。这里无需也不声称不可约、混合或
prefix 频率等于原条件律。

## FK 校正

从有效父词出发，强制集合为：rank 1 用 V；1<r<tau 或 r>tau 用当前
boundary；r=tau 用未揭示非-boundary 点。记大小为 a_r。本文使用计数势
G_r=a_r：先按父态 a_r 无偏 multinomial resampling，再让每个子代独立从
该父态的允许集合均匀 mutation；空集为零势死亡。
每条允许父子边满足

\[
 \gamma_{r-1}(x)M_r(x,y)G_r(x,y)=1=\gamma_r(y).
\]

因此计数 normalizer Z_N=sum gamma_N=A_(2,tau)。若相对于所有无重复词
记概率 normalizer，则每层势为 a_r/(N-r+1)，即
Z_prob,r=Z_count,r/(N)_r。在 r<N 做上述 resample、mutation 后再执行
预定 L 次 K_r；最后 r=N 只按父态势加权并各自 mutation，不再 resample
或 move，终端函数也以最后势加权。bootstrap FK 的条件期望递推及
K_r^L eta_r=eta_r 给出 E Zhat_N=Z_N；不需要 MCMC 混合。

当 tau<=N-1 或 tau=N+1（一峰约定）时，连通方格图中存活的 N-1 前缀所余唯一点必邻接已揭示集，
所以它在 boundary，终端势为 1、不死亡。可是 tau=N 要求最后点为
non-boundary，强制集合为空，终端全死且 A_(2,N)=0。早期空强制集合同样
必须作为零贡献死亡保留。终端比值 Dhat/Zhat 仍只是有限样本数值诊断。
