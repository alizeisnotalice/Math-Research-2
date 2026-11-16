# 原始双峰条件律的顺序重要性抽样设计

## 对象与带权恒等式

把一个双射标号按从高到低写作顶点序列
\(X_1,\ldots,X_N\)，令 \(S_t=\{X_1,\ldots,X_t\}\)。
在时刻 \(t\) 出生，当且仅当 \(X_t\notin\partial S_{t-1}\)，即它没有
更高标号邻点；出生点恰好是最终的严格峰。

考虑下列提议核 \(q\)：初始从 \(V\) 均匀取点；已经有一峰时从
全部未揭示点均匀取；已经有两峰时从当前边界 \(\partial S\) 均匀取。
令 \(a_t\) 是该时刻实际可取集合大小，定义
\[
 W=\prod_{t=1}^N a_t.
\]
对每一个有至多两峰的完整下降序列 \(x\)，提议概率正好为
\(q(x)=W(x)^{-1}\)；反之每条提议路径都恰有一或两峰。因此对
\(k=1,2\)，以及任何峰对函数 \(f\)，严格有
\[
 A_k=\mathbb E_q[W\mathbf1_{\{K=k\}}],\qquad
 \sum_{K=2}f=\mathbb E_q[W\mathbf1_{\{K=2\}}f]. \tag{1}
\]
这里没有近似：它只是有限路径树上的重要性采样恒等式。既有记号的
\(Q_n\) 取 \(f=\|X-Y\|_2^2\)，故得到 \(A_2Q_nn^2\)。原题中的
距离 \(\rho\) 仍是方格图的曼哈顿距离；\(\rho^2\le2\|X-Y\|_2^2\)
正是收尾尾界带系数2的原因。审计同时记录 \(\rho^2\) 总和。

`original_law_tree_audit.py` 逐叶用 `Fraction` 枚举该提议树并验证
\(qW=1\)，没有调用此前独立集或包含排除程序。实际运行的结果为：

| n | \(A_1\) | \(A_2\) | \(Q_n\) |
|---:|---:|---:|---:|
|2|16|8|1/2|
|3|29568|143808|1100/2247|

因此该直接原模型提议的核与既有完整标号结果吻合。输出还按第二出生的
下降秩保存精确计数和距离平方和。

## 用于大尺寸的分层版本

原提议在第二峰很早或很晚时权重可能严重集中。以第二出生下降秩
\(\tau\in\{2,\ldots,N\}\) 分层可避免把罕有层当作频率估计：

1. 取第一个点；在位置 \(2,\ldots,\tau-1\) 从边界均匀取，强制尚无第二峰；
2. 位置 \(\tau\) 从未揭示的非边界点均匀取，强制第二峰；
3. 以后始终从边界均匀取，强制不再产生峰。

空的强制集合表示该层没有支持。对有支持的一条路径，仍以每一步实际
强制集合大小的乘积作 \(W_\tau\)。它逐点满足
\(q_\tau(x)W_\tau(x)=1\)，故
\[
 A_{2,\tau}=\mathbb E_{q_\tau}W_\tau,\qquad
 D_{2,\tau}=\mathbb E_{q_\tau}(W_\tau\|X-Y\|_2^2),
\]
并且 \(A_2=\sum_\tau A_{2,\tau}\)、\(D_2=\sum_\tau D_{2,\tau}\)。
某些强制前缀会在 \(\tau\) 前走入无可用点的死端；把它记为零贡献样本，
不能在 GPU 中丢弃后重新按幸存样本归一化。所有实际有支持的完整序列仍逐点
满足上式，故零贡献处理不改变无偏恒等式。
这给出可并行的 GPU 接口：每个 \((\tau,\text{root stratum})\) 批次保存
`sum_W`, `sum_W_euclidean2`, `sum_W_manhattan2`, `sum_W2`、零支持数、种子和命令缓冲区状态。

单峰可另行以“首点任取、以后边界均匀取”的连通序提议计算；其权重
期望恰为 \(A_1\)。这比在同一两峰提议中等待 \(K=1\) 终点更稳定。

一个直接可实现的每线程接口如下；这里 `rank` 从1开始，`born` 由新点
是否接触当前 `S` 决定。

```text
state = (S=0, births=0, peak1, peak2, logW=0, dead=false)
for rank = 1..N:
    allowed = V\S                    if births < 2 else boundary(S)
    choose v uniformly from allowed; logW += log(|allowed|)
    if v not in boundary(S): births += 1; record the new peak and its rank
    S += v
at rank N: emit exp(logW) * 1[births==2],
           exp(logW) * 1[births==2] * euclidean2(peak1,peak2),
           exp(logW) * 1[births==2] * manhattan2(peak1,peak2)
```

分层版本只把 `allowed` 改为：rank小于 \(\tau\) 时 `boundary(S)`，等于
\(\tau\) 时 `unrevealed\boundary(S)`，之后仍为 `boundary(S)`；第一个
rank始终用全部顶点。若 `allowed` 为空，写零贡献记录并停止该线程。
实现必须保存的最小输出为每个层的 `M, dead, sumW, sumW_e2, sumW_rho2,
sumW2, max_logW, max_normalized_weight, seeds, commandBuffer_status`。为
校准，先对n=2、n=3运行相同随机流的单轨迹CPU/GPU比较；不是只比较最终
浮点均值。

## 统计边界

对每个固定层，样本均值是理想随机数下上述分子、分母的无偏估计；各层
独立运行并相加时仍无偏。\(\widehat Q=\widehat D_2/(N\widehat A_2)\)
是自归一化比值，有限样本通常有偏，且随机流、有限精度、GPU 实现错误和
SMC 重采样谱系会改变其理论条件。故此算法的数值输出只能用于发现构造或
诊断权重退化，不能被报告为 \(\mathbb P(\cdot\mid K=2)\) 的认证估计。

若加入粒子重采样，每次 resample 后必须为每个子代刷新其未来随机状态；
复制父粒子的 RNG 会令相同父系后代确定性地共享整个未来，不可把粒子数
解释为有效独立样本数。保存祖先数、最大归一化权重和跨重复的 Z 权重诊断，
但这些仍非严格误差证书。

实用的失败诊断是：层内死端比例、\(\sum W\)^2/\sum W^2\) 的权重ESS、
最大归一化权重、跨独立重复的层和、以及重采样后的不同祖先数。某层只有
少数非零样本、ESS接近1、最大权重接近1，或跨重复的Z加权和不稳定时，
该层不能外推到更大n；应改变分层/提议，而非增加同一退化配置的宣传性样本。

## 实际复现

```sh
python3 -B original_law_tree_audit.py 2
python3 -B original_law_tree_audit.py 3
```

Python 3.15 标准库即可；没有后台进程或随机种子。GPU 实现前应先让单条
CPU/GPU轨迹在相同整数随机流下匹配 \(a_t,W,\tau,\rho^2\)，并检查每个
command buffer 的完成状态。
