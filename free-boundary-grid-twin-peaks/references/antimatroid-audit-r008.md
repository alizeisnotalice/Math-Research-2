# r008：图搜索 antimatroid 与峰数多项式的文献/结构核验

**任务。** 令 \(G\) 为连通图，\(S\) 为非空的指定最高标号顶点集。\(h_G(S)\) 是在先固定 \(S\) 为最高 \(|S|\) 个标号后，余下顶点按降标号加入、每一步均与较早顶点相邻的顺序数。对 \(|S|=2\) 且 \(S\) 非邻，\(h_G(S)\) 正是“全局最高、次高恰为两峰”的扇区权重。本记录只审计能否由现有 greedoid/antimatroid 理论给出方格所需的完成核估计或空间比较；不把有限计算当作一般结论。

## 已核实的结构同一性（严格）

将 \(S\) 收缩为一个根 \(r\)（平行边/环不影响“下一顶点有早邻居”的可行性）。在 \(G/S\) 上，令
\[
 \mathcal F_r=\{X\subseteq V(G/S)\setminus\{r\}:(G/S)[\{r\}\cup X]\text{ 连通}\}.
\]
它是 rooted vertex-search antimatroid：空集可达，且若 \(X,Y\) 各与根连通，则 \(X\cup Y\) 仍与根连通。一个 basic word 恰是从根开始、每步加入与已访问集相邻顶点的全排列。因此
\[
 h_G(S)=\#\{\text{basic words of }\mathcal F_r\}.
\]
这里 \(S\) 原先可不连通；收缩只是一种计数编码，不能把原图曼哈顿距离改成收缩图距离。

这与 Björner--Ziegler 给出的 rooted graph vertex-search greedoid 的定义及其 basic-word/search-ordering 对应完全一致。它只提供**准确的语言/可行集模型**，不是完成数的比较不等式。

## 已查来源及可用范围

|来源（实际阅读范围）|可以严格取用|本问题中不能推出|
|---|---|---|
|A. Björner, G. M. Ziegler, *Introduction to Greedoids*, 1992，PDF 的 antimatroid 与 vertex-search 段落（全文相关段）|有限 antimatroid 当且仅当 accessible 且 union-closed；rooted 图中 \(X\cup\{r\}\) 连通的集合形成 vertex-search antimatroid，basic words 是搜索顺序。|basic-word 总数随根位置、根集几何或边界长度的单调性；任何 \(p(S)\) 上界。|
|D. Eppstein, *Antimatroids and Balanced Pairs*, arXiv:1302.5967（摘要、node-search 定义、Theorem 19 段）|对**distance-hereditary** 图的 node-search antimatroid，存在某一对元素，其在均匀 basic word 中的先后概率落在 \([1/3,2/3]\)。|此为“存在一对”的平衡结论，不是固定前缀 \(S\) 的 completion probability，也不是 \(h_G(\{u,v\})\) 的平均/距离比较。\(P_n\square P_n\;(n\ge3)\) 的外边界是诱导 \(C_{4(n-1)}\)（长度至少 8），故不是 distance-hereditary；该定理不适用。|
|S. H. Chan, I. Pak, *Log-concave poset inequalities*, arXiv:2110.10740（摘要、Introduction/Theorem 1.31 与 WeakLoc 讨论）|对 interval greedoids 的**满足文中 k-admissibility 的权重**，给出加权可行 word 层数的对数凹型不等式。|其 \(L_q(k)\) 是可行 word 的长度层统计，而非固定 \(S\) 后的 basic-word completion count；权重须先验验证 k-admissible。文中还明确 WeakLoc 并非所有 antimatroid 成立。因此不能直接套到均匀 basic words，不能导出 \(p(S)\le4/|\partial S|\)。|
|V. Cardinal et al., *Optimal antimatroid sorting*, ESA 2025 / arXiv:2507.13994（摘要、vertex-search 应用段）|vertex-search orders 可作为受限排序对象并可优化比较策略。|优化查询/排序算法没有给均匀 basic-word 的计数或 completion-probability 比较。|
|A. Diaz-Lopez, L. Everham, P. E. Harris, E. Insko, V. Marcantonio, M. Omar, *Peaks on Graphs*, arXiv:1708.08493（摘要、定义和算法主张）|固定峰集的图标号可递归构造；定义与本问题（方格所有顶点度至少 2）兼容。|该文未给任意图或方格“按峰数”生成多项式的全实根、普通对数凹或能控制 \(A_2/A_1\) 的渐近定理。|

## 对候选边界不等式的结论

令 \(p(S)\) 表示给定已完成连接核 \(S\) 后的“两步/下一步完成”条件概率（采用主项目的精确定义）。目前查到的 antimatroid 文献没有一个定理的结论或假设能化为
\[
 p(S)\le \frac4{|\partial S|}
\]
或关于 \(h_G(\{u,v\})\) 随两点位置/距离的比较。一般 bounded-degree broom 图已有反例，说明单凭 antimatroid 公理或有界度不够；研究其他加权或限制形态版本时，仍必须核验额外几何结构。该一般候选随后已由本轮n8反例严格否定；它不能作为已知graph-search不等式的推论。

同理，边增加会增加固定根集的可行 basic words（每个原搜索序仍可行），但不同 \(\{u,v\}\) 的收缩图有不同顶点识别；该单调性不能比较不同两点对的 \(h_G\)。小网格已见“更远的 pair 权重更小”的实例，故不存在未经附加条件的逐对“距离越远权重越大”原则。

## 峰数单变量多项式：一个严格的格点反例

设 \(F_G(z)=\sum_k A_kz^k\)，其中 \(A_k\) 是恰有 \(k\) 个（邻域严格）局部极大的标号数。对 \(P_3\square P_3\)，逐一枚举 \(9!\) 个排列得到
\[
 F_3(z)=29568z+143808z^2+160704z^3+24000z^4+4800z^5.
\]
令 \(Q(z)=F_3(z)/z=\sum_{j=0}^4q_jz^j\)。普通对数凹在 \(j=3\) 已失败：
\[
 q_3^2=24000^2=576000000
 <160704\cdot4800=771379200=q_2q_4.
\]
更强地，若 \(Q\) 有全实根，由正系数知根皆非正，Newton 不等式要求 \(q_j/\binom4j\) 对数凹；但
\[
 (q_3/4)^2=6000^2=36000000
 <(q_2/6)q_4=26784\cdot4800=128563200.
\]
故 \(Q\)，从而 \(F_3\)，并非全实根。脚本 `verify_peak_polynomial_counterexample.py` 用与 independent-set 递推无关的 \(9!\) 穷举核验这一直方图和两条严格整数不等式。

这排除了“对任意方格峰数多项式套用 real-rootedness/ULC/ordinary log-concavity，进而控制 \(A_2/A_1\)”这一路。它不排除别的、带权或截断后的结构性质；截至本轮，未查到适用于本 \(F_G\) 的此类性质。

## 对主问题的实际增量

* 可信的转换是 \(h_G(S)\leftrightarrow\) 收缩根图的 vertex-search antimatroid basic-word 数，允许采用二维割线、反射或注入直接研究计数。
* 所检索的平衡/对数凹/排序结果均不提供 top-2 扇区远峰的常数下界，也不提供正向 \(p(S)\) 估计。
* 生成多项式的通用实根/对数凹路线有已认证的 \(3\times3\) 反例；今后不应以此为候选主引理。

## 可复现性

* 检索日期：2026-11-15；检索式包括 `vertex search antimatroid basic words rooted connected graph`, `antimatroids balanced pairs node search`, `log-concave poset inequalities interval greedoids`, `Peaks on Graphs`。
* 本轮严格有限核验：`python3 -B verify_peak_polynomial_counterexample.py`，范围为 \(3\times3\) 的全体 \(9!\) 标号。
* 本轮不声称检索穷尽，也未从上述文献推出任何方格渐近结论。

## 访问链接

* Björner--Ziegler PDF: <https://page.math.tu-berlin.de/~felsner/Lehre/SemMatS/Literatur/BjoernerZiegler%3AGreedoids.pdf>
* Eppstein: <https://arxiv.org/abs/1302.5967>
* Chan--Pak: <https://arxiv.org/abs/2110.10740>
* Cardinal et al.: <https://arxiv.org/abs/2507.13994>
* Diaz-Lopez et al.: <https://arxiv.org/abs/1708.08493>
