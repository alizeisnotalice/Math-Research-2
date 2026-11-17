# R055：\(\Psi\) 的时间加权多尺度路线

## 范围

令 \(G=P_n\square P_n\)、\(N=n^2\)，并令
\(\sigma=(v_1,\ldots,v_N)\) 是一峰降序词。沿 R036/R038，根为 \(r=v_1\)，
对每个 eligible 顶点 \(p\notin\{r\}\cup N(r)\)，记最早邻点秩为 \(l_p\)，
first singleton-child 截止秩为 \(b_p\)，以及
\[
 R_p=b_p-l_p-1,\qquad
 \Psi(\sigma)=\sum_p\frac{l_p-1}{R_p}.
\tag{1}
\]

R054 已严格给出前缀占用恒等式
\(\Psi=\sum_t\sum_{p:l_p>t}1/R_p\)。本轮不重复它；目标是考察它遗漏的部分：
若连通前缀很早支配全图，如何从该前缀自身的形成历史支付额外的对数因子。

## 两个严格的逐词对照

令 \(\Delta=4\)，\(M=N-1-\deg(r)\) 是 eligible 顶点数。将 \(l_p-1\) 从小到大
排为 \(c_1\le\cdots\le c_M\)。每个已出现顶点最多是 \(\Delta\) 个点的最早邻点，
且 eligible 点的 \(l_p\ge2\)，所以
\[
 c_j\ge\left\lceil\frac j\Delta\right\rceil.
\tag{2}
\]
又 \(R_p\le N-l_p=N-1-(l_p-1)\)。函数 \(x/(N-1-x)\) 单调递增，故严格有
\[
 \Psi(\sigma)\ge
 \sum_{j=1}^M
 \frac{\lceil j/\Delta\rceil}{N-1-\lceil j/\Delta\rceil}
 \ge\frac{M(M+\Delta)}{2\Delta(N-1)}.
\tag{3}
\]
对自由方格 \(M\ge N-5\)，这给出 \(\Psi=\Omega(N)\) 的显式点态底界。它只使用
最大度，不会产生 \(\log n\)，因而是多尺度候选必须严格超过的基线。

令 \(\tau\) 为首次支配时间，即最小的 \(t\) 使 \(U_t=\{v_1,\ldots,v_t\}\)
的闭邻域覆盖 \(V\)。所有 eligible \(p\) 都有 \(l_p\le\tau\)，故由同一容量计数
\[
 \boxed{\quad \tau\ge\left\lceil\frac{M}{4}\right\rceil
 \ge\left\lceil\frac{N-5}{4}\right\rceil.\quad}
\tag{4}
\]
这严格排除了 \(o(N)\) 时间的支配；但 \(\tau=\Theta(N)\) 仍可能使 R054 的
未支配集在后段完全消失。

## 一个可复用的链调度引理

称一组 vertex-disjoint 序列
\[
 x_{i,0},x_{i,1},\ldots,x_{i,L_i}\qquad(1\le i\le h)
\]
为 clean chain family，若相邻项在方格中相邻、其秩严格递增，且所有内部点
\(p=x_{i,j}\)（\(1\le j<L_i\)）eligible，
其最早邻点恰为 \(x_{i,j-1}\)，且 \(x_{i,j+1}\) 在加入时以 \(p\) 为唯一早邻点。
令 \(M_{\mathcal C}=\sum_i(L_i-1)\)。

写 \(a_{i,j}=\operatorname{rank}(x_{i,j})\)。则
\[
 l_{x_{i,j}}=a_{i,j-1},\qquad
 R_{x_{i,j}}\le a_{i,j+1}-a_{i,j-1}-1.
\tag{5}
\]
每一条链的右端分母和小于 \(2N\)，所以总和小于 \(2hN\)。又
\(a_{i,j-1}-1\) 是 \(M_{\mathcal C}\) 个不同的正整数，故
\[
 \sum_{i,j}\sqrt{l_{x_{i,j}}-1}
 \ge\sum_{u=1}^{M_{\mathcal C}}\sqrt u
 \ge\frac23 M_{\mathcal C}^{3/2}.
\]
Engel--Cauchy 因而给出严格的调度下界
\[
 \boxed{\quad
 \sum_{p\in\mathcal C}\frac{l_p-1}{R_p}
 \ge\frac{2M_{\mathcal C}^3}{9hN}.
\quad} \tag{6}
\]
这不是 comb 族的重证；它是任何可抽取 clean chain family 都可调用的有限图引理。

## 一个必须避免的假加强

不能把 R036 的 \(\Gamma=\sum1/R_p\) 全量质量乘以面积来声称
\(\Psi\ge cN\Gamma\)。逐行词已经严格反驳这种形式：取
\(\operatorname{rank}(x,y)=yn+x+1\)。顶行的 eligible 点
\((x,0)\), \(2\le x\le n-2\)，满足 \(l=x,b=x+2,R=1\)，故
\(\Gamma\ge n-3\)，但它们对 \(\Psi\) 的总贡献只有 \(O(N)\)。
首列和内部点的直接分母计算给出
\[
 \Psi\le N(1+H_{n-1})+1,
\]
从而 \(\Psi/(N\Gamma)\to0\)。因此有用的多尺度论证必须选取受秩控制的
\(\Gamma\) 子质量，不能使用全量 \(\Gamma\)。

## 明确候选：早支配的链打包原则

下面是一个足以推出所需点态结论的二维候选，而不是已证命题。

**候选 P（时间加权链打包）。** 存在绝对 \(c_0>0\)，使每个 \(n\ge n_0\) 的
一峰词都含有限个两两 vertex-disjoint 的 clean chain families
\(\mathcal C_\alpha\)，并满足
\[
 \sum_\alpha\frac{M_{\mathcal C_\alpha}^3}{h_\alpha}
 \ge c_0N^2\log n.
\tag{P}
\]
由 (6) 的逐族相加，P 会立即给出
\[
 \Psi(\sigma)\ge\frac{2c_0}{9}N\log n.
\tag{7}
\]

该候选的含义是：当前缀在某些尺度已经形成早支配骨架时，不能只看静态的
dominating set；必须从它的**形成顺序**提取 predecessor--singleton-child 链，
并以 (6) 的秩--间隔账本累积。它与 R054 的 \(R_p\) 依赖完整词这一警告相容，
因为 (5) 只把实际的 singleton child 当作上界证据，不对残余图重算 \(R\)。

## 第一处未证的二维步骤

目前没有从任意连通、早支配方格前缀抽取 P 中链包的定理。两个直接替代都不足：

1. 前缀占用只在未支配点上计 \(1/R_p\)。一旦前缀支配全图，该部分精确为零，
   不能再直接套 R036 的根壳层 \(\Gamma\) 下界。
2. 静态连通支配或边界长度不控制 (5) 中的 first singleton-child 时间。
   一个新点可以有两个早邻点，从而不为任何已有点提供这样的 child。

这个缺口不是 \(R_p\) 的索引问题，而是从二维早支配的几何形态推出足量 clean
chain score 的结构提取问题。证明 P 或构造 \(\Psi=o(N\log n)\) 的一峰词，才会
改变当前路线的状态。

## 一个严格的反诊断：不能要求晚时刻仍有未支配质量

为排除一个过强但自然的候选，取任意奇数 \(n\ge5\)，先沿第 \(0\) 列加入全部顶点，
随后加入每条偶数行的其余点，最后才加入奇数行的其余点。这是合法一峰词。前两段
已经组成连通支配前缀，大小至多
\[
 T=n+\frac{n+1}{2}(n-1)=\frac{n^2+2n-1}{2}.
\tag{8}
\]
每个剩余的奇数行点同时邻接两条已出现的偶数行，故该时刻以后不存在未支配点。
也就是说 R054 的 \(g(t)\) 在 \(t\ge T\) 恒为零（实际首次支配可能更早）。
这个例子不是
\(\Psi\) 小的反例；它只否定“对线性比例的所有晚时刻仍有
\(g(t)\gtrsim\log n\)”之类的候选。因此 P 若为真，贡献必须来自早支配骨架自身，
正是 (6) 试图记录的历史成本。

## 有限审计

脚本 <code>audit_domination_time.py</code> 用整数/有理数检查一个既有 \(8\times8\)
低 \(\Psi\) 搜索词及上述 \(9\times9\) comb 词。它核验一峰条件、所有 \(R_p>0\)、
首次支配时间，以及所有 eligible \(l_p\le\tau\)。这只是定义与边界的低内存审计；
(2)--(6) 的严格性来自上述计数和不等式，P 仍是待证候选。

## 后续状态修订（R058）

候选 P 已被逐行一峰词**严格否定**，不应再作为待证路线。该词中具有
sole-parent child 的 eligible 父点恰有 \(2(n-3)\) 个，因此任意彼此顶点不交的
clean-chain families 都满足
\[
 \sum_\alpha M_{\mathcal C_\alpha}^3/h_\alpha\le8(n-3)^3
 =o(N^2\log n).
\]
详见独立材料
<code>/Volumes/data/CodexTemp/T2/a022-chain-packing-obstruction-r058/ROW_EXCLUDES_CHAIN_PACKING.md</code>。
局部链调度引理 (6) 仍然正确；失效的是以 clean chain 覆盖所有 \(\Psi\) 代价的
全称提取步骤。后续必须另行计入 \(b_p=N+1\) 的无 sole-child 项。
