# r015：递归 H-骨架树的 hook 下界与固定密度 level boundary

设 (G_n=P_n\square P_n)，(N=n^2)。本轮构造的是一个**确定的**有根生成树，
只把其 heap 标号当作原方格单峰标号的子族；没有把树模型或均匀边界生长误作原条件律。

## 1. 明确的角根递归树

对一个坐标矩形 (Q=[x_0,x_0+w-1]\times[y_0,y_0+h-1])，指定其左上角为根。
若 (w\le2) 或 (h\le2)，取从左上角出发的逐行蛇形 Hamilton 路径为树。

若 (w,h\ge3)，令 (c=\lfloor w/2\rfloor,d=\lfloor h/2\rfloor)。先取 H-骨架：
左列 (x=x_0)、中央列 (x=x_0+c) 和中央行 (y=y_0+d)。**只选择**两列各自的
纵边及中央行的横边；不取这些点集之间的其余方格横边（尤其 (w=3) 时两列相邻，
若取诱导子图会有环）。这是一棵树，顶点数为

\[
 |H_Q|=w+2h-2. \tag{1}
\]

其余顶点至多分成四个矩形：

\[
\begin{array}{c|c|c}
\text{子块}&\text{尺寸}&\text{子根及其到 H 的附着边}\\ \hline
\mathrm{TL}&(c-1)\times d &(x_0+1,y_0)\sim(x_0,y_0)\\
\mathrm{BL}&(c-1)\times(h-d-1)&(x_0+1,y_0+d+1)\sim(x_0,y_0+d+1)\\
\mathrm{TR}&(w-c-1)\times d &(x_0+c+1,y_0)\sim(x_0+c,y_0)\\
\mathrm{BR}&(w-c-1)\times(h-d-1)&(x_0+c+1,y_0+d+1)\sim(x_0+c,y_0+d+1).
\end{array} \tag{2}
\]

空子块省略；对其余子块递归。每个子块只用一条附着边接到 H，故归纳得到一棵
覆盖 (Q) 的树，根保持为指定左上角。起始 (Q=[0,n-1]^2) 的树记为 (T_n)。

## 2. 全 (n) hook 熵界

令 (s_v) 是以左上角根定向后 (v) 的后代子树大小。以下是完全确定性的较强界：

\[
 \sum_{v\in T_n}\log s_v\le17N. \tag{3}
\]

**证明。** 对 (2) 中每个非空子块，直接由四个尺寸公式有

\[
 w_{\rm child}\le\lfloor w/2\rfloor,
 \qquad h_{\rm child}\le\lfloor h/2\rfloor. \tag{4}
\]

这避免了任何 ceiling 或薄条带的额外累计。令一个递归节点的最大边长为 (L)。
若两边都至少 3，H 的点数至多 (3L)，每个 H 点的子树至多 (L^2)，所以它的
hook 对数贡献至多 (6L\log L)。若节点薄，即一边至多 2，则蛇形路径面积至多
(2L)，其 hook 对数至多 (2L\log(2L)\le6L\log L)；(L=1) 时贡献为零。
因此每一个递归节点的贡献一律至多 (6L\log L)。

令 \(J=\lfloor\log_2n\rfloor\)，\(\alpha=n/2^J\in[1,2)\)。由 (4)，第
(j) 层至多有 (4^j) 个节点，且其最大边长至多 \(\lfloor n/2^j\rfloor\)。于是

\[
\begin{aligned}
 \sum_{v\in T_n}\log s_v
 &\le6n\sum_{j=0}^{J}2^j\log(n/2^j)\\
 &\le\frac{12n^2}{\alpha}\log(2\alpha)
 \le24(\log2)N<17N.
\end{aligned} \tag{5}
\]

这里第二行使用 \(\sum_{r\ge0}2^{-r}=\sum_{r\ge0}r2^{-r}=2\)，最后一步用
\(\log2<0.7\)（例如 \(e^{0.7}>2\)）。因而还得到整数形式
\(\prod_vs_v\le2^{24N}\)。这证明 (3)。
此前的 128 常数来自把终止薄块另行粗估的版本；(4) 将它们纳入同一层和，故不是
本轮最终常数。\(\square\)

有根树的 hook-length 公式给出根最高的 tree-heap 标号数

\[
 \frac{N!}{\prod_vs_v}\ge N!e^{-17N}. \tag{6}
\]

每个非根在树中有一个更高父点，这条父边也是原方格边；故它在原方格中不是峰。
因此 (6) 真正给出自由方格单峰计数的全 (n) 下界：

\[
 \boxed{A_1(n)\ge N!e^{-17N}.} \tag{7}
\]

`audit_recursive_hook_tree.py` 实际构造 (1)--(2)，逐项检查树有 (N-1) 条边、
从根连通、整数 hook 商、子块 floor 缩小关系，并精确检验 \(\prod_vs_v\le2^{24N}\)
及 (3) 于 (1\le n\le128)。所有检查通过；最大实测
\(N^{-1}\sum_v\log s_v\) 为约 2.260，远低于安全常数 17。这是计算交叉核验，
(3) 的一般性来自上面的层和证明。

## 3. 固定阈值、正比例高度的双峰粗边界结论

令 (D_j) 为双峰标号的最高 (j) 个标签集，外顶点边界为 \(\partial D_j\)。
高到低揭示中，除首点外恰有一次非边界加入（第二峰）。故对任意 \(s\ge1\) 和
\(1\le r\le N-1\)，至少 (r) 个时刻满足 \(|\partial D_j|\le s\) 的双峰标号数至多

\[
 M_2(s,r)\le N^3\binom Nr s^{r-1}N^{N-r+1}. \tag{8}
\]

计数理由：先选至少r个小边界时刻中的r个，至多binom(N,r)种；再固定第二出生的秩，至多N种。
除第二出生一步之外，选定时刻的下一点至多有s种选择；至少r−1次获得这个限制。
其余各步和首点均以N种粗计，所得界甚至小于(8)；多取N²因子保持统一的宽松形式。
该逐前缀分支计数允许边界随过去变化，不假设均匀边界生长。

此外，稳定第二峰提升的逐输出权重（同样不假设均匀边界生长）给

\[
 R_n=\frac{A_2(n)}{A_1(n)}\ge\frac N8-\frac32\ge\frac N{32}
 \qquad(n\ge4). \tag{9}
\]

这里 (9) 可由 \(R_n=\mathbb E_1S/\mathbb E_2\beta\)、
\(0\le\beta\le1\) 及逐源提升下界推出；也可直接用加权提升的输出纤维账本。
由 (7)、(9)、(N!\ge(N/e)^N)，设

\[
 \delta=e^{-50},\qquad r=\lceil N/2\rceil,
\]

当 \(\delta N<1\) 时坏事件为空；其余情形由计数界得

\[
\begin{aligned}
 \mathbb P_2\left\{\#\{j\le N-1:|\partial D_j|\le\delta N\}\ge r\right\}
 &\le32N^2e^{(18+\log2)N}\delta^{r-1}\\
 &\le32e^{50}N^2e^{-(7-\log2)N}.
\end{aligned} \tag{10}
\]

所以存在绝对 (n_0)，使全部 (n\ge n_0) 有

\[
 \boxed{
 \mathbb P_2\left\{\#\{j\le N-1:|\partial D_j|>e^{-50}N\}
 \ge\lfloor N/2\rfloor\right\}\ge1-e^{-5N}.
 } \tag{11}
\]

这比 r014 的 (n^{2-f(n)}=o(N)) 阈值更强：现在阈值是固定正比例 \(\delta N\)，
且至少固定正比例高度的边界大。常数 \(e^{-50}\) 很保守，只为使 (10) 的所有
常数方向显式；没有声称优化。

## 4. 对原题的范围

(11) 是关于双峰标号的严格几何粗糙性结论，但仍不是峰距结论。一个标号可以有
多数大的上层边界而两峰近，也可在远峰时如此；还没有从 (W(S,w)) 的坏空间矩形态
到 (11) 的反向蕴含。树只用于产生单峰标号子族和分母下界，并没有近似自然森林或
改变原始条件律。
