# R036：自由边界方格的 dyadic 壳层回插下界

## 结论及范围

令 $G=P_n\square P_n$，$N=n^2$。记 $A_1(n)$ 为一峰降序词数，$B(n)$ 为
恰有两个峰、且最高两个标号恰为这两个峰的降序词数（即本项目的 top-two 双峰
切片）。则对每个
$n\ge2$，

$$
\boxed{\quad
 \frac{B(n)}{A_1(n)}\ \ge\
 \frac{\bigl[\lfloor\log _2\lfloor n/2\rfloor\rfloor-1\bigr]_+}{32}
 =\frac{\bigl[\lfloor\log _2n\rfloor-2\bigr]_+}{32}.
\quad} \tag{D1}
$$

这里 $[x]_+=\max\{x,0\}$。因此，严格地，

$$
\frac{B(n)}{A_1(n)}\ge\frac{\log _2 n}{32}-O(1),\qquad n\to\infty. \tag{D2}
$$

下界从 $n\ge8$ 起非平凡。它只比较一峰切片和 top-two 双峰切片；本附件
**不**断言原问题的双峰条件律已经集中，也不把此式替换为 $A_2/B$ 的界。

## 回插的调和计数引理

以下秩一律以 $1$ 起算。令
$\sigma=(v_1,\ldots,v_N)$ 是一峰降序词，根为 $r=v_1$。对于
$p\notin\{r\}\cup N(r)$，写

$$
 l_p=\min_{u\sim p}\operatorname {rank}_\sigma(u),
$$

并令 $b_p$ 是原词中最早的顶点 $v$，满足它晚于 $p$ 出现而其唯一早邻点为
$p$；若没有这样的 $v$，令 $b_p=N+1$。将 $p$ 稳定地移至第二位所得的词
记为 $\tau(\sigma,p)$。删去该第二位后，恰能在原词秩

$$
 l_p+1,\ l_p+2,\ \ldots,\ b_p-1
$$

回插并仍得到一峰词。因此该双峰词在这类提升下的逆像大小是

$$
 R_p=b_p-l_p-1. \tag{D3}
$$

的确，回插若早于 $l_p+1$，$p$ 没有早邻点；若不早于 $b_p$，定义 $b_p$ 的
顶点失去其唯一早邻点。两端之间，被越过的每个顶点原本不是仅以 $p$ 为早邻点，
故删去或后移 $p$ 后仍有另一早邻点。这也证明所有中间位置均合法。

对每个提升对赋予权 $1/R_p$。固定一个位于提升像中的双峰词，其全部逆像贡献
正好为 $1$，而像外双峰词贡献为零。故

$$
 \frac {B(n)}{A_1(n)}\ge
 \mathbb E_{\sigma\in\mathcal A_1}
 \left[\sum_{p\notin\{r\}\cup N(r)}\frac1{R_p}\right]. \tag{D4}
$$

这是正权计数，既不需要满像，也不需要对双峰词作均匀性假设。

## 壳层首次到达的唯一父点

固定这样的 $\sigma$。写

$$
 d(v)=\|v-r\|_\infty,\qquad D=\max_{v\in G}d(v),
 \qquad T_k=\min\{t:d(v_t)=k\}\quad(0\le k\le D).
$$

一峰词的每个前缀诱导连通：归纳地，词中每一个非根顶点都有一个更早邻点，
故加入该点保持前缀连通。每个 $T_k$ 存在：取一条从 $r$ 到距离 $D$ 的顶点的
格路，其 Chebyshev 距离逐边至多变化 $1$，必经过所有中间壳层。前缀连通性又给出

$$
 t<T_k\ \Longrightarrow\ d(v_t)\le k-1. \tag{D5}
$$

否则首个出现于壳层 $\ge k$ 的顶点以前的前缀，连同该点仍连通；从该点到 $r$ 的
此前缀路径在首次越过从距离 $<k$ 到距离 $\ge k$ 时给出一个更早的壳 $k$ 顶点，
矛盾。

令 $x_k=v_{T_k}$。若 $x_k$ 的两个坐标都达到 Chebyshev 距离 $k$，则它没有
任何邻点位于严格内球 $\{d\le k-1\}$：每次只改变一个坐标，另一个坐标的距离
仍为 $k$。这与 (D5) 及 $x_k$ 必有早邻点矛盾。故 $x_k$ 恰有一个坐标达到
距离 $k$。于是它恰有一个邻点 $p_k$ 位于壳 $k-1$，其余邻点均位于壳 $k$ 或
更外。由 (D5)，$p_k$ 是 $x_k$ 的唯一早邻点。特别地，

$$
 d(p_k)=k-1,
 \qquad k\ge3\ \Longrightarrow\ p_k\notin\{r\}\cup N(r). \tag{D6}
$$

不同的 $k$ 给出不同的 $p_k$，因为它们在不同的 Chebyshev 壳上。

现在 $p_k$ 的每个邻点均有距离至少 $k-2$，故由 (D5)

$$
 l_{p_k}\ge T_{k-2}.
$$

而 $x_k$ 本身证明 $b_{p_k}\le T_k$。代入 (D3) 得到关键的逐壳不等式

$$
 R_{p_k}\le T_k-T_{k-2}-1 \qquad (3\le k\le D). \tag{D7}
$$

右侧自动为正，因为 $R_{p_k}\ge1$。

## Dyadic 求和

取 $m=2^j\ge2$ 且 $2m\le D$。把 (D7) 对
$k=m+1,\ldots,2m$ 相加，得到

$$
\begin{aligned}
 \sum_{k=m+1}^{2m}R_{p_k}
 &\le T_{2m-1}+T_{2m}-T_{m-1}-T_m-m\\
 &\le 2T_{2m}. \tag{D8}
\end{aligned}
$$

在自由边界方格中，半径 $s$ 的 $\ell_\infty$ 球不超过无限方格中相应球的大小，
即

$$
 \bigl|\{v:d(v)\le s\}\bigr|\le(2s+1)^2.
$$

首次进入壳 $k$ 以前只可能出现此球内的点，故

$$
 T_k\le(2k-1)^2+1.
$$

于是

$$
 \sum_{k=m+1}^{2m}R_{p_k}
 \le2\bigl[(4m-1)^2+1\bigr]
 =32m^2-16m+4\le32m^2. \tag{D9}
$$

Cauchy--Schwarz 给出每一个这类块的确定性贡献

$$
 \sum_{k=m+1}^{2m}\frac1{R_{p_k}}
 \ge\frac{m^2}{\sum_{k=m+1}^{2m}R_{p_k}}
 \ge\frac1{32}. \tag{D10}
$$

对于 $m=2,4,8,\ldots$，这些 $k$ 区间两两不交，且相应的 $p_k$ 也两两不同。
可用的块数为

$$
 \bigl[\lfloor\log_2D\rfloor-1\bigr]_+.
$$

因而每个一峰词逐词满足

$$
 \sum_{p\notin\{r\}\cup N(r)}\frac1{R_p}
 \ge\frac{\bigl[\lfloor\log_2D\rfloor-1\bigr]_+}{32}. \tag{D11}
$$

最后，任意 $r\in P_n\square P_n$ 到四个角之一的 Chebyshev 距离至少为
$\lceil(n-1)/2\rceil=\lfloor n/2\rfloor$，所以
$D\ge\lfloor n/2\rfloor$。恒等式
$\lfloor\log_2\lfloor n/2\rfloor\rfloor-1=\lfloor\log_2n\rfloor-2$
（在正部内理解）给出 (D1)。将 (D11) 代入 (D4) 完成证明。

## 一个严格的早首次合并推论

这段只记录 (D1) 与已有的首次合并编码之间的直接接口，仍不涉及峰距的收缩。
令 $B_q$ 是 $B$ 中首次合并前两支总大小为 $q$ 的词数，且令
$\mathbb P_B$ 为 $B$ 上的均匀律。对任意有限连通图，父路径编码给出

$$
 B_q\le3^{q-2}A_1\qquad(q\ge2). \tag{D12}
$$

故若 $L_n$ 表示 (D1) 的右端，则在 $n\ge8$（此时 $L_n>0$）及任意
$2\le k\le N-1$ 下有

$$
 \mathbb P_B(q\le k)
 \le\frac{\sum_{q=2}^k3^{q-2}A_1}{B}
 \le\frac{3^{k-1}-1}{2L_n}. \tag{D13}
$$

特别地，对 $0<\delta<1$、$n\ge64$，令

$$
 k=\left\lfloor(1-\delta)\log_3(\log_2n)\right\rfloor .
$$

若 $k<2$，事件为空；否则 $L_n\ge(\log_2n)/64$，从而

$$
 \mathbb P_B(q\le k)\le\frac{32}{3}(\log_2n)^{-\delta}. \tag{D14}
$$

这里 $L_n\ge(\log_2n)/64$ 使用
$\lfloor\log_2n\rfloor-2\ge(\log_2n)/2$（$n\ge64$）。在完整双峰律下，
首次合并纤维的逐点关系 $m\le q/2$ 还给出

$$
 \mathbb P_2(q\le k)\le\frac{k}{2}\,\mathbb P_B(q\le k), \tag{D15}
$$

因为 $A_2\ge B$。所以同一 $k$ 下完整双峰律的早合并质量也趋于零；这仅排除
首次合并秩很小的部分，**并不**闭合原问题的峰距估计。

## 有限核验（不承担证明）

`verify_dyadic_shell_samples.py` 对固定种子生成的
$4\times4,5\times5,6\times6,8\times8,16\times16,32\times32$ 一峰生长词逐条
断言 (D5)--(D9) 及回插秩不等式；输出
`dyadic-shell-sample-check.json`。这只是对索引和边界情况的轻量核验，严格结论由
上述逐词证明承担。
