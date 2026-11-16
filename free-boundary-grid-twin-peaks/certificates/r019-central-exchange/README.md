# r019：中央切有效镜像交换的精确计数

对偶数 \(n=2m\)，令 \(L=\{x<m\}\)、\(R=\{x\ge m\}\)。每个同侧双峰
标号有唯一的“低峰” \(p\) 与“高峰/全图最大点” \(q\)。令 \(r\) 为 \(p\) 关于
中央竖轴的镜像，并令 \(t=(p\ r)\)。交换这两个点的标号后，输出在原图 \(G\)
有峰集 \(\{r,q\}\)，当且仅当原高到低 word 同时在 \(G\) 与
\(H=t^{-1}Gt\) 中有峰集 \(\{p,q\}\)。

`central_exchange_dp.cpp` 以该等价性进行子集 DP。它从 \(q\) 的单点前缀开始；
对待加入点 \(v\)，在 **两个** 图中同时检查

\[
(v=p)\Longleftrightarrow(N_G(v)\cap S=\varnothing),\qquad (v=p)\Longleftrightarrow(N_H(v)\cap S=\varnothing).
\]

因此 \(p\) 是唯一的第二出生点，其他点在两图中都必须已有较早邻点。所有 word
都恰被一次计数。对 \(n=2,3\)，程序把每个不同的 \((p,q,r)\)（\(p,q\) 非邻且
三点不同）与固定首点 \(q\) 的完整排列直接枚举交叉核验。

记 `same` 为两峰同侧的全部 \(K_2\) 标号数，`split` 为两侧各一峰的全部数，
而 \(E\) 为上述双图条件下的同侧源标号数。因为高峰 \(q\) 和镜像目标 \(r\) 都可
从成功输出恢复，固定“低峰”规则下这个有效域交换是单射，故预期并由程序断言
\(E\le\texttt{split}\)。输出也独立以普通单图 DP 重得 \(A_2(4)\)。

额外的必要条件说明有效域为何不是自动的：若 \(u\) 是 \(p\)-分支合并前的首个
后继，且 \(p\) 是 \(u\) 唯一较高邻点，又 \(u\ne r\) 且 \(r\not\sim u\)，则有效交换后 \(u\)
不得成为新峰，强制

\[
\lambda(r)>\lambda(u).
\]

这是直接比较交换后 \(u\) 的所有邻点得到的；其余邻点标号不变，而 \(p\) 的
新标号正是 \(\lambda(r)\)。它是一个严格的必要条件，不是对有效域比例的估计。

## 已运行的 \(4\times4\) 精确结果

程序通过 \(n=2\) 的 8 个和 \(n=3\) 的 336 个不同 \((p,q,r)\) 三元组的完整
排列交叉核验；所有双图 DP 计数相等。\(n=4\) 的输出为

\[
\begin{array}{c|r}
A_2&779299174080\\
\texttt{same}&191585447680\\
\texttt{split}&587713726400\\
E&88561436528.
\end{array}
\]

于是

\[
\theta:=E/\texttt{same}=\frac{5535089783}{11974090480}
\approx0.4622555502,
\qquad
\frac{\theta}{1+\theta}\approx0.3161250098.
\]

每个同侧有序源对都有正的有效词数，且程序逐对断言双图计数不超过原单图有序
计数，并总计断言 \(E\le\texttt{split}\)。独立调用上游 `pair_dp.cpp` 的
`pair_core_crosscheck.py` 重得同一 \(A_2\)、`same`、`split` 三数。实际运行、
资源和完整整数输出见 `n4-central-exchange*.json` 与
`n4-central-exchange-resource.txt`。

在本目录用 `clang++ -O3 -std=c++17 central_exchange_dp.cpp -o /Volumes/data/CodexTemp/T2/a022-central-exchange` 编译，并运行该临时可执行文件；stdout给完整精确JSON。`pair_core_crosscheck.py` 参数为由 `../../computation/pair_dp.cpp` 编译的临时可执行文件路径。所有验证进程已退出。
