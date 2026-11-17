# 两行释放的精确纤维公式

以下固定 \(n\ge3\)；\(n=2\) 时 \(M=0\) 是不含 birth 的单独退化情形，不在本记录的陈述内。令 \(T\) 为第一行，\(B\) 为第二行，\(R=T\cup B\)，余下方格记作 \(H\)，\(M=|H|=n(n-2)\)。固定一个从 \(R\) 出发的完成词 \(\sigma=(\sigma_1,\ldots,\sigma_M)\)。在 \(H\) 内无较早邻点的 birth 位于第三行，位置为
\[
1=r_1<r_2<\cdots<r_p.
\]
它们相互不邻接，故 \(1\le p\le a=\lceil n/2\rceil\)。

删去第二行后等于此 \(\sigma\) 的全部 \(T\)-完成词数恰为
\[
 C(\sigma)=\frac{(M+n)!}{(M+p)!}\prod_{i=1}^p(r_i+i-1).\tag{1}
\]
事实上，第 \(i\) 个 birth 的正上方第二行点 gate 必须在 \(\sigma_{r_i}\) 前，且这是唯一约束；其余第二行点都邻接 \(T\)。先将 gates 逐一插入链 \(\sigma\)。此前 gates 都在 \(\sigma_{r_i}\) 前，因此第 \(i\) gate 有 \(r_i+i-1\) 个可插入缝隙。其余 \(n-p\) 点自由交错，数为 \((M+n)!/(M+p)!\)。删去第二行反向恢复 \(\sigma\)、birth 集、gate 顺序及自由交错，故无重数。

因此若 \(H_T=h(T),H_R=h(R)\)，则
\[
H_T=\sum_{\sigma\in\mathcal C_R}C(\sigma).\tag{2}
\]

设 \(E_p=(2p-1)!!\)。由 \(r_i\ge i\)，
\[
C(\sigma)\ge E_p\frac{(M+n)!}{(M+p)!}.
\]
右侧从 \(p\) 到 \(p+1\) 的比值为 \((2p+1)/(M+p+1)<1\)，故
\[
\frac{H_T}{H_R}\ge E_a\frac{(M+n)!}{(M+a)!}.\tag{3}
\]
另一方面，
\[
\frac{C(\sigma)}{(M+n)_n}=\prod_{i=1}^p\frac{r_i+i-1}{M+i}
\le\frac1{M+1},
\]
因 \(r_1=1\)，其余因子至多一。因此
\[
\frac{H_T}{H_R}\le\frac{(M+n)!}{(M+1)!}.\tag{4}
\]
对 \(n=3\)，四个 \(p=1\) 词各有纤维 30，两个 \(p=2\) 词各有纤维 18，故 \(H_T=4\cdot30+2\cdot18=156\)。

这是完成数的精确分解，不比较 \(H_T\) 与全体一峰词，因而不解决原双峰问题。
