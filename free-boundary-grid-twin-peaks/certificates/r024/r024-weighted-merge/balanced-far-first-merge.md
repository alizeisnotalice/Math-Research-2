# 平衡首次合并纤维仍可有宏观根距

令 (n\ge3)，(k=\lfloor(n-1)/2\rfloor)。在 (P_n\square P_n) 的顶行取

\[
A=\{(x,0):0\le x<k\},\quad w=(k,0),\quad
B=\{(x,0):k+1\le x\le2k\}.
\]

两个支都是诱导路径 (P_k)，彼此没有边，而 (w) 分别邻接两个端点。因此它们是合法的首次合并三元组；合并后集合是一段连通顶行，故至少存在一个逐边界补全词，(W(A,B,w)>0)。

路径 (P_k) 的唯一峰根律是位置参数 (\operatorname{Bin}(k-1,1/2))。故两个根独立时

\[
\mathbb E\bigl[\|X-Y\|_2^2\mid A,B,w\bigr]
=(k+1)^2+\frac{k-1}{2}.
\]

另一方面，固定首次合并数据的 top-two / all-two interleaving balance 为

\[
\theta(A,B)=\frac{2k^2}{(2k)(2k-1)}=\frac{k}{2k-1},
\qquad
m(A,B)=\theta^{-1}=2-\frac1k<2.
\]

因此这是一族同时满足 (m<2)、\(\theta>1/2\)，但条件归一化二阶根距趋于 (1/4) 的真实正权重方格状态。

这不说明该族的总 (W/A_2\) 质量不消失，故不构成原题反例。它严格排除的只是一个潜在捷径：仅证明首次合并不平衡权 (m) 一致有界或可积，不能单独推出 (Q_n\to0)。还必须控制平衡纤维自身的根几何，或者证明这类宏观几何的总 (W)-质量趋零。

`verify_balanced_far_arms.py` 以 Fraction 实际检查 (3\le n\le100) 的图形条件及上述公式；输出 `balanced_far_arms.json`。
