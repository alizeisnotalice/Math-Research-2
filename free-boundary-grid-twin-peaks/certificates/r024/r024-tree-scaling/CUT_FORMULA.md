# 树切边公式及其 hook-length 比式

令 \(B(T)\) 为恰有两个峰、且两个最高标签都是峰的标号数。若树边 \(e=uv\)
删去后给出 \(A\ni u\)、\(C\ni v\)，\(|A|=s\)、总顶点数为 \(N\)，则

\[
B(T)=\sum_{e=uv}\binom{N-2}{s-1}
\left(A_1(A)A_1(C)-a_u(A)a_v(C)\right). \tag{C1}
\]

其中 \(a_r(H)\) 是以 \(r\) 为唯一峰的标号数，\(A_1(H)=\sum_r a_r(H)\)。
在树上，两个峰的升邻居定向恰有一个入度二的汇合点；切掉其两条入边所得两个
单峰分量，按交错恢复，给出 (C1)。若两单峰根正为边端 \(u,v\)，它们相邻，不能
同时是两最高峰，故要减去第二项。该式随后由 `verify_cut_formula.py` 以独立hook
计数与积分消息DP逐值比对。

若全树单峰根取在 \(C\) 中的 \(r\)，不能把两个 rooted-forest 的 hook
分母直接相乘：\(u\) 的祖先的子树大小会改变。把 \(T\) 以 \(r\) 为根，令
\(v\) 为 \(u\) 的父点，且令 \(t_w\) 是全树中路径 \(r\rightsquigarrow v\) 上
各非根点 \(w\) 的 rooted subtree 大小。正确的 hook-length 恒等式是

\[
a_r(T)=\binom{N-1}{s}a_u(A)a_r(C)G_e(r),\qquad
G_e(r)=\prod_{w\in(r\rightsquigarrow v],\,w\ne r}\frac{t_w-s}{t_w}. \tag{C2}
\]

例如三点路径以端点为根、切去另一端叶时，左边为 \(1\)，而未含 \(G\) 的右边为
\(2\)；此时 \(G=1/2\)。它是本公式的回归例。

因此 (C1) 的非负主项除以 \(a_r(T)\) 是

\[
\frac{s}{N-1}\frac{A_1(A)}{a_u(A)}\frac{A_1(C)}{a_r(C)}\frac1{G_e(r)}. \tag{C3}
\]

这是后续控制满二叉桥树 \(B(T)\) 的解析入口。本文件只证明恒等式及有限交叉核验；
所有高度的因子统一上界及由此导出的渐近结论已在独立证明
[`../r024-global-normalizer/double-binary-b-upper.md`](../r024-global-normalizer/double-binary-b-upper.md)
给出，不能从本文件本身自动推出。
