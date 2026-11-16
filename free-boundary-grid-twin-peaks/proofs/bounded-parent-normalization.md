# 有界父容量的提前桥构造

状态：严格单分割计数工具。它给双峰标号数的一个下界，不能自动控制所有分割的重数或原题概率。

令有限连通图 \(G\) 的顶点分割为非空连通且互不邻接的 \(A,B\)，写
\(a=|A|,b=|B|,M=a+b\)，其余点集为 \(C\)，大小 \(g\)。固定 \(C\) 中的
\(k\ge1\) 个点，并为每点指定一个在 \(A\cup B\) 中的邻父。假定每个父最多
指定 \(\Delta\) 个孩子。图最大度至多 \(\Delta\) 时，任意这种指定自动满足容量。

令 \(q=\lceil k/\Delta\rceil\)，\(r=k-\Delta(q-1)\in\{1,\ldots,\Delta\}\)，并假定
\(q\le M-1\)。取两条内部一峰词，先把两根放在全局词的头两个位置。这样的 A/B 交错数是
\[
 2\binom{M-2}{a-1}.
\tag{1}
\]

按父在这个 A/B 词中的秩非降地排列指定孩子，记父秩为 \(r_i\)。容量给出
\[
 k-i+1\le\Delta(M-r_i+1),\qquad
 r_i\le t_i:=M-\left\lfloor\frac{k-i}{\Delta}\right\rfloor.
\tag{2}
\]
前一式因为末尾 \(k-i+1\) 个孩子的父都在秩 \(r_i,\ldots,M\) 中，而每点容量至多 \(\Delta\)。

强制 A/B 词的前 \(M-q\) 点先出现。把最后 q 个连续 A/B 点记为 gates
\(G_1,\ldots,G_q\)。令前 r 个孩子在 \(G_1\) 后，随后每 \(\Delta\) 个孩子分别在
\(G_2,\ldots,G_q\) 后；每一孩子的父秩由 (2) 不超过其 gate 秩，所以实际父已出现。
gate 与孩子形成一棵根为 \(G_1\) 的有根树：第一个 gate 挂 r 个叶及 \(G_2\)，以后每个
gate 挂 \(\Delta\) 个叶及下一个 gate。hook-length 公式给其线性扩张数
\[
 E_{k,\Delta}=\frac{(k+q-1)!}{(\Delta+1)^{q-1}(q-1)!}.
\tag{3}
\]
这里 \(G_1\) 是该子词第一点；\(q\le M-1\) 确保它不会位于两根之前。

每个选中孩子在其父之后出现，故不会形成新峰；它不必同时邻接两支。令 \(C'\) 为选中
孩子集，\(Z=A\cup B\cup C'\)，并令 \(h_G(Z)\) 是从 Z 出发、每步附着到当前已揭示
集合的 R=\(C\setminus C'\) 完成词数。即使 Z 有两个分量，连通 G 也保证只要未完成就存在
它们并集的外邻点，故 \(h_G(Z)>0\)。允许所有这 \(h_G(Z)\) 个后缀。固定分割、孩子、父
和输出词后，可恢复 A/B 子词、root-first 交错、gate 树线性扩张和完整 R 后缀，故无重数。
得到
\[
 D_{\rm new}=a(A)a(B)\,2\binom{M-2}{a-1}E_{k,\Delta}h_G(Z)\le A_2(G).
\tag{4}
\]

令 \(F_{\rm num}=\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\})\) 为原首次合并账本中的
完成数和。任一被 Fnum 计的 C 词删除固定 C' 后，R 相对词在 Z 已全揭示时仍合法：原较早
邻点若删去便已经在 Z，若保留则仍较早。给定 R 词，插回 k 个有区别的 C' 点的方式至多
\(\binom gk k!=(g)_k\)。故
\[
 F_{\rm num}\le(g)_k h_G(Z).
\]
由 \(|\mathcal F|=a(A)a(B)\binom MaF_{\rm num}\)，有
\[
 \frac{|\mathcal F|}{D_{\rm new}}
 \le\frac{M(M-1)}{2ab}\frac{(g)_k}{E_{k,\Delta}}.
\tag{5}
\]

若 \(X\sim\operatorname{Bin}(k+q-1,1/(\Delta+1))\)，则 \(q-1\) 是一个模式，
故 \(\mathbb P(X=q-1)\ge1/(k+q)\)。由 (3)，这实际上给出稍强于所提形式的
\[
 \frac{E_{k,\Delta}}{k!}
 \ge\frac{(1+1/\Delta)^k}{k+q}
 \ge\frac{(1+1/\Delta)^{k-\Delta}}{k+q}.
\tag{6}
\]
代入 (5) 并用 \((g)_k/k!=\binom gk\)，得到更强的界
\[
 \frac{|\mathcal F|}{D_{\rm new}}
 \le\frac{M(M-1)}{2ab}\binom gk(k+q)(1+1/\Delta)^{-k}.
\tag{7}
\]
使用题述较弱的 (6) 的第二个不等式，可把最后指数改写为 \(\Delta-k\)。这避免了旧式
\(g^{g-k}\) 损失，但仍不控制分割总数。

## 方格上的可用常数推论

在自由方格中取 \(\Delta=4\)。对固定 \(\eta>0\)，若一个固定分割及其指定父满足
\[
 k\ge\frac{63}{64}g,\qquad q\le M-1,\qquad ab\ge\eta M^2,
\tag{8}
\]
则由 (7) 和 \(D_{\rm new}\le A_2\)，有
\[
 \frac{|\mathcal F|}{A_2}
 \le\frac{|\mathcal F|}{D_{\rm new}}
 \le\frac{g}{\eta}e^{-g/10}.
\tag{9}
\]
确实，令 \(s=g-k\le g/64\)。若 \(s=0\)，二项式项等于 \(1\)，结论直接成立；否则
\[
 \binom gk=\binom gs\le\left(\frac{eg}{s}\right)^s
 \le\exp\left\{\frac{1+\log64}{64}g\right\}.
\]
又 \(\log(5/4)\ge1/5\)、\(\log2<7/10\)，所以
\[
 \frac{1+\log64}{64}-\frac{63}{64}\log(5/4)
 \le\frac{26-63}{320}=-\frac{37}{320}<-\frac1{10}.
\]
最后 \(q=\lceil k/4\rceil\le k\)，故 \(k+q\le2g\)，而
\(M(M-1)/(2ab)\le1/(2\eta)\)。代入 (7) 即得 (9)。

这是**单一固定分割及固定指定父**的指数压低。要用于一类 boundary-dense 形态，仍须另行
控制可选分割、选点和父指定的数量或做一个无重数的全局聚合；(9) 本身不能对所有这类形态求和。

当 \(\Delta=1\) 时，\(q=k,r=1\)，(3) 化为
\[
 E_{k,1}=\frac{(2k-1)!}{2^{k-1}(k-1)!},
\]
正是每个 gate 附一叶的 matching gate-tree 计数。

`audit_bounded_parent.py` 穷举 \(1\le k\le6,1\le\Delta\le4\) 的 gate-tree 全部线性扩张，
核对 (3)、(6)；它还在完整 3×3 方格的 T 分割上允许所有 \(h_G(Z)\) 后缀，实际生成
所有构造词，独立核对每词恰有两次出生、无重复和限制词界。
同一审计另取底行左端这一只邻 A、不邻 B 的选中点，仍实际生成 384 个互异双峰词；这核对
了本命题不需要选中点是共同桥。
