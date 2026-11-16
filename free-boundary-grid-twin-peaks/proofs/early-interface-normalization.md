# 提前界面根的晚合并归一化引理

设有限连通简单图 \(G\) 的顶点不交分割为
\[
V=A\mathbin{\dot\cup}B\mathbin{\dot\cup}C,\qquad
|A|=a\ge1,\ |B|=b\ge1,\ |C|=g\ge1,
\]
其中 \(G[A]\)、\(G[B]\) 连通且 \(A,B\) 无边。写
\[
\mathcal W=\{w\in C:N_G(w)\cap A\ne\varnothing,\ N_G(w)\cap B\ne\varnothing\}.
\]
令 \(h_A(C)\) 是已经揭示全部 \(A\) 后，无新出生地揭示全部 \(C\) 的排列数；先假定
\(h_A(C)>0\)。令 \(\mathcal F\) 是先完整揭示 \(A,B\)，再由某 \(w\in\mathcal W\)
首次合并的双峰词族。首次合并账本给出
\[
|\mathcal F|=a(A)a(B)\binom{a+b}{a}
\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\}).
\tag{1}
\]

## 提前出现 \(B\) 根

有严格界
\[
\frac{|\mathcal F|}{A_2(G)}
\le\frac{\binom{a+b}{a}\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\})}
{h_A(C)\left[\binom{a+b+g}{b}-\binom{b+g-1}{b}\right]}.
\tag{2}
\]

固定 \(A,B\) 的两条一峰词与一条由 \(h_A(C)\) 计数的 \(C\)-词。令 \(H\) 是完整
\(A\)-词后接该 \(C\)-词。将 \(H\) 与 \(B\)-词交错，保留 \(B\) 根在 \(H\) 首个
\(C\) 点之前的交错。全部交错数为 \(\binom{a+b+g}{b}\)。坏交错中，第一个 \(B\)
点在首个 \(C\) 点之后；开头遂已放置全部 \(a\) 个 \(A\) 点和一个 \(C\) 点。删去该
强制前缀，余下 \(b\) 个 \(B\) 点与 \(g-1\) 个 \(H\) 点任意交错，故坏交错恰为
\[
\binom{b+g-1}{b}.
\tag{3}
\]
每个保留词恰有两个出生：两支根都先于所有 \(C\) 点，且两支无边；支内非根各有
更高的同支早邻点；\(C\) 点沿其 \(h_A(C)\)-词已有 \(A\) 或早到 \(C\) 邻点。
限制到 \(A,B,C\) 的三个相对词和二元交错词唯一恢复数据。因此
\[
D_B=a(A)a(B)h_A(C)
\left[\binom{a+b+g}{b}-\binom{b+g-1}{b}\right]\le A_2(G),
\tag{4}
\]
以 (1) 除之即得 (2)。

## 双向选择

交换 \(A,B\) 得
\[
D_A=a(A)a(B)h_B(C)
\left[\binom{a+b+g}{a}-\binom{a+g-1}{a}\right].
\tag{5}
\]
所以
\[
\frac{|\mathcal F|}{A_2(G)}
\le\frac{\binom{a+b}{a}\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\})}
{\max\{R_A,R_B\}},
\tag{6}
\]
\[
\begin{aligned}
R_A&=h_B(C)\left[\binom{a+b+g}{a}-\binom{a+g-1}{a}\right],\\
R_B&=h_A(C)\left[\binom{a+b+g}{b}-\binom{b+g-1}{b}\right].
\end{aligned}
\tag{7}
\]
真正待估的是 \(\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\})\) 相对于
\(\max\{R_A,R_B\}\) 的大小。该表述包含最优的两侧交换，单向失败不能否定它。

## 完整 \(P_n\square P_n\) 上的单向阶乘例

在完整自由边界方格 \(P_n\square P_n\)（\(n\ge3\)）中，取
\[
A=\{(0,0)\},\quad
C=\{(x,0):1\le x\le n-2\}\cup\{(0,1)\},\quad
B=V\setminus(A\cup C).
\tag{8}
\]
这里 \(G[A]\)、\(G[B]\)、\(G-A\)、\(G-B\) 都连通，且 \(A,B\) 无边。令
\(p_i=(i,0)\) 及 \(q=(0,1)\)。从 \(A\) 生长时，路径
\(p_1,\ldots,p_{n-2}\) 的次序被强制，\(q\) 可插入任一位置，故
\[
g=n-1,\qquad h_A(C)=g.
\tag{9}
\]
两座桥仅为 \(p_1,q\)。完整 \(B\) 已揭示后，\(C\) 的每一点都邻接 \(B\)，所以每座
桥的后缀完成数为 \((g-1)!\)，从而
\[
\Gamma(A,B,C):=\frac{\sum_{w\in\mathcal W}h_G(A\cup B\cup\{w\})}{h_A(C)}
=\frac{2(g-1)!}{g}.
\tag{10}
\]
这是双方补集连通时的阶乘级单向增益。然而交换 \(A,B\) 后，每个 \(C\) 点都邻接
新的 \(A=B\)，故
\[
h_B(C)=g!,\qquad \Gamma(B,A,C)=\frac2g.
\tag{11}
\]
所以该例不反驳 (6)，而是要求后续论证选择较优方向。

## 可用条件及范围

若每个 \(C\) 点都邻接 \(A\)，则
\[
h_A(C)=g!,\qquad
h_G(A\cup B\cup\{w\})=(g-1)!\ (w\in\mathcal W),\qquad
\Gamma(A,B,C)=|\mathcal W|/g.
\tag{12}
\]
方格的一列隔离带与同心宽度一环带满足此条件。一般首次合并形态是否能使 (6) 的
双向完成数比衰减，仍是本路线的第一未证量。

## 核验

computation/audit_early_interface.py 枚举小的 \((a,b,g)\) 交错词，逐项验证 (3) 及好词数，
并在七顶点图逐词检查好词恰有两个出生。它还对 (8) 的完整方格例在
\(n=3,\ldots,8\) 逐个枚举 \(C\)-排列，核对 \(h_A(C)=g\)、\(h_B(C)=g!\) 和两座桥的
\((g-1)!\) 后缀数。对同心参数
\((a,b,g,c_0)=(16,1,8,4)\)，晚合并因子为
\[
4\binom{17}{16}7!=342720,
\]
而 (4) 的比较族因子为
\[
8!\left[\binom{25}{1}-\binom81\right]=685440.
\]
所以加强比为 \(1/2\)，旧的“\(B\) 根全局第一”比较只给 \(17/2\)。

在双向式中允许某方向完成数为零，只要求两方向至少一个为正。因 \(A\) 连通，\(h_A(C)>0\) 当且仅当诱导图 \(G[V\setminus B]\) 连通；这由存在从 \(A\) 出发的生成树顺序与反向连通性直接得到。若两方向均为零，本引理不适用。不得把所有单分割上界相加后忽略分割数量或比较族间重数。
