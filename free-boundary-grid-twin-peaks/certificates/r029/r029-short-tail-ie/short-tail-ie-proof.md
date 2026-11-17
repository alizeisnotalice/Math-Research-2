# R029：穿孔纤维短间隙的独立峰容斥系数

## 记号和结论

令 $G=P_n\square P_n$，固定删除的最高点 $v$，并设

$$H=G-v,\qquad M=|V(H)|=n^2-1,\qquad N_v=N_G(v),\qquad Z=V(H)\setminus N_v.$$

记 $\mathcal D_v$ 为 $H$ 上根在 $Z$、且该根是唯一 $Z$-峰的词。对
$\zeta\in\mathcal D_v$，令 $p$ 是首点，$J$ 是 $p$ 的第一个 $H$-邻点的
位置，并令 $L=J-1$。对 $I\subseteq Z$ 独立，设

$$ C_H(I)=\#\{H\text{-标号}:I\text{ 中每点为峰}\},
 \quad U_I=V(H)\setminus\overline N_H(I),\quad P_I=Z\cap U_I. $$

因为 $I\subset Z$，$\overline N_H(I)=\overline N_G(I)$。下文所有边均为
有向选择的第一端点在左，写

$$e(A,U)=\sum_{p\in A}|N_H(p)\cap U|.$$

定义 $h_t=\#\{\zeta\in\mathcal D_v:L(\zeta)=t\}$。则正确的短系数为

$$
h_1=\sum_{I\subseteq Z\atop I\text{ independent}}
(-1)^{|I|}\frac{C_H(I)}{(M)_2}e(P_I,U_I),\tag{1}
$$

$$
h_2=\sum_{I\subseteq Z\atop I\text{ independent}}
(-1)^{|I|}\frac{C_H(I)}{(M)_3}
\sum_{q\in N_v\cap U_I}e(P_I\setminus\overline N_H(q),U_I).\tag{2}
$$

更一般地，令 $\mathcal E_t(U)$ 是长度 $t-1$ 的有序、互异顶点词
$\eta=(y_1,\ldots,y_{t-1})$，其顶点均在 $U$，并满足

$$y_1\in N_v,qquad y_j\in N_v\ \text{或}\ N_H(y_j)\cap\{y_1,\ldots,y_{j-1}\}\ne\varnothing\quad(j\ge2).$$

令 $\overline N_H(\eta)$ 是其顶点集的闭邻域；空词闭邻域为空。则

$$
h_t=\sum_{I\subseteq Z\atop I\text{ independent}}
(-1)^{|I|}\frac{C_H(I)}{(M)_{t+1}}
\sum_{\eta\in\mathcal E_t(U_I)}e(P_I\setminus\overline N_H(\eta),U_I).
\tag{3}
$$

式 (1) 与 (2) 正是题设候选。式 (3) 中“$\eta\subset U$”必须理解为顶点
集包含于 $U$ 且词内顶点互异；漏掉“有序/互异”会少计或重复计。

## 证明

先固定 $p\in Z$，并考虑 $L=t$ 的合法词前缀

$$ (p,y_1,\ldots,y_{t-1},x)=(p,\eta,x).\tag{4} $$

这里 $x\in N_H(p)$ 是第一个 $p$-邻点。由于除 $p$ 外的任何 $Z$ 点都不可
成为峰，$y_1$ 若在 $Z$ 将没有早邻点，故必在 $N_v$；对后续 $y_j$，若不在
$N_v$ 就必须邻接早先的某个 $y_i$。这正是 $\mathcal E_t$ 的生长条件。反过来
这个条件保证前缀中没有额外的 $Z$-峰。条件“$x$ 是首个 $p$-邻点”等价于

$$p\notin\overline N_H(\eta).\tag{5}$$

给定 $I$，若前缀所有点在 $U_I$，则删去该固定有序前缀、再插回其指定的最高
$t+1$ 个位置，严格给出

$$
\#\{I\text{ 均峰且前缀为 (4)}\}=\frac{C_H(I)}{(M)_{t+1}}.\tag{6}
$$

这里无遗漏的分母：因为前缀避开 $\overline N_H(I)$，插入任何这些最高标号都不
改变 $I$ 的峰约束；反向删除给 $(M)_{t+1}$ 个可能的有序最高标号插入。特别地，
只要一个这样的固定前缀存在，(6) 是整数。若对应的边/词和为零，公式中的该项
本来就是零，无须要求 $C_H(I)/(M)_{t+1}$ 本身可除：例如 $n=3$ 角孔、
$I=\{4\}$、$t=3$ 时 $C_H(I)=8064$ 而 $(M)_4=1680$，但 $U_I$ 不足以产生
有效前缀，故其边和为零。这个零项细节在模实现中无碍（先乘边和再除或在有限域
乘逆元都给零），但整数解释不得声称对所有 $I$ 单独整除。若某前缀点落在
$\overline N_H(I)$，该固定高标号会破坏对应的 $I$ 峰约束，所以不能省去
$U_I$ 限制。

固定 $\eta$ 后，(5) 和 (6) 所允许的有向末边 $(p,x)$ 数正是

$$e(P_I\setminus\overline N_H(\eta),U_I).$$

它自动排除 $p$ 与 $\eta$ 相同或相邻；从而 $x$ 也不可能等于已出现的
$y_j$（若 $x=y_j$，则 $p\in\overline N_H(\eta)$）。没有额外的“$x\notin\eta$”
因子或互邻修正。

最后对所有 $I\subseteq Z$ 独立作 $(-1)^{|I|}$ 容斥。首点 $p$ 已自动是峰；
上述前缀已经保证其余前缀 $Z$ 点不是峰。后缀中的全部额外 $Z$ 峰恰被这个容斥
消去，因为固定词的贡献是

$$
\sum_{I\subseteq M_H\cap Z\setminus\{p\}}(-1)^{|I|}
=1_{\{|M_H\cap Z|=1\}}.
$$

在 (3) 中，$U_I$ 的限制正好强制 $I$ 与所有前缀点及其邻域相容，所以可把和
写为所有 $I\subseteq Z$ 独立而不需另列“$I\cap\overline N(\text{prefix})=\varnothing$”。
这证明 (3)，并给出 (1)、(2)。

## 有限独立核验

`short_tail_ie_check.py` 由 $C_H(\varnothing)=M!$ 和精确递推

$$C_H(I)=\frac{\sum_{u\in I}C_H(I\setminus\{u\})}{|\overline N_H(I)|}$$

直接计算 (1)--(3)，完全不调用 R028 prefix DP 的状态递推。对每一个孔逐项与
R028 的 `tail-n3.json`、`tail-n4.json` 比较：$h_1,h_2,h_3$ 全部相等。

| $n$ | 孔数 | 结果 |
|---:|---:|:--|
| 3 | 9 | 全部通过 |
| 4 | 16 | 全部通过 |

例如中心孔：

| $n$ | $(h_1,h_2,h_3)$ |
|---:|:--|
| 3 | $(2688,1600,1024)$ |
| 4 | $(10102900080,5724050640,3607092720)$ |

这是恒等式及有限计算核验，不提供对 $\sum_{t\ge4}h_t$ 的控制，也不证明任何
全局尾收缩或 $B/A\ge1/5$。
