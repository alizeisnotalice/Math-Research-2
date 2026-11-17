# R037：稳定提升像的精确 singleton 刻画

## 结论

设 $G$ 为有限连通简单图。令 $\mathcal A_1$ 为一峰降序词集合，令 $B$ 为恰有
两个峰且这两个峰占据最高两个标号的词集合。对

$$
 \mathcal D=\{(\sigma,p):\sigma\in\mathcal A_1,
  \ r=\sigma_1,\ p\notin\{r\}\cup N(r)\},
$$

令 $\Phi(\sigma,p)$ 是把 $p$ 从 $\sigma$ 中稳定删出、再放到根 $r$ 后第二位所得
的词。令 $B_{\mathrm{low},1}$ 为 $B$ 中满足下述条件的子集：第二位的低根 $p$ 在
**第一次合并之前的低根分支恰为 singleton $\{p\}$**。则

$$
 \boxed{\quad\operatorname{Im}\Phi=B_{\mathrm{low},1}.\quad} \tag{1}
$$

若 $l_p,b_p$ 及 $R_p=b_p-l_p-1$ 按 R036 的回插定义，则每个
$\tau\in B_{\mathrm{low},1}$ 的实际逆像数恰为 $R_p$，其中 $p=\tau_2$。故严格地

$$
 \boxed{\quad
 \mathbb E_{\sigma\in\mathrm{Unif}(\mathcal A_1)}
 \sum_{p\notin\{r\}\cup N(r)}\frac1{R_p}
 =\frac{|B_{\mathrm{low},1}|}{A_1(G)}.
 \quad} \tag{2}
$$

因此 R036 的 dyadic 论证实际上给出方格的更精确下界

$$
 \frac{|B_{\mathrm{low},1}(P_n\square P_n)|}{A_1(P_n\square P_n)}
 \ge\frac{[\lfloor\log_2 n\rfloor-2]_+}{32}. \tag{3}
$$

## 像包含于低根 singleton 切片

取 $(\sigma,p)\in\mathcal D$。在 $\sigma$ 中，令 $w$ 为 $p$ 的最早出现邻点。
因为 $p$ 不是根而 $\sigma$ 一峰，$w$ 在 $p$ 之前出现。又 $p\not\sim r$，所以
$w$ 不是根；一峰性质给出 $w$ 在 $\sigma$ 中有一个早邻点 $h$。该 $h$ 不可能是
$p$，故属于根 $r$ 的早连通分支。

稳定提升保持除 $p$ 之外的一切相对次序。因此在
$\tau=\Phi(\sigma,p)$ 中，$w$ 仍是 $p$ 的首个邻点；在 $w$ 前，$p$ 没有邻点，
而 $h$ 已在根分支中。于是 $w$ 正是第一次把 singleton 低根分支 $\{p\}$ 与高根
分支相连的顶点。故 $\tau\in B_{\mathrm{low},1}$。

这里也验证了前向词确实恰有两峰：提升 $p$ 只新造峰 $p$；其余顶点原有的早邻点
不因 $p$ 被提前而消失。

## 低根 singleton 切片包含于像

反之，取 $\tau\in B_{\mathrm{low},1}$，写根为 $r=\tau_1$、低根为
$p=\tau_2$，并令 $w$ 为 $p$ 的首个邻点。由于低根分支在第一次合并前为 singleton，
$w$ 有一个早邻点 $h\ne p$。又在 $w$ 之前没有 $p$ 的邻点；所以删去 $p$ 后，
$w$ 之前的所有顶点仍逐步连到根 $r$，且 $w$ 仍借 $h$ 连到该根分支。

从 $\tau$ 稳定删去 $p$，并将它**紧接在 $w$ 后**插回，得到 $\sigma$。在插入
之前的前缀刚已说明是一峰连通前缀；插入 $p$ 时它邻接较早的 $w$；其后每个顶点
若原来唯一的早邻点为 $p$，则 $p$ 现在仍较早，其他顶点的早邻点不受影响。故
$\sigma$ 是一峰词。又 $p\not\sim r$（否则 $p$ 不会是第二峰），所以
$(\sigma,p)\in\mathcal D$，并且 $\Phi(\sigma,p)=\tau$。

这也说明至少有一个合法回插槽：在该 $\sigma$ 中
$l_p<\operatorname{rank}_\sigma(p)<b_p$，从而 $R_p\ge1$。R036 的逐槽论证给出
同一 $\tau$ 的全部逆像槽恰为 $R_p$ 个。对每个输出纤维以 $1/R_p$ 加权即贡献一，
从而得到 (2)。

## 一个必要的措辞修正

“直接删去 $p$ 后整词必是一峰词”并不正确；正确的反向操作是上节的**在 $w$ 后
回插 $p$**。例如在 $P_3\square P_3$（按行编号 $0,\ldots,8$）中，

$$
 \tau=(0,2,1,3,5,4,6,7,8)
$$

属于 $B_{\mathrm{low},1}$，第一次合并前的两个分支为 $\{0\},\{2\}$；但删去
$p=2$ 后，$5$ 成为额外 birth，所得词的 births 为 $(0,5)$。这不反驳 (1)，因为
将 $2$ 插到首个邻点 $w=1$ 之后即可恢复一峰词。

## 与已有 singleton / $\beta$ 账本的关系

这不是重证 r026 的完整双峰律恒等式。已归档的
`certificates/r026/r026-singleton-slice/SINGLETON_SLICE_IDENTITY.md`
（SHA-256 `8bebf55c1708392f488f93f1b5ab4c19c0c459f39e67000ce935660b4455b99e`）
已证明恰一 singleton 纤维的交错账本。结合两侧皆 singleton 的 $q=2$ 纤维，若
$B_{\mathrm{s}}$ 表示 $B$ 中至少一侧为 singleton 的词数，则

$$
 \boxed{\quad B_{\mathrm{s}}=2|B_{\mathrm{low},1}|-B_{q=2}.\quad} \tag{4}
$$

理由是：恰一 singleton 的 top-two 两种首位交错中恰一词的 singleton 为低根；
而 $q=2$ 时两根均为 singleton，须在两份有向计数中扣回一次。

另据既有 `proofs/singleton-beta.md` 第 1 节
（SHA-256 `fc6e2973b6a8334cb9b4a362b7424ff5421017d9c1b52eb0b56684d29f08644d`），在
一般 $K=2$ 词中 $\beta>0$ 当且仅当首次合并时低根分支 singleton。因此
$B_{\mathrm{low},1}$ 是该 $\beta>0$ 支持与 top-two 切片 $B$ 的交集。式 (2) 是
**计数**恒等式；它不使用 $\beta$ 的数值、不等同于全体 $A_2$ 上的 $\beta$ 支持或
$\mathbb E_2\beta$，也不提供 $A_2/B$ 的一致界。

## 独立有限核验（仅作实现交叉检查）

`audit_singleton_image_n3.py` 穷举 $P_3\square P_3$ 的全部排列，实际构造每个
$(\sigma,p)\mapsto\tau$ 并以 `Counter` 计其逆像；逐纤维断言该数等于 $R_p$，再独立
枚举全部 $B$ 并按首次合并分类。结果为

$$
 A_1=29568,\quad |B|=81984,\quad
 |\operatorname{Im}\Phi|=|B_{\mathrm{low},1}|=35264,
$$

$$
 B_{q=2}=14784,\qquad B_{\mathrm{s}}=55744
 =2\cdot35264-14784.
$$

此有限枚举不承担一般图证明；一般结论由前两节的双向构造承担。
