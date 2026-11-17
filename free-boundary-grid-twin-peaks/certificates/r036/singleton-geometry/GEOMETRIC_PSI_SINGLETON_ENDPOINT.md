# R042：几何加权回插与 singleton 端点乘积质量

## 精确几何回插恒等式

对一峰降序词 $\sigma$，根为 $r$，对 R038 的 eligible 顶点 $p$ 定义

$$
 e(r,p):=\frac{\|r-p\|_2^2}{n^2},\qquad
 \Psi_E(\sigma):=
 \sum_{p\notin\{r\}\cup N(r)}\frac{l_p-1}{R_p}\,e(r,p). \tag{1}
$$

令 $B_+^E$ 是所有恰两峰且 $\beta>0$ 的标号的 $e(r_1,r_2)$ 总和，其中
$r_1$ 是全局高根、$r_2$ 是另一峰。则对每个有限连通图（在方格时使用上式的归一化），

$$
 \boxed{\quad A_1\mathbb E_1\Psi_E=B_+^E.\quad} \tag{2}
$$

证明与 R038 的标量纤维相同，但根位置不随稳定移动或两支交错而变。固定一个
top-two 的 low-singleton 输出 $\tau$，合并前高支大小为 $q-1$。其所有 $R_p$ 个
一峰逆像都有相同的根对 $(r,p)$ 和 $l_p=q$，所以该提升纤维对左端贡献

$$
 R_p\frac{q-1}{R_p}e(r,p)=(q-1)e(r,p).
$$

而同一固定首次合并数据恰有 $q-1$ 个一般 $K=2$ 的 low-singleton、即
$\beta>0$ 输出，根对不变。逐纤维相加即得 (2)。这不是只对标量 $B_+$ 的替换，
而是保留原图平方欧氏根距的精确恒等式。

## 主定理：high/low singleton 的精确方向恒等式

本节只处理支大小恰为 $1$ 的方向比较，不涉及另行研究的一般小支
\(a\le k\)。在完整 \(K=2\) 条件律下，令 \(\mathsf H_{\rm low}\) 表示首次合并时
**低根**所在支为 singleton 的事件；按既有 singleton--\(\beta\) 刻画，
\(\mathsf H_{\rm low}=\{\beta>0\}\)。令 \(\mathsf H_{\rm high}\) 表示同一时刻
**全局高根**所在支为 singleton 的事件。

对任何只依赖 \(q\) 和无序根对的非负对称可观测量 \(F\)（特别是
\(1,\rho,E,\mathbf1_{\{\rho>\varepsilon n\}}\) 及它们的乘积），有精确有限图恒等式

$$
 \boxed{\quad
 \mathbb E_2\!\left[F\mathbf1_{\mathsf H_{\rm high}}\right]
 =\mathbb E_2\!\left[\frac{F}{q-1}\mathbf1_{\mathsf H_{\rm low}}\right].
 \quad} \tag{H1}
$$

**证明。** 固定两支内部次序、首次合并数据和公共后缀。在恰一支为 singleton、
合并前总大小为 \(q\) 的纤维中，singleton 在一处交错位置时为高支，在其余
\(q-1\) 处时为低支。无序根对和 \(q\) 不随此交错改变，故高支的贡献 \(F\) 正好
等于低支的总贡献 \((q-1)(F/(q-1))\)。若 \(q=2\)，双方均为 singleton；两个根
方向在两边各出现一次，且 \(q-1=1\)。故该边界纤维也满足同一等式。对所有纤维
求和并除以 \(A_2\)，即得 (H1)。\(\square\)

令 \(d=\rho(r_1,r_2)\)。两条合并前支和合并顶点给出一条根间路径，故 \(d\le q\)；
不同峰不相邻，故 \(d\ge2\)。在方格中，若两根坐标差为 \((s,t)\)，则

$$
 E=\frac{s^2+t^2}{n^2}
 \le\frac{(n-1)(s+t)}{n^2}
 =\frac{(n-1)d}{n^2}. \tag{H2}
$$

在 (H1) 中分别取 \(F=d,E\)，得到完全不需 UI 假设的原条件律结论

$$
 \mathbb E_2[d\mathbf1_{\mathsf H_{\rm high}}]\le2,\qquad
 \boxed{\quad
 \mathbb E_2[E\mathbf1_{\mathsf H_{\rm high}}]
 \le\frac{2(n-1)}{n^2}.\quad} \tag{H3}
$$

又若 \(\varepsilon n\ge2\)，则 \(d>\varepsilon n\) 时整数关系
\(q-1\ge d-1\ge\lfloor\varepsilon n\rfloor\) 和 (H1) 给出

$$
 \boxed{\quad
 \mathbb P_2(\mathsf H_{\rm high},\rho>\varepsilon n)
 \le\frac1{\lfloor\varepsilon n\rfloor}.\quad} \tag{H4}
$$

双方 singleton 的 \(q=2\) 扇区只有 \(d=2\)，故对固定正比例的 far event 在大
\(n\) 时无贡献。(H3)--(H4) 是 high-root singleton 原双峰扇区的真实空间控制；
它们不控制 low-root singleton，亦不控制非 singleton endpoint 支。

## 与 joint kernel 中的 singleton 端点的精确比较

沿 `JOINT_GEOMETRY_KERNEL.md` 的记号，令 $\lambda_n$ 是 top-two $B$-law，
$Z_n=\int\kappa\,d\lambda_n=A_2/B$。令 $\mathsf S$ 表示首次合并时至少一支为
singleton 的 top-two 词，并令

$$
 X_+:=B_+^E=A_1\mathbb E_1\Psi_E,\qquad
 X_{\mathsf S}:=B\int_{\mathsf S}E\kappa\,d\lambda_n. \tag{3}
$$

把第一次合并纤维分成“恰一 singleton”和“双方 singleton ($q=2$)”两类。对恰一
singleton 的固定纤维，设去除交错的公共完成数为 $P$、高支大小为 $q-1$、几何为
$E$。则

$$
 X_+\text{ 的该纤维贡献}=(q-1)PE,\qquad
 X_{\mathsf S}\text{ 的该纤维贡献}=qPE. \tag{4}
$$

后式来自两个 top-two 根方向各出现一次，且各自
$\kappa=m=q/2$。对 $q=2$ 的双方 singleton 纤维，两边贡献均为 $2PE$。因此

$$
 X_{\mathsf S}=X_++C_E,\qquad
 0\le C_E\le\frac12X_+, \tag{5}
$$

其中 $C_E=\sum_{\text{恰一 singleton fibres}}PE$。特别地，

$$
 \boxed{\quad
 \frac{A_1\mathbb E_1\Psi_E}{A_2}
 \le \frac{\int_{\mathsf S}E\kappa\,d\lambda_n}{Z_n}
 \le\frac32\frac{A_1\mathbb E_1\Psi_E}{A_2}.
 \quad} \tag{6}
$$

左端正是原 $K=2$ 条件律下

$$
 \frac{B_+^E}{A_2}=\mathbb E_2\bigl[E\,\mathbf1_{\{\beta>0\}}\bigr]. \tag{7}
$$

所以 singleton 端点的 $E\kappa/Z_n$ 质量趋零，当且仅当

$$
 \mathbb E_2\bigl[E\,\mathbf1_{\{\beta>0\}}\bigr]\longrightarrow0. \tag{8}
$$

这把一个 endpoint product 问题精确地转成原条件律下的、带空间因子的
$\beta$-正支持质量；没有遗漏 singleton 的 $\kappa=q/2$ 放大。

若只取真正端点 $C_\delta=\{\Theta<\delta\}\cup\{\Theta>1-\delta\}$，对 singleton
纤维其条件正是 $q>1/\delta$（端点等号可按闭开约定处理）。因此 (4)--(6) 同样逐
$q$ 截断成立：大 $q$ singleton 端点质量与 $B_+$ 中同样大 $q$ 的几何质量在
$[1,3/2]$ 的绝对因子内等价。

## 这没有自动控制全部 endpoint

式 (6) 只覆盖一支大小为 $1$ 的端点。$\Theta<\delta$ 还包括
$2\le a\le\delta q$ 的非 singleton 小支；当前回插恒等式没有把这些纤维归约为
$B_+$。所以 (8) 是 singleton endpoint 的严格等价判据，不是 r041 式 (17) 的完整证明。

并且现有结果未证明 (8)。有限 $n$ 的 $\beta>0$ 支持计数不是渐近上界，且
$B_+\le A_2$ 只能给平凡界。不能以标量 $\mathbb E_1\Psi=B_+/A_1$ 的任何增长或
衰减结论代替带 $E$ 的比值 (7)。

## 一个具体的点态失败：端点 singleton 可以远

不存在“在 endpoint union 上 $E=o(\min\{\Theta,1-\Theta\})$”的逐词界。对 $n\ge5$，取

$$
 v=(n-1,n-1),\quad x=(n-2,n-1),\quad y=(n-1,n-2),
$$

并如既有 `proofs/singleton-beta.md` 第 4 节，令高支为除
$v,x,y,(n-2,n-2),(n-1,n-3)$ 外的连通集。取从 $r=(0,0)$ 开始的高支连通次序，
构造词时先写 $r,v$，再按该次序写其余高支点，最后写
$x,y,(n-2,n-2),(n-1,n-3)$。这是一个
top-two、low-singleton 首次合并词，且

$$
 q=N-4,\qquad \Theta=\frac{N-5}{N-4},\qquad
 \kappa=\frac{N-4}{2},\qquad
 E=\frac{2(n-1)^2}{n^2}\longrightarrow2. \tag{9}
$$

所以 $\min\{\Theta,1-\Theta\}=1/(N-4)$，并且
$E/\min\{\Theta,1-\Theta\}\asymp N$。该族只说明 endpoint union 上的点态局部化不可能；它没有给出相对 $A_2$ 或
$Z_n$ 的正质量，故不是原题反例，也不否定 (8)。

`audit_far_singleton_family.py` 对 $n=5,6,8,16$ 实际生成此词，逐点直接判定恰有
两个峰、首次合并顶点为 $x$，并核对合并前两支大小为 $(N-5,1)$；它只复核这一显式
构造，不给出任何计数或渐近质量结论。

## $3\times3$ 独立核验（不承担一般证明）

`audit_geometric_psi_n3.py` 完整枚举全部排列，独立计算 (1) 的有理和、直接筛选
$\beta>0$ 双峰词并计算几何和、再计算 top-two singleton 的原始 $E\kappa$ 和。结果
为

$$
 B_+^E=\sum\Psi_E=\frac{322784}{9},\qquad
 X_{\mathsf S}=\frac{136768}{3},
$$

并逐项通过 (2)、(5) 的 $[1,3/2]$ 界。有限枚举只核对定义、方向和 $q=2$ 边界；
一般结论由上述纤维计数承担。

<code>audit_high_singleton_identity_n3.py</code> 另以不同于回插的直接枚举方式，
逐个双峰词识别全局高根/低根的首次合并支。它对
\(F=1,d,E,\mathbf1_{\{d\ge3\}}\) 分别核对 (H1)：例如 high-singleton 原始词数为
\(35264\)，而 low-singleton 的 \((q-1)^{-1}\) 加权数亦为 \(35264\)；两边的
\(d\) 和均为 \(83808\)，\(E\) 和均为 \(125344/9\)。这只是 \(3\times3\) 的独立
定义和边界审计，不承担一般证明。
