# 核接口与两项路线障碍：独立严格附审

2026-11-15（Asia/Shanghai）。对应父任务第9轮的一般上界核接口审查，本附审不另起研究轮次。审查工作区为独立 clone `pcw-kernel-audit`；没有改主仓、旧 T 或其他 Agent 文件，没有提交或推送。

**结论：两项障碍均成立，新增固定核族地板也已独立严格核验。** 从准确旧 K₄ 出发，只减原点样条的额外半积分收益严格小于 `5×10⁻⁶`。对固定 `p=47/40` 的奇异幂，只将生成元阶数改为任意 `p<a<2`，无论增加怎样的一个固定有限 C² 紧支撑样条组合，都会在轴旁违反 C4。后者是解析结论，不依赖有限网格。此外完整 1182 变量族的成本严格大于 `3.589342`。以上均不构成新的 C₂ 上界或 C₂ 下界；没有发现反例，但原稿中压缩的 O(1) 说明应保留下文的一致性与量词细节。

## 1. 准确对象与来源核

沿用实际旧检查器的归一化

\[
 A_a f(x)=\sum_{j=1}^2\int_{\mathbb R}
 [f(x+te_j)+f(x-te_j)-2f(x)]|t|^{-1-a}\,dt.
\]

这是符号为 `−(|ξ₁|ᵃ+|ξ₂|ᵃ)` 的生成元的正倍数。积分遍及整个实轴；把它误换成半轴会错一个因子 2。

我回读了旧 `round004/proofs/kernel.md`、准确 `round004/kernel/candidate.json`、旧上游 provenance 及实际 `hl_core.py`，而非以历史 PASS 摘要代替。准确 K₄ 是

\[
 K_4(u,v)=c\bigl((|u|+|v|)^{-p}-R^{-p}\bigr)_++S_4(u,v),
\quad p=47/40,\quad R=33/20,\quad h=1/32,
\]

其中 `c=154500107513046875/85907875725049856>0`。准确逐系数比较证明 K₄ 与旧 scaled 核 K₀ 的唯一区别是 `a["0,0"]` 减少 `1/200`。因此

\[
 P(u,v)=B(32u)B(32v),\qquad K_4=K_0-P/200,
 \qquad \tfrac12\int P=1/2048.
\]

这里 B 是积分为 1 的中心三次样条。旧输入与当前 Tomar 公开核的奇异系数并不相同，不能把二者视为同一证书。源文件 SHA256 详列在 [comparison.json](../computation/results/kernel-interface-review/comparison.json)；准确旧 candidate 哈希为 `88af41536e1117d953bbc4594ee99f1319b79b75d98d740c5d3c82caaf1f306f`，旧证明哈希为 `f6180f44d4dee187f90da2c0e208cb85da0470a18062fadeedd284586023efd2`。

旧完整 C4 记录仅作为历史实际运行来源；本附审没有重放其 443757 个叶，亦不重新认证整个旧一般上界。

## 2. 原点射线：严格小于 5×10⁻⁶

固定 `q=(35/64,5/256)`，考察 `K_4−δP`，`δ≥0`。q 位于正象限、`u+v<R`，避开原点、坐标轴及截断界面。对每个固定 δ，其生成元在 q 附近连续，并代表该处的分布。因此 C4 必须蕴含

\[
 0\le A_p(K_4-\delta P)(q)=A_pK_4(q)-\delta A_pP(q).
\]

`A_pP(q)>0` 也有直接几何理由：q 在 P 支撑之外，第二坐标仍位于样条正区间；水平跳跃命中 P 的正值区间，而所有跳跃值非负。

独立程序没有导入 `hl_core`，没有使用其已聚合的第四差分 `D_m`；它逐一展开张量基函数，使用三次样条的一维生成元公式分别求和，并独立写出奇异项闭式。Python 3.12.13、python-flint 0.9.0、192-bit Arb 的实际成功运行给出

\[
\begin{aligned}
 A_pK_4(q)&=0.0009179316247032744083392495348\ldots,\\
 A_pP(q)&=0.0928146348097201943608385332476\ldots.
\end{aligned}
\]

严格结论来自 JSON 中的准确有理向外端点，而非这两个打印小数。程序实际认证了较粗而充分的准确不等式

\[
 0<A_pK_4(q)<23/25000,\qquad A_pP(q)>9/100.
\]

故每个仍满足 C4 的射线成员都必须满足

\[
 \frac\delta{2048}
 \le\frac{A_pK_4(q)}{2048A_pP(q)}
 <\frac{23}{4608000}<\frac1{200000}.
\]

较紧球为 `4.829074660946747192388849…×10⁻⁶`。因此原稿的“额外收益不足 5×10⁻⁶”是严格的，方向和符号正确。这里 δ 是从 K₄ 再减去的量，不是从 K₀ 开始的总减量。不能由此声称 cap 端点可行，也不能排除联动改变其他系数的更大改进。

独立半积分也复得 `3.5908653129321965375303012628776…`。与作者输出中的五项准确区间均相交，见 comparison.json；球的不同包含关系是两种求和次序的舍入半径差异，没有被误记为相同区间。

## 3. 固定奇异幂、提高生成元阶数：完整轴旁证明

下述结论比特定 K₄ 更一般。固定

\[
 1<p<a<2,\quad c>0,\quad R>0,\quad
 F(u,v)=\bigl((|u|+|v|)^{-p}-R^{-p}\bigr)_+,
 \quad K=cF+S,
\]

其中 S 是**一个固定的** C² 紧支撑函数；有限个三次张量样条的任意有限实系数组合均属于此类。不需要 S 非负或有任何系数大小上限。

取任意紧区间 `I=[u₋,u₊]⊂(0,R)`，选 η>0 足够小，使 `4η<min(u₋,R−u₊)`。下述所有 O(1) 都对 `u∈I`、`0<v<η` 一致，常数允许依赖 `I,p,a,c,R,S`。

### 3.1 竖直负尖点与截断远尾

沿竖线设

\[
 f_u(s)=F(u,s),\quad b_u=-p u^{-p-1},\quad
 g_u(s)=f_u(s)-f_u(0)-b_u|s|.
\]

在 `|s|<2η` 内，g_u 为 C²，且 `|g_u''(s)|≤p(p+1)u₋^{-p-2}`。在全实轴上 `|g_u(s)|≤C(1+|s|)`，C 对 `u∈I` 一致。于是：小跳跃部分由二阶差分 `O(t²)` 控制；`|t|≥η` 的部分由线性增长控制，其尾积分是有限的，因为 `a>1`。因此 `A_a^v g_u(v)=O(1)`，已经包括截断点和远尾，不是在无截断核上形式删掉尾项。

直接积分给出

\[
 \int_{\mathbb R}(|v+t|+|v-t|-2v)|t|^{-1-a}dt
 =4\int_v^\infty(t-v)t^{-1-a}dt
 =\frac4{a(a-1)}v^{1-a}.
\]

故

\[
 A_a^v(cF)(u,v)
 =-\frac{4cp}{a(a-1)}u^{-p-1}v^{1-a}+O(1). \tag{V}
\]

### 3.2 水平正奇异项、角度系数与混合余项

水平小跳跃 `|t|<η` 时，`u±t` 离原点及截断界面有一致正距离，因此二阶差分部分有界。其余跳跃中，除 `u±t` 落在 `[-η,η]` 的两个区间外，核值有一致界，跳跃权也离 t=0 有正距离，所以总贡献 O(1)。两个原点附近区间合成

\[
 2c\int_{-\eta}^{\eta}
 \bigl((|s|+v)^{-p}-R^{-p}\bigr)|u-s|^{-1-a}ds+O(1).
\]

对 `|s|≤η`，一致展开

\[
 |u-s|^{-1-a}=u^{-1-a}+(1+a)u^{-2-a}s+O(s^2).
\]

奇项在对称区间准确抵消。剩余误差由
`∫₀^η s²(s+v)^(−p)ds≤∫₀^η s^(2−p)ds<∞` 控制；常数项 `R^(−p)` 亦只给 O(1)。因此

\[
 A_a^u(cF)(u,v)
 =\frac{4c}{p-1}u^{-1-a}v^{1-p}+O(1). \tag{H}
\]

这里两个跳跃方向及整个实轴积分都已计入，系数不是 2c。所有权展开与奇异项的乘积余项都已控制；由于生成元本身是两个坐标方向之和，不存在额外的交叉导数项。

将 (H) 的正首项除以 (V) 负首项的绝对值，得到准确的角度比例

\[
 \frac{a(a-1)}{p(p-1)}\left(\frac vu\right)^{a-p}\longrightarrow0.
\]

因此“负项主导”不是只比较 v 的幂却忽略 u 系数；它对上述 u 紧区间一致成立。

### 3.3 每个固定有限样条修补确实只有有界贡献

对任意 `S∈C_c²`，直接拆分 `|t|≤1` 与 `|t|>1` 即得全空间界

\[
 |A_aS(x)|\le
 \frac2{2-a}(\|S_{uu}\|_\infty+\|S_{vv}\|_\infty)
 +\frac{16}a\|S\|_\infty<\infty.
\]

有限组合、固定 h>0、每个有限系数保证右端有限。可以取极大的系数，但其贡献仍不能抵消趋于负无穷的轴旁主项。这里并未断言右端对所有样条系数组合有共同界；**所得到的负条带宽度允许依赖 S**。

### 3.4 从一致渐近到分布 C4 失败

合并 (V)、(H) 及 S，可取有限 C，使

\[
 A_aK(u,v)\le-N_0v^{1-a}+H_0v^{1-p}+C,
\]

其中

\[
 N_0=\frac{4cp}{a(a-1)}u_+^{-p-1}>0,
 \qquad H_0=\frac{4c}{p-1}u_-^{-1-a}.
\]

除以 `v^(1−a)>0` 后，右侧为
`−N₀+H₀v^(a−p)+Cv^(a−1)`，后两项趋于零。故存在 ε>0，使 `A_aK<0` 在开条带 `u∈int(I),0<v<ε` 上成立。

在此开集的任意紧子矩形上，核在小跳跃附近 C²，沿两坐标线与奇异原点都有正距离；小跳跃可用二阶界、其余跳跃可用可积的尾权统一控制。截断 Fubini 及支配收敛遂直接证明此处的点态积分代表分布 `A_aK`，并且连续。取支撑于该开集的非零非负光滑测试函数，其分布配对严格为负。因此分布 C4 失败，绝非仅在零测轴线失败。

这一步不借用只写在 `p=a` 核类中的原 Lemma 3.1；对失配指数所需的局部分布表示已在这里另证。

对于旧 `p=47/40` 与 `a=6/5`，独立 Fraction 运算实际核验

\[
 a-p=1/40,\quad -4p/(a(a-1))=-235/12,
 \quad 4/(p-1)=160/7.
\]

取 a=p 时两首项准确抵消。这说明所审稿件的符号与系数正确。结论仅排除保持 p 不变而提高 a 的方案；不排除同步改变奇异幂、生成元、非光滑补偿结构或核族。

## 4. 当前原文接口与量词

本附审实际打开了两位作者当前网页并通过页面链接读取 PDF，而非只复述作者 Agent 的摘录。具体支持位置如下。

| 当前原文 | 核对位置 | 结论 |
|---|---|---|
| [Tomar PDF](https://www.mefyi.com/papers/weak-type-squares/weak-type-squares.pdf) | 式(2.1)、Theorem 2.2、Lemma 3.1、§4、§5、Appendix A/B | 同一算子归一化与分布接口；有限样条的点态转分布；一般 L¹、全部正实尺度的传递。 |
| [Lin 文章](https://coolrmal.github.io/articles/two-partial-balayage-principles/) | Principle 2、Corollary 3.1、§4.3 | 同一共同半径中心正方形和半积分目标。 |
| [Tomar 页面](https://www.mefyi.com/papers/weak-type-squares/) | Abstract、Formal statement | 当前公开上界 3.59123；不能替代旧 K₄ 证书。 |

更具体地，Tomar PDF 的 Theorem 2.2（p.3）要求 K 偶、L¹、全空间非负、D 上≥1、原点外分布 `AK≥0`；Lemma 3.1 及其证明（pp.3–4）控制轴线/界面的可积主项。Appendix A（pp.9–10）用有理尺度的共同零测集及逼近传至全部尺度、全部 L¹。Appendix B（pp.10–11）明确区别旧“有限测度减原点原子”版本。§4（p.6）也说明指数失配的轴旁障碍。当前 Lin §4.3 的措辞与 Tomar 对较早 Lin 工作的描述有版本差异，应保留，不合并成一次已核验的形式化运行。本附审未构建 Lean。

菱形平均分母为 `2r²`；变换 `T(u,v)=(u+v,u−v)` 的 `|det T|=2` 同时作用于输入质量和接收面积，弱常数不变。核取 K≥1_D 时比较得到 `C₂≤(1/2)∫K`。这里没有额外丢失或重复一个 2，也没有改成偏心矩形、离散接收集或有限半径族。

若继续做固定 1182 参数核族的有限行对偶，认证的成本地板仅是该族充分证书方法的下界。无论其大于 3 或 3.5，均不得写成 C₂ 的下界，也不得写成一般核类已经最优。

## 5. 本次实际运行与交付

独立代码：[kernel_interface_review_check.py](../computation/kernel_interface_review_check.py)。主要输出：[independent-point.json](../computation/results/kernel-interface-review/independent-point.json)、[comparison.json](../computation/results/kernel-interface-review/comparison.json)。最终代码 SHA256 `29a214de0316a83328911ada082a7d5b925178844ec3e63baad6386ec905da9c`。最后一次实际运行约 0.097 秒、退出码 0。临时准确输入副本与日志位于 `/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-independent-audit/`。

复现命令（在本独立 clone 根目录）：

```sh
/Volumes/data/CodexProjects/T/planar-centered-square-weak11/computation/kernel-baseline/venv/bin/python -B \
  planar-centered-square-weak11/computation/kernel_interface_review_check.py \
  /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-independent-audit \
  planar-centered-square-weak11/computation/results/kernel-interface-review/independent-point.json
```

原作者材料给出的结果路径应读取为 `computation/results/kernel-interface/point-check.json`，而不是课题根下的 `results/kernel-interface/point-check.json`。这是定位修正，不影响数学结论。

本次实际验证：输入准确差分、192-bit 单点 C4 必要行、半积分、严格收益上限、有理渐近系数、来源量词；解析证明补全了失配指数的两方向一致余项和负开集。本次未重跑全域 C1/C2/C4、未求解浮点 LP、未验证完整 Lean 工程，也未核验备份；新增完整族地板的独立核验见下一节。提交哈希与推送核验不适用：父任务要求本子 Agent 仅交付独立文件。尚未解决的是一般输入锐上界、核方法的完备性及真正 C₂；两项路线障碍没有解决这些问题。

## 6. 新增完整固定族地板：独立全矩阵 Neumann 核验

本节为父任务追加的审查。作者随后提供 `floor-certified.json`，声称冻结 `p=α=47/40,R=33/20,h=1/32` 和准确旧输入的全部 1181 个样条系数类，令 `c≥0`、全部样条系数任意取实数时，每个满足连续 C1/C2/C4 的核都满足成本 `J=(1/2)∫K>3.589342`。这一**固定族**下界已经通过以下独立核验。

### 6.1 全变量覆盖、行必要性与弱对偶

独立标准库程序 [kernel_floor_scope_check.py](../computation/kernel_floor_scope_check.py) 没有导入原 LP、矩阵构造器或 hl_core，实际逐项认证：

- 证书 1181 个互异类恰好等于准确 K₄ 输入的全部键，没有删去弱响应列、零特征列或负系数列；加 c 后共 1182 变量。这里并未允许增加原键表以外的类。
- 对称轨道互不重叠，轨道大小分布为 `1×1+99×4+1081×8`。每个样条类的半积分成本准确是 `h²×轨道大小/2`；奇异系数列成本是 `pR^(2−p)/(2−p)`。
- 1182 个活动行互异，全部属于原 10657 点菜单，构成为 15 个 `K≥1` 行、465 个 `K≥0` 行、702 个 `AK≥0` 行。所有 K 行均不在原点，右端准确等于 `1_{u+v≤1}`；所有 A 行均有 `u≥v>0` 且 `u+v≠R`。整个10657点菜单也逐项通过了相同的准确位置/右端检查。

K 在原点外连续，故 C1/C2 的几乎处处形式也蕴含上述每个点的不等式，包括 D 的边界。A 点处积分连续并代表分布，因此分布 C4 蕴含这些点值非负。故全部都是原连续问题的必要条件；不必证明10657点菜单足以保证连续可行性。

写 `θ=(c,a₁,…,a₁₁₈₁)`、必要约束 `Gθ≥b`、成本 `mᵀθ`。令 A 为活动矩阵的转置：`A=G_activeᵀ`。若认证

\[
 y\ge0,\qquad Ay=m,
\]

将其他菜单行权重取准确 0，则对每个连续可行 θ，直接有

\[
 m^T\theta=y^TG_{\rm active}\theta\ge y^Tb.
\]

这里只用弱对偶，不使用浮点 LP 的最优状态、强对偶或一个未经验证的 primal 候选。

**负系数问题已检查：** 1181 个样条系数为自由实变量，这些列必须保持等式。仅将有理 midpoint 权重缩放到 `Gᵀy≤m` 不能给该完整族地板。实际认证的是全部 1182 个等式，c 列也用等式，故不存在这一漏洞。

### 6.2 独立重建，而非重跑原 solve

独立程序 [kernel_floor_neumann_review.py](../computation/kernel_floor_neumann_review.py) 已自带活动1182行的必要性预检：K/A标签、`u≥v≥0`、非原点、K右端 `1_{u+v≤1}`、A的 `v>0,u+v≠R,rhs=0`，并检查活动行互异。因此最小独立重放不依赖未归档的10657点全菜单；原scope checker与证据仍保留。该程序不导入 `hl_core` 或作者的 `kernel_family_floor.py`，不调用 solve/inverse、SciPy 或 NumPy。它从证书的准确有理点与 candidate 重建整个 `1182×1182` 的 A：样条项逐轨道、逐基函数直接展开；奇异项的菱形内/外闭式在此程序中独立写出。所有超几何整数差提示逐式核对，外侧第一原函数用 `a−c=−1`，其余相关原函数用 `b−c=−1`。

预条件 R 的来源可以是不可信浮点运算；独立程序将 `floor-preconditioner.bin` 每个 little-endian binary64 项解释成**准确 dyadic 有理数**。y₀ 是 `floor-certified.json` 每个准确有理包围端点的算术平均，不取打印小数。随后只用 Arb 的向外矩阵乘法。

具体而言，设 `E=I−RA`，认证 `q≥‖E‖∞` 且 q<1，再认证

\[
 \eta\ge\|R(m-Ay_0)\|_\infty,\qquad
 \delta\ge\frac\eta{1-q}.
\]

由 Neumann 级数，RA 可逆，故 A 和 R 均可逆。存在唯一真实解 `Ay*=m`，且

\[
 \|y^*-y_0\|_\infty\le\delta.
\]

因此 `min_i(y₀)_i>δ` 保证真实 y* 逐项严格正。本活动集恰有15个 rhs=1，故 `bᵀy*≥bᵀy₀−15δ`。这证明了等式真解的存在和正性；不是只检查一个跨零残差。

### 6.3 实际严格结果与复现边界

Python 3.12.13、python-flint 0.9.0、192-bit 的独立运行实际成功，最终最小包重放重建约9.129秒，总计约15.324秒、退出码0。保存的准确有理端点认证

\[
\begin{aligned}
 q&\le0.00075005778788188737027228764268923\ldots<1,\\
 \eta&\le5.2702019018594278531\ldots\times10^{-47},\\
 \delta&\le5.2741578250105952466\ldots\times10^{-47},\\
 \min_i(y_0)_i&>\delta,\\
 b^Ty^*&\ge3.58934235879042629824443857909495788423660003729\ldots\\
 &>1794671/500000=3.589342.
\end{aligned}
\]

这里小数只展示量级；严格方向由 [floor-neumann-independent.json](../computation/results/kernel-interface-review/floor-neumann-independent.json) 中向外有理端点承担。完整几何/类检查输出为 [floor-scope.json](../computation/results/kernel-interface-review/floor-scope.json)。Neumann 独立代码 SHA256 为 `a478c0dc81836cadabb55543985da780ff82980fe19f6542c42bbc35e1c173c4`。

此次输入绑定：`floor-certified.json` SHA256 `3264f2e5ad8690fb20651651be4dd19e089d80e6e50339de30d76dba403d0173`；R 二进制 SHA256 `2c6e75a224177c81875f78dce3d5d314f25e1da03a1e6eb518c99821f87f2dfc`。此处R的归档问题已解决：作者用可逆 byte-shuffle8+gzip 保存为 `floor-preconditioner.byteshuffle8.gz`（9,822,059字节，SHA256 `6765a557999b826e7910f14cfd9be485b1ce3864544f2a10869a66d61952e937`）。本附审实际运行归档 `unpack_kernel_preconditioner.py`，经压缩/原始双哈希核对恢复11,176,992字节R到独占CodexTemp，再从归档candidate与floor-certified重放本独立checker。最新代码相对首次数值核验仅增加scope断言与记录，数学矩阵及Neumann公式未变。R的数值近似正确性无需信任，任意R都必须经过上述乘法验证才可接受。

本结果证明冻结族的成本地板。它不证明一般核类成本地板，不证明此核法对锐常数完备，也不把 `3.589342` 变成 C₂ 下界。数学范围内未发现阻断缺口；一般 C₂ 仍开放。

最终最小包重放只需归档candidate、floor-certified、压缩R及其storage元数据、解码脚本和本独立checker，不读取旧T的数据或完整LP菜单。复现示例（用已有含python-flint的解释器运行）：

```sh
python -B computation/unpack_kernel_preconditioner.py \
  computation/results/kernel-interface/floor-preconditioner.byteshuffle8.gz \
  /Volumes/data/CodexTemp/T2/kernel-floor-review-R.bin
python -B computation/kernel_floor_neumann_review.py \
  computation/results/kernel-interface/floor-certified.json \
  computation/inputs/kernel-audit/candidate.json \
  /Volumes/data/CodexTemp/T2/kernel-floor-review-R.bin \
  computation/results/kernel-interface-review/floor-neumann-independent.json
```
