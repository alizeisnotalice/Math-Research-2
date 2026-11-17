# 第37轮：低根小支的有序纤维比较及其精确损失

## 目标和记号

固定一个有限连通简单图，按降序揭示顶点。对一个恰有两个峰的词，第一次合并前
全局高根和低根所在连通支的大小分别记为 \(a,b\)，并令 \(q=a+b\)。这里
“高/低”始终指两峰标签的高低，而不是以固定顶点编号作的规范定向。

对 \(b\ge1\)，令 \(\mathsf L_b\) 为低根支恰有 \(b\) 个顶点的事件，令
\(\mathsf H_b\) 为高根支恰有 \(b\) 个顶点的事件。它们在 \(a=b\) 的纤维可以同时
发生，这没有歧义。下文的 \(F\) 只依赖无序根对和 \(q\)，从而在固定纤维的交错
中不变；可取 \(F=1,\rho,\|r_1-r_2\|_2^2/n^2\)，或乘以对称距离事件。

目标是控制方格上
\[
 \mathbb E_2[\rho\,\mathbf1_{\{b\le k\}}]
 \quad\hbox{或}\quad
 \mathbb E_2[E\,\mathbf1_{\{b\le k\}}],
\]
其中 \(b\) 此处是低根支大小。第36轮 已控制相应的高根支事件；以下核查这种控制
能否经可逆/分数插入转回低根侧。

## 精确方向变换

在固定无向首次合并资料中，固定两支内部词、合并点和后缀，支大小为 \(a,b\)。
若全局高根在大小为 \(a\) 的支，完整 \(K=2\) 词的交错数为
\[
 \binom{q-1}{a-1}=\binom{q-1}{b}.
\]
若交换两根方向、令大小为 \(b\) 的支成为全局高根，则交错数为
\[
 \binom{q-1}{b-1}.
\]
两者之比是 \(a/b\)。因此，对每个整数 \(b\ge1\)，有严格的有限图恒等式
\[
 \boxed{\quad
 \mathbb E_2[F\,\mathbf1_{\mathsf L_b}]
 =\mathbb E_2\!\left[
       \frac{q-b}{b}\,F\,\mathbf1_{\mathsf H_b}\right].
 \quad} \tag{1}
\]

这里右边的 \(q-b\) 正是高支大小为 \(b\) 时的低支大小。等式也覆盖 \(a=b\)：
系数为 \(1\)，且两事件都在该纤维中发生。对 \(1\le b\le k\) 求和得
\[
 \mathbb E_2[F\,\mathbf1_{\{b_{\rm low}\le k\}}]
 =\mathbb E_2\!\left[
       \frac{q-b_{\rm high}}{b_{\rm high}}\,
       F\,\mathbf1_{\{b_{\rm high}\le k\}}\right]. \tag{2}
\]

这是“整支删除—改变根方向—再分数插入”给出的自然无损账本：一侧每个交错
对应另一侧，计重系数就是交错计数的精确比值。没有使用混合、渐近或图的特殊性。

当 \(b=1\)，(1) 化为
\[
 \mathbb E_2[F\,\mathbf1_{\mathsf H_1}]
 =\mathbb E_2\!\left[\frac{F}{q-1}\mathbf1_{\mathsf L_1}\right], \tag{3}
\]
这与 第36轮 singleton-geometry 的 high/low-singleton 恒等式一致。并且 第36轮 singleton-geometry 的几何回插式严格给出
\[
 \mathbb E_2[E\,\mathbf1_{\mathsf L_1}]
 =\frac{A_1}{A_2}\,\mathbb E_1\Psi_E. \tag{4}
\]
故 \(b=1\) 的问题完全等价于控制这个已准确归一化的一峰加权量；(3) 只控制它
额外乘 \(1/(q-1)\) 的版本。

已证的 early-singleton 计数仍给出一个有限但真实的原律切分：对
\(2\le K\le N-1\)，方格 \(n\ge4\) 上
\[
 \mathbb P_2(\mathsf L_1,\ q\le K)
 \le \frac{16K(K-1)}{N},
 \qquad
 \mathbb E_2[E\,\mathbf1_{\mathsf L_1}\mathbf1_{\{q\le K\}}]
 \le\frac{32K(K-1)}{N}. \tag{5}
\]
第二式仅用 \(0\le E\le2\)。所以任何 \(K=o(n)\) 时，singleton 的早合并远峰
质量确实消失；未控制的正是 \(q>K\) 的 late-singleton 完成数加权质量。

## 与高根小支尾界的唯一可得传递

第36轮 已证明，对任意 \(k\ge1\)，
\[
 \mathbb E_2[\rho\,\mathbf1_{\{b_{\rm high}\le k\}}]\le2k,
 \qquad
 \mathbb E_2[E\,\mathbf1_{\{b_{\rm high}\le k\}}]
 \le\frac{2k(n-1)}{n^2}. \tag{6}
\]
直接代入 (2) 会产生因子 \((q-b)/b\)，它可达 \(N-1\)。因此只能得到
\[
 \mathbb E_2[\rho\,\mathbf1_{\{b_{\rm low}\le k\}}]
 \le2k(N-1), \tag{7}
\]
这甚至弱于平凡量级界，且不能导出固定比例远峰概率的衰减。

也可把 (2) 写成确实有界的“分数”版本：
\[
 \mathbb E_2\!\left[
   \frac{\rho}{q-b_{\rm low}}\mathbf1_{\{b_{\rm low}\le k\}}\right]
 \le2k. \tag{8}
\]
证明是逐 \(b\) 用 (1)，右侧成为
\(\mathbb E_2[(\rho/b)\mathbf1_{\mathsf H_b}]\)，再以
\(\rho\le q\) 或直接逐纤维比较
\[
 \frac{a}{q}\frac{\rho}{a}
 =\frac{\rho}{q}
 \le\frac{b}{q}\frac{\rho}{b}.
\]
对 \(b\le k\) 求和后，因 \(1/b\le1\)，第36轮 的整体高支界直接给出右侧至多
\(2k\)。但在 \(q\asymp N\)、\(\rho\asymp n\)
的低根 singleton 纤维上，左端每词权重仅为 \(O(1/n)\)，所以 (7) 不控制所需的
无权 \(\rho\) 或 \(E\) 质量。

## 精确失败点：大交错比

对于 \(b=1\)，第36轮 singleton-geometry 的显式方格词族（\(n\ge5\)）有
\[
 a=N-5,\quad q=N-4,\quad \rho=2n-2,\quad
 E=\frac{2(n-1)^2}{n^2}.
\tag{9}
\]
在该固定首次合并纤维中，低-singleton 方向有 \(a=N-5\) 个交错，而
high-singleton 方向仅有一个；这正好实现 (1) 的系数 \(a/b=N-5\)。因此在这个
纤维层级，由 high-singleton 的 \(O(n^{-1})\) 空间矩界反推 low-singleton 的空间矩
会损失一个 \(\Theta(N)\) 因子；这个损失不是松估计或某个非可逆删除造成的。

该族本身没有给出相对 \(A_2\) 的正质量，故不是原题反例。它只严格排除了下述路线：
仅凭 第36轮 的 high-small 尾界代入保持首次合并纤维的根交换/整支插入式 (2)，
得不到 \(\mathbb E_2[E\mathbf1_{\{b_{\rm low}\le k\}}]=o(1)\)。

结合 (4)，\(b=1\) 的第一未证步骤可准确表述为
\[
 \frac{A_1}{A_2}\,\mathbb E_1\Psi_E\longrightarrow0. \tag{10}
\]
它需要新的完成数或几何加权估计；高/低方向对称、分数插入和 第36轮 的 high-small
控制均不能替代它。

## 有限整数/有理审计

脚本 <code>audit_low_small_orientation_n3.py</code> 完整枚举 \(P_3\square P_3\) 的
\(9!\) 个排列，直接识别两根、第一次合并前两支和 \(q\)。它对每一个实际出现的精确
支大小 \(b\)，以 \(F=1,\rho,E\) 核对 (1)。这是定义、方向和 \(a=b\) 边界的有限
审计；(1)--(4) 的一般性来自上述交错计数。

## 状态

严格得到 (1)--(5) 及由已证 第36轮 导出的弱分数传递 (8)。严格定位的失败是 (9)
中的不受控交错比。没有证明低根小支的远峰质量趋零，也没有证明或否定原题。
