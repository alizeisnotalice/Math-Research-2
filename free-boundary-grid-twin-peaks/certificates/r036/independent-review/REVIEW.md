# r036 独立数学审查

审查日期：2026-11-17（Asia/Shanghai）。本记录只读审查 staged r036 的三份解析材料；没有运行枚举、修改 staged 文件、SSH 或启动重计算。

审查对象及当次 SHA-256：

| 文件 | SHA-256 |
| --- | --- |
| `ordered-branch/ORDERED_BRANCH_KERNEL.md` | `1ed51fc95f12c97c0bb64200aa71bf71b064b4add48d56bd48ebad1512fe79c7` |
| `singleton-geometry/GEOMETRIC_PSI_SINGLETON_ENDPOINT.md` | `924e4fd34e5f008efed5d8dbdf60b4714188e36b1743304dde2ddba6fe3453ba` |
| `early-low-branch/BETA_SUPPORT_ANALYSIS.md` | `b3bf3c95716d8b8e20cf107d710fb69bd24e83dc07a62af5719074e9aa4aeb82` |
| `singleton-geometry/audit_high_singleton_identity_n3.py` | `c36825610bb17cfe35e711b3bedd55166bf779c3c7aee754c9f98c61ea06af56` |

## 判定

**解析部分通过本次审查；没有发现推翻所列有限图恒等式或不等式的数学漏洞。** 该判定不证明原 twin-peaks 问题，也不补足 low-small endpoint 或 balanced interior 的空间加权质量控制。

初读时的旧措辞“几何加权结论既不能上界 (7)，也不能由 (7) 得到”不正确：因为 0≤E≤2，确有单向蕴含
\[
\mathbb P_2(\beta>0)\to0
\quad\Longrightarrow\quad
\mathbb E_2[E\mathbf1_{\{\beta>0\}}]\to0.
\]
反向不成立，原因正是该支持上没有 E 的统一正下界。现 staged §2 已明确写出这项充分而非必要的单向关系，且不再含旧句；因此它不再是条件性阻塞。该修正不影响 §1 的早秩计数或 §3 的局部反例。

## 1. 有序首次合并核：通过

固定未定向纤维，支大小为 (a,b)，(q=a+b)，并指定全局高根在 (A)。完整双峰交错数与 top-two 交错数分别为
\[
\binom{q-1}{a-1},\qquad\binom{q-2}{a-1}.
\]
它们之比严格为
\[
\kappa_{\rm ord}=\frac{q-1}{b}.
\]
这里分母是低根支大小 (b)，而不是对称核，也不是 (a)；文件作了正确区分。原 (K=2) 条件律中高根落在 (A) 的概率为 (a/q)，所以
\(\mathbb E_2(q/a)=2\) 是逐纤维的两个方向各贡献 1 的恒等式。该恒等式没有被误当作低根小支或内部区的控制。

对根交换不变的 (H(D,E,U))，两个有序方向平均确实给
\[
\frac12\left(\frac{q-1}{a}+\frac{q-1}{b}\right)
=\frac{q(q-1)}{2ab},
\]
故仅在这种对称可观测量上可回到 r035 的对称核。这个量词限制是必要的，文稿已写明。

若 (0<\delta<1/2)，每一纤维至多一个支满足大小比例不超过 \(\delta\)；高根进入该支的条件概率就是该比例，故
\(\mathbb P_2(\Theta\le\delta)\le\delta\)。再用 (0\le E<2) 即给相应空间矩端点界。若高根支大小至多 (k)，两支和合并点连通给出 \(\rho\le q\)；逐方向平均为
\[
\rho\left[\frac aq\mathbf1_{a\le k}+\frac bq\mathbf1_{b\le k}\right]
\le a\mathbf1_{a\le k}+b\mathbf1_{b\le k}\le2k.
\]
因此 \(\mathbb E_2[\rho\mathbf1_{\{\mathrm{high\ size}\le k\}}]\le2k\) 的常数、方向和量词均正确。

## 2. 几何回插与 high/low singleton：通过

几何权 \(e(r,p)\) 在稳定移动、支内次序及交错中保持根对不变，因此标量回插纤维式可逐项乘该权。每个低-singleton top-two 输出的 (R_p) 个逆像共有 \(q-1\) 倍的加权贡献，与一般 (K=2) 的低-singleton 插入位数完全相符，得到
\(A_1\mathbb E_1\Psi_E=B_+^E\)。

形式措辞需注意：当前 (1) 的 (e(r,p)) 是方格的欧氏权，故“每个有限连通图”应理解为“方格，或任意有限图上预先指定、仅依赖根对且在该纤维中保持不变的权”。这不影响方格结论；若保留全图说法，建议显式加入该抽象权的定义。

在“恰一 singleton”纤维中，singleton 占据第一交错位置时为高根支，其余 (q-1) 个位置为低根支，所以对任何只依赖 (q) 和无序根对的非负 (F)，
\[
 \mathbb E_2[F\mathbf1_{\mathsf H_{\rm high}}]
 =\mathbb E_2\!\left[\frac F{q-1}\mathbf1_{\mathsf H_{\rm low}}\right].
\]
当 (q=2) 时每个方向同时 high/low singleton，且 (q-1=1)，故没有漏项或需除以二的例外。以 \(d\le q\) 代入得到 \(d/(q-1)\le2\)，结合方格不等式 \(E\le(n-1)d/n^2\)，文稿的 (H3) 和 (H4) 均正确。

对于 top-two singleton endpoint 的 \(E\kappa\) 质量，恰一 singleton 纤维中
\(X_+=(q-1)PE\)、\(X_{\mathsf S}=qPE\)。其余项 (C_E=PE\) 且此时 (q\ge3)，所以逐纤维 \(C_E\le X_+/2\)。双方 singleton 的 (q=2) 纤维在两边均为 (2PE)。故
\[
X_+\le X_{\mathsf S}\le\tfrac32X_+
\]
及其除以 \(A_2=BZ_n\) 后的版本成立。这是 singleton endpoint 的精确比较，不延伸到支大小 \(2,\dots,\lfloor\delta q\rfloor\)；文稿没有越界声称。

## 3. 早 singleton 合并与局部动作：通过

低 singleton 的两步双射中高支大小为 \(a\)、\(q=a+1\)，给定一峰词和 \(a\) 有恰好 \(a\) 个插入槽。因此限制 \(q\le K\) 后，逐一峰词至多
\[
\sum_{a=1}^{K-1}a=K(K-1)/2
\]
个像，故 \(B_{+,\le K}\le K(K-1)A_1/2\)。再用既有方格输入 \(A_2/A_1\ge N/32\) 得到常数 16；范围 \(2\le K\le N-1\) 与代数相符。该论证不涉及 late singleton，文稿也正确保留该缺口。

给出的 \(3\times3\) 词中，交换角点 singleton 的另一个邻点 (y) 至 (w) 前时，(y) 同时接触 singleton 与已形成高支，因而仍是首次合并而非低支扩张。它只否定该逐词紧前交换策略，不否定更复杂的完成数加权注入；该范围表述正确。

## 4. 证据和范围记录

* 当前 r036 目录中存在 `audit_n3_ordered_kernel.py`、`audit_geometric_psi_n3.py`、`audit_far_singleton_family.py` 及 `audit_corner_swap_n3.py`；它们可作为有限定义/边界的交叉核验，但不承担一般纤维证明。
* `audit_high_singleton_identity_n3.py` 现已加入 staged r036，并写出 `audit-high-singleton-identity-n3.json`。该实际 3×3 输出记录 high 词数与加权 low 词数均为 35264，距离和均为 83808，欧氏平方和均为 125344/9，且 `all_identities_passed=true`；故先前的证据文件缺失问题已消除。它仍只是有限定义/方向的交叉核验，不承担 (H1) 的一般证明。
* 已审材料明确留下 balanced interior 以及 low-root-small/late-singleton 的加权空间质量未控；本审查没有把高根小支消失、有限 n 审计或 singleton 比较外推为 \(Q_n\to0\) 或原远峰概率结论。
