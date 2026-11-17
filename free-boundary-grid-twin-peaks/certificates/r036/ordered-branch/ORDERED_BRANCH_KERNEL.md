# r047：有序高根支的精确核与单侧端点消去

## 结论与范围

本文件修正首次合并比例变量的一个必要区分。若记录全局高根所在的支，则从 top-two
\(B\)-law 到完整 \(K=2\) 条件律的核不是对称核
\(m=q(q-1)/(2ab)\)，而是有序核 \((q-1)/b\)。对根交换不变的量（例如峰距与平方
欧氏距离），两者积分相同；但对“高根支很小”之类非对称事件，必须使用有序核。

严格结果是：完整 \(K=2\) 条件律下高根支占比的零端点自动可忽略；比例端点中只剩
低根支很小的一侧。除此以外，平衡内部纤维仍可有宏观根距且其总加权空间几何也未控制。
本文件不解决任一剩余问题，因而不解决原题。

## 1. 有序首次合并纤维

固定一个未记录根方向的首次合并纤维：两支 \(A,B\)（大小 \(a,b\ge1\)）、各自内部
一峰词、合并点和合并后缀均已固定，记 \(q=a+b\)。讨论全局高根落在固定分支 \(A\) 的方向时，定义
\[
\Theta=\frac aq,\qquad U=\frac1q.
\tag{1}
\]
峰距及平方欧氏根距记为
\[
D=\frac{\rho}{n},\qquad E=\frac{\|r_1-r_2\|_2^2}{n^2}.
\tag{2}
\]
它们在这个纤维的所有交错及根方向间不变。

完整双峰词的交错总数为 \(\binom qa\)。全局高根在 \(A\) 的交错数为
\(\binom{q-1}{a-1}\)，故条件概率为 \(a/q\)。在 top-two 扇区，根必须占前两位；
高根在 \(A\) 或 \(B\) 的交错数各为 \(\binom{q-2}{a-1}\)，所以两方向恰各占一半。

于是，对有序 \(B\)-law，固定纤维且高根在 \(A\) 时，从 \(B\) 到完整 \(K=2\) 的
计数比严格为
\[
\kappa_{\rm ord}
=\frac{\binom{q-1}{a-1}}{\binom{q-2}{a-1}}
=\frac{q-1}{b}
=\frac{1-U}{1-\Theta}.
\tag{3}
\]
令 \(Z_n=A_2(n)/B(n)\)。若 \(\lambda_n^{\rm ord}\) 是
\((\Theta,D,E,U)\) 的有序 \(B\)-law，则对任意有界可测 \(F\)，
\[
\boxed{\quad
\mathbb E_2F(\Theta,D,E,U)
=\frac1{Z_n}\int F\,\kappa_{\rm ord}\,d\lambda_n^{\rm ord}.
\quad}
\tag{4}
\]

这包括 \(a=b=1,q=2\)：此时 \(\kappa_{\rm ord}=1\)。若 \(b=1,q>2\)，则
\(\kappa_{\rm ord}=q-1\)，亦无例外。

## 2. \(L\)-提升的有序条件均值

在高根 \(A\)、低根 \(B\) 的方向，低根首个支内邻点等待的精确交错计算给出
\[
\mathbb E_B[L\mid\mathcal D,\ {\rm high}=A]
=\frac{q-1}{b}
=\kappa_{\rm ord}.
\tag{5}
\]
若 \(b=1\)，低根直到合并点才有邻点，所以 \(L=q-1\) 确认 (5)；若
\(a=b=1\)，\(L=1\)。因而 (4) 也直接来自
\(\mathbb E_B[LH\mid\mathcal D,{\rm high}=A]=\kappa_{\rm ord}H\)，其中 \(H\)
是固定根对的任意函数。

对换 \(A,B\) 后，有序核变为 \((q-1)/a\)。在有序 \(B\)-law 下根交换把
\(\Theta\) 变为 \(1-\Theta\)，并保持 \(D,E,U\)。所以对于根交换不变的 \(H(D,E,U)\)，
\[
\int H\kappa_{\rm ord}\,d\lambda_n^{\rm ord}
=\int H\,\frac12\left(\frac{q-1}{a}+\frac{q-1}{b}\right)
d\lambda_n^{\rm ord}
=\int H\,m(a,b)\,d\lambda_n^{\rm ord}.
\tag{6}
\]
其中
\[
m(a,b)=\frac{q(q-1)}{2ab}.
\tag{7}
\]
故 r041 对 \(Q_n\)、对称距离事件等使用 \(m\) 的式子仍正确；限制仅在于不能把
该对称核用于依赖“哪一支是高根”的事件。

## 3. 高根支小的端点在完整双峰律下自动消失

固定一个未定向纤维，完整 \(K=2\) 交错下高根在大小为 \(a\) 的支的概率为 \(a/q\)。
若 \(0<\delta<1/2\)，至多一支可满足大小占比不超过 \(\delta\)。故逐纤维有
\[
\mathbb P_2(\Theta\le\delta\mid\mathcal D)\le\delta,
\]
并且无条件严格有
\[
\boxed{\quad\mathbb P_2(\Theta\le\delta)\le\delta.\quad}
\tag{8}
\]
又 \(0\le E\le2\)，所以
\[
\mathbb E_2[E{\bf1}_{\{\Theta\le\delta\}}]\le2\delta.
\tag{9}
\]
这不需要 \(A_2/B\) 的界或任何一致可积性。

因此，有序比例变量下 \(\Theta\downarrow0\) 的端点没有可能承载非零的 \(Q_n\) 极限。
与 r041 的对称比例不同，尚未控制的比例端点是
\(\Theta\uparrow1\)，即低根支 \(b/q\downarrow0\)。这不处理
\(\Theta\) 有界离开端点的平衡内部区；该区的空间几何同样仍未证。

还有一个只诊断高根支的精确恒等式：
\[
\mathbb E_2\frac q a=2.
\tag{10}
\]
确实，在任一未定向纤维中，高根在 \(A\) 的概率为 \(a/q\)，对应
\((a/q)(q/a)=1\)；高根在 \(B\) 的方向同样贡献 \(1\)。它与 (8) 一致，但不提供
低根小支或平衡内部区的空间控制。

## 4. 固定小高根支的定量控制

对整数 \(k\ge1\)，仍在每个未定向纤维中条件化。若高根在 \(A\) 且 \(a\le k\)，
则该方向的概率为 \(a/q\)，且合并前两支加上合并点构成含 \(q+1\) 点的连通子图，
故原图 \(\rho\le q\)。因此该方向
对 \(\mathbb E_2[\rho{\bf1}_{\{a\le k\}}]\) 的贡献至多 \(a\)。交换两支相加，得到
\[
\rho\left[
\frac aq{\bf1}_{\{a\le k\}}+\frac bq{\bf1}_{\{b\le k\}}
\right]
\le a{\bf1}_{\{a\le k\}}+b{\bf1}_{\{b\le k\}}\le2k.
\]
对所有纤维平均即有
\[
\boxed{\quad
\mathbb E_2[\rho{\bf1}_{\{\text{high-root branch size}\le k\}}]\le2k.
\quad}
\tag{11}
\]

所以对每个 \(\varepsilon>0\)，
\[
\mathbb P_2(\rho>\varepsilon n,\ \text{high branch}\le k)
\le\frac{2k}{\varepsilon n}.
\tag{12}
\]
并且写坐标差为 \((s,t)\)，有
\(s^2+t^2\le(n-1)(|s|+|t|)=(n-1)\rho\)，从而
\[
\mathbb E_2\!\left[E{\bf1}_{\{\text{high branch}\le k\}}\right]
\le\frac{2k(n-1)}{n^2}.
\tag{13}
\]
任取 \(k=o(n)\)，(12)--(13) 都趋于零。

## 5. 剩余问题

(8)--(13) 排除了高根支小造成的端点风险。它们不约束低根支小、即
\(b/q\to0\) 的词；singleton low branch 正属于这一端点，正是 r042 与 r045 所隔离的
部分。即使两支平衡，已有真实方格纤维仍可有宏观根距，而其总加权空间质量也未控制。
故本结果没有证明 \(Q_n\to0\)，也没有证明
\(\mathbb P_2(\beta>0)\to0\)。

## 6. \(3\times3\) 完整审计

audit_n3_ordered_kernel.py 完整枚举 \(9!\) 个排列，独立判定两个 birth、第一次合并、
top-two 前缀和低根首邻点。它精确验证：

- \(A_2=143808,\ B=81984\)；
- 每个非零 \((a,b)\) 桶有
  \(\overline L_B=(q-1)/b\)；
- \(A_{2,a,b}=B_{a,b}(q-1)/b\)；
- \(\sum_{B\text{ words}}\kappa_{\rm ord}=A_2\)；
- \(\sum_{K=2\text{ words}}q/a=2A_2\)，即 (10)；
- 对 \(\delta=1/4,1/3,2/5\) 验证 (8)，并对 \(1\le k\le7\) 验证 (11)。

这是有限整数/有理交叉核验；(3)--(13) 的一般性来自前述交错计数。
