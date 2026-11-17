# r037：top-two 正规化下的首次合并比例—间隙混合极限

## 对象、记号与严格结论

令 \(G_n=P_n\square P_n\)、\(N=n^2\)，并在 top-two 双峰词集合 \(B_n\) 上取均匀律
\(\mathbb P_{B,n}\)。词的两个根是前两点；第一个根带全局最高标号，第二个根带次高标号。
令 \(a,b\ge1\) 是首次合并前分别从这两个根长出的两支大小，\(q=a+b\)，并令

\[
 \Theta_n:=\frac{a}{q}\in[0,1].
\tag{1}
\]

交换词的前两个根给出 \(B_n\) 上的一个无不动点对合：后续每个点出现时，已见集合
不变，且两个根原本不相邻。因此它交换 \(a,b\)，保留 \(q\)。所以
\(\Theta_n\) 的律 \(\mu_n\) 精确满足

\[
 \mu_n(E)=\mu_n(1-E).
\tag{2}
\]

等价地，也可在每个无序首次合并纤维上使用
\(\tfrac12(\delta_{a/q}+\delta_{b/q})\) 的对称化；两种写法对下述对称函数给出同一积分。

令 \(L\ge1\) 是第二根首次有邻点的间隙，\(h_t(n)\) 是 \(B_n\) 中 \(L=t\) 的计数。
本轮证明：

**定理。** 对任一趋于无穷的 \(n\) 子列，都可再取子列 \(n_k\)，使
\(\mu_{n_k}\Rightarrow\nu\) 于 \([0,1]\)，其中 \(\nu\) 是对称概率测度。对任意固定
\(t\ge2\)，

\[
 \boxed{\quad
 \frac{h_t(n_k)}{B(n_k)}
 \longrightarrow
 \int_0^1
 \frac{\theta(1-\theta)^{t-1}+(1-\theta)\theta^{t-1}}2\,d\nu(\theta).
 \quad}
\tag{3}
\]

并且

\[
 \frac{h_1(n)}{B(n)}\longrightarrow\frac12.
\tag{4}
\]

这是完整 \(B\)-law 的固定间隙结论，未附加远峰事件或非 singleton 假设。它不控制
\(\mathbb E_{B,n}L=A_2(n)/B(n)\)，也不推出 \(A_2/B\) 有界。

## 1. \(B\)-law 下首次合并秩趋于无穷

既有的 dyadic 回插定理给出

\[
 L_n:=\frac{B(n)}{A_1(n)}
 \ge\frac{[\lfloor\log_2n\rfloor-2]_+}{32}\longrightarrow\infty.
\tag{5}
\]

而父路径编码给出 \(B_q(n)\le3^{q-2}A_1(n)\)。故每个固定整数 \(K\ge2\) 都有

\[
 \mathbb P_{B,n}(q\le K)
 \le\frac{\sum_{q=2}^{K}3^{q-2}A_1(n)}{B(n)}
 =\frac{3^{K-1}-1}{2L_n}\longrightarrow0.
\tag{6}
\]

即 \(q\to\infty\) 依 \(B\)-law 概率收敛。这一步只利用 \(B/A_1\) 的发散，绝不使用
也不蕴含 \(A_2/B\) 的任何上界。

## 2. 每个纤维的精确律及固定 \(t\) 近似

先固定未记录根方向的首次合并数据和两支内部词，只平均哪一支根是第二根的两个方向。
若 \(a,b\ge2\)，置 \(x=a-1,y=b-1,m=q-2\)，则交错计数给出

\[
 p_t(a,b):=\mathbb P(L=t\mid a,b)
 =\frac{x(y)_{t-1}+y(x)_{t-1}}{2(m)_t},
 \qquad1\le t\le m.
\tag{7}
\]

超过此范围的 \(p_t\) 取零。若一支 singleton 且 \(q>2\)，则精确地

\[
 p_1(a,b)=p_{q-1}(a,b)=\frac12,\qquad
 p_t(a,b)=0\quad(t\notin\{1,q-1\});
\tag{8}
\]

\(a=b=1\) 时 \(p_1=1\)。这些式子正是有序根标号的条件律再作根交换平均；
(2) 说明它们可以在整个 \(B\)-law 下使用。

令

\[
 g_t(\theta):=
 \frac{\theta(1-\theta)^{t-1}+(1-\theta)\theta^{t-1}}2.
\tag{9}
\]

对 \(t\ge2\)，直接逐个比较 (7) 的有限乘积和相应幂，得到一条足够的统一界：

\[
 \sup_{a+b=q}\left|p_t(a,b)-g_t(a/q)\right|
 \le\frac{32t^2}{q},
 \qquad q\ge4t.
\tag{10}
\]

以下给出 (10) 的显式粗估计。令 \(k=t-1\)、\(M=m-1\)，并把 (7) 的第一项写作
\[
U=\frac{x}{m}\cdot\frac{(y)_k}{(M)_k}.
\]
若 \(y\ge k\)，每个因子满足
\[
\left|\frac{y-s}{M-s}-\frac yM\right|
=\frac{s(M-y)}{M(M-s)}\le\frac{s}{M-s}\qquad(0\le s<k).
\]
乘积望远镜相减，且各因子属于 \([0,1]\)，于是
\[
\left|\frac{(y)_k}{(M)_k}-\left(\frac yM\right)^k\right|
\le\frac{k^2}{M}.
\]
若 \(y<k\)，左侧第一个分式为零，而 \((y/M)^k\le k/M\le k^2/M\)，同一界仍成立。
又有
\[
\left|\frac xm-\frac{x+1}{m+2}\right|\le\frac1m,\qquad
\left|\frac yM-\frac{y+1}{m+2}\right|\le\frac6m
\]
（\(q\ge4t\)）。幂函数在 \([0,1]\) 上为 \(k\)-Lipschitz，故
\[
\left|U-\frac aq\left(1-\frac aq\right)^k\right|
\le\frac{6t^2}{m}\le\frac{12t^2}{q}.
\]
交换 \(x,y\) 并取二者平均后，(10) 随即成立，事实上上式已给出更小常数。
singleton 时固定 \(t\ge2\) 的 (8) 左侧为零，而
\(g_t(q^{-1})\le q^{-1}\)，故也被 (10) 支配。该常数不求最优。

对 \(t=1\)，(7)--(8) 给出 \(p_1=1/2\)（唯一例外 \(q=2\) 时为 \(1\)），故由 (6)
立刻得到 (4)。

## 3. 弱极限与混合公式

\([0,1]\) 紧，故任一 \(n\to\infty\) 子列均可再取 \(n_k\) 使
\(\mu_{n_k}\Rightarrow\nu\)。由 (2)，极限 \(\nu\) 仍对称。

对固定 \(t\ge2\)，(10)、(6) 蕴含

\[
 \frac{h_t(n)}{B(n)}
 =\mathbb E_{B,n}p_t(a,b)
 =\int_0^1g_t(\theta)\,d\mu_n(\theta)+o(1).
\tag{11}
\]

这里 \(o(1)\) 的证明不交换无限和：先以 (6) 丢弃 \(q<4t\)，其余部分由
\(32t^2/q\) 控制。由 \(q\to\infty\) 依概率收敛和 \(0\le1/q\le1\)，
\(\mathbb E_{B,n}(1/q)\to0\)。函数 \(g_t\) 连续且有界，
故弱收敛将 (11) 提升为 (3)。

## 4. 有限间隙逃逸质量与 tightness 的精确判据

为统一记号，令 \(g_1\equiv1/2\)。对 \(0<\theta<1\)，几何级数给出

\[
 \sum_{t\ge1}g_t(\theta)=1.
\tag{12}
\]

而在端点 \(\theta=0,1\)，有 \(g_1=1/2\) 且 \(g_t=0\)（\(t\ge2\)），所以和为
\(1/2\)。由单调收敛，

\[
 \sum_{t\ge1}\lim_{k\to\infty}\frac{h_t(n_k)}{B(n_k)}
 =1-\frac{\nu(\{0\})+\nu(\{1\})}{2}.
\tag{13}
\]

因而沿该子列，有限 \(L\) 桶中逃向 \(L=\infty\) 的质量严格为

\[
 \boxed{\quad
 \lim_{K\to\infty}\lim_{k\to\infty}
 \mathbb P_{B,n_k}(L>K)
 =\frac{\nu(\{0\})+\nu(\{1\})}{2}.
 \quad}
\tag{14}
\]

因此，\(L\) 在 \(\mathbb P_{B,n}\) 下紧，当且仅当每个 \(\mu_n\) 的子序列弱极限
\(\nu\) 都满足 \(\nu(\{0,1\})=0\)。必要性由 (14) 直接给出。反之，若 \(L\) 不紧，
可取 \(n_k\) 和 \(K_k\uparrow\infty\)，使 \(\mathbb P(L>K_k)\ge\varepsilon\)；
再令 \(\mu_{n_k}\Rightarrow\nu\)。若端点质量为零，则 (13) 的有限部分质量为 \(1\)，
选固定 \(K\) 使其极限累计质量 \(>1-\varepsilon/2\)，再用
\(\{L>K_k\}\subseteq\{L>K\}\) 矛盾。

这是概率 tightness，而非 \(L\) 的一致可积或有界均值。即使 (14) 的端点缺失为零，
也仍可能有很小概率上的很大 \(L\) 使 \(\mathbb E_{B,n}L=A_2/B\) 发散；因此本定理
既不证明也不暗示 \(A_2/B\) 有界。

## 复核与范围

解析输入为：

- a022-reinsertion-bootstrap-r036/DYADIC_BOX_PROOF.md 的 (D1)、(D12)；
- a022-small-merge-normalizer-r033/EXPONENTIAL_PATH_ENCODING.md 的 \(B_q\) 编码；
- a022-archive-final-r023/tau2-lift-theory-r023/tau2-lift-theory.md 的固定纤维交错律。

audit_fixed_gap_mixture.py 用标准库 fractions.Fraction 穷尽有限 \(2\le q<80\) 的
纤维，交叉核验 (7)--(8) 的归一化、singleton 边界及 (10) 的非锐常数。在此范围内的
通过只是索引/常数的有限精确核验；(3)--(14) 的全 \(n\) 结论由上面的解析论证承担。
