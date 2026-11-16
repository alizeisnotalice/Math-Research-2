# r023：\(\tau=2\) 提升中 \(L\) 的首次合并尾部、二阶矩与真正缺口

设 \(\sigma\) 在 \(\tau=2\) 双峰词上均匀。固定首次合并状态、两侧内部单峰词和
合并后后缀，记两侧大小为 \(a,b\ge1\)，\(q=a+b\)。下述条件期望只对尚未固定的
AB/BA 方向与两侧交错取。

## 交错内的完整律

先取 AB 方向，即 A 根先出现、B 根为低峰。令 \(R=L-1\)。对于
\(1\le\ell\le a\)，有精确生存函数

\[
 \mathbb P(L\ge\ell\mid\mathrm{AB})
 ={\binom{a-1}{\ell-1}\over\binom{q-2}{\ell-1}}
 ={(a-1)_{\ell-1}\over(q-2)_{\ell-1}}, \tag{T1}
\]

而 \(\ell>a\) 时它为零。这里即使 \(b=1\) 也成立：首个 B 邻点是 merge
点 \(w\)，且 \(L=a=q-1\) 确定。对于 \(b\ge2\)，(T1) 只是在剩余交错的前
\(\ell-1\) 位全为 A 的概率。

以负超几何的 factorial moment（或由 (T1) 求和）得

\[
 \mathbb E[R\mid\mathrm{AB}]={a-1\over b},\qquad
 \mathbb E[R(R-1)\mid\mathrm{AB}]
 ={2(a-1)(a-2)\over b(b+1)}, \tag{T2}
\]

从而

\[
 \mu_{a,b}:=\mathbb E[L\mid\mathrm{AB}]={q-1\over b},
\]

\[
 s_{a,b}:=\mathbb E[L^2\mid\mathrm{AB}]
 =1+{3(a-1)\over b}+{2(a-1)(a-2)\over b(b+1)}. \tag{T3}
\]

BA 的所有式子只需交换 \(a,b\)。在不记录方向的 \(\tau=2\) 纤维中，两方向
等概率，故令

\[
 m(a,b)={\mu_{a,b}+\mu_{b,a}\over2}
 ={q(q-1)\over2ab}=\theta(a,b)^{-1}, \tag{T4}
\]

\[
 s(a,b)={s_{a,b}+s_{b,a}\over2}. \tag{T5}
\]

这给出完整的条件尾和二阶矩，而不只是平均的 \(1/\theta\)。

## 一个全 \(n\) 的非平凡二阶界

对 \(F=a-1\)，直接整理 (T3) 得

\[
 2\mu_{a,b}^2-s_{a,b}
 ={2F^2+Fb^2+3Fb+b^3+b^2\over b^2(b+1)}\ge0. \tag{T6}
\]

因此每一方向 \(s_{a,b}\le2\mu_{a,b}^2\)。平均方向后，利用
\(u^2+v^2\le(u+v)^2\)，有

\[
 s(a,b)\le4m(a,b)^2. \tag{T7}
\]

记 \(\mathcal D\) 为固定首次合并数据、两侧内部词和后缀（不记录方向及交错）。
于是对完整 \(\tau=2\) 律，

\[
 \mathbb E_{\tau=2}[L^2]\le4\,\mathbb E_{\tau=2}[m(\mathcal D)^2]. \tag{T8}
\]

这是所有自由方格、甚至上述任意首次合并纤维的严格界。它不是对
\(\mathbb E L^2/(\mathbb EL)^2\) 的统一界：不同纤维的 \(m=1/\theta\) 仍可相差
到 \(\Theta(N)\)。确有确定界

\[
 1\le m(a,b)\le q/2\le(N-1)/2, \tag{T9}
\]

因为 \(ab\ge q-1\)，但这在渐近上仍过弱。

### 常数 \(4\) 在纤维参数层面不能降低

这里的 \(4\) 不是两次粗略不等式造成的可随意损失。令
\(a=t\)、\(b=t^3\)。直接代入 (T3)--(T5) 得

\[
 m(t,t^3)=\frac{(t+t^3)(t+t^3-1)}{2t^4}
 =\frac{t^2}{2}+1+O(t^{-1}),
\]

而

\[
 s_{t,t^3}=O(1),\qquad
 s_{t^3,t}=2t^4+O(t^3),\qquad
 s(t,t^3)=t^4+O(t^3).
\]

从而

\[
 \frac{s(t,t^3)}{m(t,t^3)^2}\longrightarrow4. \tag{T9a}
\]

因此，对所有正整数参数 \(a,b\) 的纯纤维不等式
\(s(a,b)\le C m(a,b)^2\) 必有 \(C\ge4\)。这只说明 (T7) 的统一常数在
抽象的首次合并大小参数上最优；它**不**声称每个 \((t,t^3)\) 都能作为自由方格
中某一实际首次合并纤维的大小，也不产生原问题的空间矩下界。

## 几何量的精确相关性位置

令 \(H\) 是只依赖两根且对交换两根不变的量，例如
\(H=\|x-y\|_2^2/n^2\) 或远峰指示 \(\mathbf1\{\rho>\varepsilon n\}\)。在固定
\(\mathcal D\) 后，根和 \(H\) 已经确定；交错只随机化 \(L\)。故 r022 的
Rao--Blackwell 恒等式具体为

\[
 \mathbb E_{\tau=2}[LH\mid\mathcal D]=m(\mathcal D)H,
 \qquad
 \mathbb E_2H={\mathbb E_{\tau=2}[mH]\over\mathbb E_{\tau=2}m}.
 \tag{T10}
\]

等价地

\[
 \mathbb E_2H=\mathbb E_{\tau=2}H+
 {\operatorname{Cov}_{\tau=2}(m,H)\over\mathbb E_{\tau=2}m}. \tag{T11}
\]

所以缺少的不是交错内 \(L\) 的尾界或方差：那部分已由 (T1)--(T8) 完全控制。
真正待证的是首次合并几何 \(H\) 与不平衡权重 \(m=1/\theta\) 跨纤维的相关性。
此外首次合并前给出 \(\rho\le q\)，但这只是一侧上界，不能规定此相关性的符号。

## 可检验而尚未证明的充分条件

以下条件把后续工作精确化，并非本轮已经证明的估计。

1. **肯定方向。** 令 \(H_n=\|x-y\|_2^2/n^2\in[0,2]\)。若
   \(\mathbb E_{\tau=2}H_n\to0\)，且族
   \(m/\mathbb Em\) 在其自身 size-weight 下统一可积，即
   \[
   \lim_{K\to\infty}\limsup_n
   {\mathbb E[m\,\mathbf1\{m>K\mathbb Em\}]\over\mathbb Em}=0, \tag{T12}
   \]
   则 (T10) 蕴含 \(\mathbb E_2H_n\to0\)，故原题肯定。证明是把 (T10) 分为
   \(m\le K\mathbb Em\) 和补集：前项不超过 \(K\mathbb EH_n\)，后项不超过
   两倍 (T12) 的尾项。这是一个真正需要同时控制几何与熵/branch-balance 的条件。

2. **证伪方向。** 对 \(E_n=\{\rho>\varepsilon n\}\)，因 \(m\ge1\)，
   \[
   \mathbb P_2(E_n)\ge{\mathbb P_{\tau=2}(E_n)\over\mathbb E_{\tau=2}m}. \tag{T13}
   \]
   所以若沿无穷子列 \(\tau=2\) 律的远峰概率有固定正下界，且
   \(\mathbb Em\) 有界（等价于 top-two 扇区有固定正占比），便严格证伪。
   更接近“几乎全远”的充分条件为
   \[
   {\mathbb E[m^2]\over(\mathbb Em)^2}\,
   \mathbb P_{\tau=2}(E_n^c)\longrightarrow0,
   \tag{T14}
   \]
   它由 Cauchy--Schwarz 控制补集的 size-weight；(T8) 可用来将实际 \(L\)
   的二阶诊断换成 merge-weight 的二阶诊断。T13--T14 只是明确的可验证目标，
   不是用一般 Cauchy 不等式替代相关性研究。

目前 (T12) 的统一可积性、\(\tau=2\) 几何矩，或 (T13) 中的两个前提均未在全
\(n\) 证明。因此本轮没有改变原题的证明/证伪状态；它将第一处缺口缩小为
\(m\) 与几何的跨纤维联合律。
