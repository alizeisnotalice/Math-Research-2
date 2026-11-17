# Top-two normalizer 的根距一致性约束

在有限连通简单图上采用均匀双射标号及原图距离。令 \(B\) 为“两个最高标签均为峰”的双峰词数，\(A=A_2\)，\(\bar m=A/B\)。首次合并状态 \(F\) 的两个合并前支大小为 \(a,b\)，记 \(q=a+b\)，并令
\[
m(F)=\frac{q(q-1)}{2ab},\qquad\theta(F)=1/m(F).
\]
既有 top-two 提权恒等式给出
\[
\frac BA=\mathbb E_2\theta(F),\qquad \bar m=A/B.\tag{1}
\]

两个根在两个连通支中。经支内路径与桥点相连，根距 \(d\le q\)。由于 \(ab\le q^2/4\)，
\[
m\ge2(1-1/q)\ge2(1-1/d),\qquad
\theta\le\frac{d}{2(d-1)}.\tag{2}
\]
这些都是同一状态上的点态关系，不需要 \(m\) 与根距独立。

## 单阈值形式

若整数 \(r\ge2\)，则 \(d>r\) 蕴含 \(d\ge r+1\)。在 \(d>r\) 上 \(m\ge2r/(r+1)\)，而全域 \(m\ge1\)。故
\[
\frac BA\le1-\frac{r-1}{2r}\mathbb P_2(d>r),
\]
即
\[
\mathbb P_2(d>r)\le\frac{2r}{r-1}(1-B/A),
\qquad
\bar m\ge\frac1{1-\frac{r-1}{2r}\mathbb P_2(d>r)}.\tag{3}
\]

## 全距离分布形式

直接将 (2) 对全双峰律平均，得到较强的线性必要条件
\[
\boxed{\quad\frac BA\le\sum_{d\ge2}\frac{d}{2(d-1)}\mathbb P_2(\rho=d).\quad}\tag{4}
\]
等价地
\[
1-\frac BA-\sum_{d\ge2}\frac{d-2}{2(d-1)}\mathbb P_2(\rho=d)\ge0.\tag{5}
\]
有限 SMC 若用同一完整 \(K=2\) 的 \(L\)-提权距离直方图与 \(\bar m\) 代入而违反 (3)/(5)，只能说明这组有限比值不能同时近似该精确结构，或字段/权重尺度有错配；不单独证明程序错误、显著性或混合失败。

## 对根对对称函数的 Rao--Blackwell 候选

在 \(B\) 律下将两个最高峰以标签次序记为 \(u_1,u_2\)。记 \(j_i\) 为 \(u_i\) 的首个邻点出现的秩，并令原提升因子为 \(L=j_2-2\)。交换降序秩 1、2（即最高标签 \(N,N-1\)）是 \(B\) 上的保测度对合；它把 \(L\) 换成 \(j_1-2\)，并保持任何只依赖无序根对的量 \(H(\{u_1,u_2\})\)。所以
\[
 \mathbb E_B[LH]=\mathbb E_B[L_{\rm sym}H],\qquad
 L_{\rm sym}=\frac{j_1+j_2-4}{2}.\tag{6}
\]
若 \(d(u_1,u_2)>2\)，两个峰没有共同邻点。秩 3 至多是一个峰的首邻点，故 \(j_1+j_2\ge7\)，从而 \(L_{\rm sym}\ge3/2\)。这给出一个对根对几何量的低成本 Rao--Blackwell 化候选；它不保证有限 SMC 的总方差下降，也不改变 normalizer 的跨重复集中问题。


## 汇总对的确定性误差下界

令 \(\widehat b\) 是报告的 \(B/A\)，\(\widehat p_d\) 是报告的距离分布（将保存的有限小数作为有理数后精确归一化），且
\[
\delta=\widehat b-\sum_{d=2}^{D}g_d\widehat p_d>0,
\qquad g_d=\frac d{2(d-1)},\quad D=2n-2.
\]
对真实的 \(b=B/A,p_d=\mathbb P_2(\rho=d)\)，由 (4) 及
\(\operatorname{TV}(p,\widehat p)=\tfrac12\sum_d|p_d-\widehat p_d|\)，有
\[
\delta\le |\widehat b-b|+(1-g_D)\operatorname{TV}(p,\widehat p).
\]
于是
\[
\max\{|\widehat b-b|,\operatorname{TV}(p,\widehat p)\}
\ge\frac{\delta}{2-g_D}.
\]
这是针对已报告数值的确定性误差下界，不依赖置信水平。脚本以 `Fraction` 将两份 n=24 汇总的小数视为有理数核验：本机 strat 的右端严格大于 0.15，MacBook multi 的右端严格大于 0.18。因此每份汇总的 normalizer 比值绝对误差与距离分布总变差误差，至少一项超过相应阈值。它不指定哪一项有误，不归因于具体程序缺陷，也不外推原题极限。
