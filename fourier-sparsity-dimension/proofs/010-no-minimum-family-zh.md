# 非 weakly-normal bent 种子的直和放大及其真实谱图提升

设 \(f_*:X=\mathbb F_2^{2k_0}\to\mathbb F_2\) 为 bent，且不在任何
\(k_0\)-维仿射平坦上仿射（non-weakly-normal）。令
\(Y=Z\oplus W\cong\mathbb F_2^{2\ell}\)，\(Q(z,w)=z\cdot w\)，并令
\[
 F_\ell(x,z,w)=f_*(x)+Q(z,w),\qquad q_\ell=2^{k_0+\ell}. \tag{1}
\]

## 来源与核验层级

本文的明确种子为
\[
 f_*(x)=\operatorname{Tr}_{\mathbb F_{2^{14}}/\mathbb F_2}(\alpha x^{57}),
 \qquad\alpha\in\mathbb F_4\setminus\mathbb F_2.
\]
其非弱正规性引用 Canteaut–Daum–Dobbertin–Leander，
*Finding nonnormal bent functions*，Discrete Applied Mathematics 154 (2006),
202–218，Fact 22（第209页），
[DOI](https://doi.org/10.1016/j.dam.2005.03.027)。该事实由原论文的穷举算法
认证，作者明确没有给出这一基例的解析证明。本轮主Agent实际核读该页、
Definition 6 和 Lemma 25；重新构造有限域并用完整整数Walsh变换核验
种子的 bent 性，**没有重跑其十四变量非弱正规性穷举**。

因此下文的一般放大、维数及谱图结论有自包含证明；采用上述明确种子时，
其“没有最小翻号”性质保留对文献计算认证的依赖。不是本轮新提供了该
基例的解析非弱正规证明。二次直和保持弱正规性也见该文 Lemma 25；下面
给出统一的辛商空间证明。

## 1. non-weak-normal 性在直和下保持

反设 \(F_\ell\) 在一个 \((k_0+\ell)\)-平坦 \(A\) 上仿射。令 \(L\) 是
其方向，\(K=L\cap Y\)、\(r=\dim K\)、\(U=\operatorname{pr}_XL\)。则
\[
 \dim U=k_0+\ell-r. \tag{2}
\]
在纯 \(K\)-方向取二阶差分可知 \(K\) 对
\[
 B_Q((z,w),(z',w'))=z\cdot w'+z'\cdot w \tag{3}
\]
完全各向同性，故 \(r\le\ell\)。选择 \(L\to U\) 的线性截面，得
\[
 A=\{(x_0+u,y_0+Tu+k):u\in U,k\in K\} \tag{4}
\]
（平移不影响以下二阶差分）。混合 \((u,k)\) 差分给
\(B_Q(Tu,k)=0\)，故 \(T(U)\subseteq K^\perp\)。令 \(k=0\) 并比较
两个 \(u\)-方向，得到
\[
 f_*|_{x_0+U}+Q\circ T\quad\text{仿射}. \tag{5}
\]
\(Q\circ T\) 的极化型经过辛商 \(K^\perp/K\)，所以其秩至多
\(2(\ell-r)\)。秩 \(2s\) 的交错型在 \(b\)-维空间中有维数 \(b-s\)
的完全各向同性子空间；于是在 \(U\) 内存在维数至少
\[
 (k_0+\ell-r)-(\ell-r)=k_0 \tag{6}
\]
的子空间，使 \(Q\circ T\) 限制为仿射。由 (5)，\(f_*\) 在一个
\(k_0\)-平坦上仿射，矛盾。

因此 \(F_\ell\) 没有重量 \(q_\ell\) 的合法翻号：这种翻号的标准平坦
刻画会给出一个 \((k_0+\ell)\)-平坦上的仿射限制。该结论依赖种子本身，
故对任意包含它的合法超空间仍成立。

## 2. 大合法空间

令 \(V_*\) 是 \(f_*\) 的任意线性合法方向空间。令 \(\mathcal T_\ell\)
由三角位移组成：
\[
 a_j=a_j(w_0,\ldots,w_{j-1})\quad(0\le j<\ell).
\]
则 \(w\mapsto w+a(w)\) 是置换且
\(\dim\mathcal T_\ell=\sum_j2^j=2^\ell-1\)。定义
\[
 \mathcal V_\ell=
 \{h_w(x)+z\cdot a(w):(h_w)_{w\in W}\in V_*^W,
 a\in\mathcal T_\ell\}. \tag{7}
\]
对固定 \(w\)，相位为
\[
 f_*(x)+h_w(x)+z\cdot(w+a(w)). \tag{8}
\]
先对 \(z\) 求 Walsh 和，立即得到所有 \(F_\ell+v\) bent。又两项对
\(z\) 的依赖不同，故
\[
 \dim\mathcal V_\ell=2^\ell\dim V_*+2^\ell-1. \tag{9}
\]

取 \(V_*=C_X\)（全体仿射函数，维数 \(2k_0+1\)），则
\[
 D_\ell:=\dim\mathcal V_\ell=(2k_0+2)2^\ell-1. \tag{10}
\]
文献给出的 CDDL 十四变量种子可取 \(k_0=7\)、\(q_*=128\)、
\(\dim C_X=15\)，从而 \(D_\ell=16\,2^\ell-1=q_\ell/8-1\)。本笔记
不重证该种子的 non-weak-normal 性；上节仅证明一旦有该性质，直和保持它。

## 3. 回到原 A28 的合法谱图提升

取 \(\mathcal V_\ell\) 的一组基，写评价向量为
\(v(\xi)\in\mathbb F_2^{D_\ell}\)（\(\xi\in X\times Y\)）。定义
\[
 H(a,b)=q_\ell^{-1}\sum_{\xi\in X\times Y}
 (-1)^{F_\ell(\xi)+a\cdot\xi+b\cdot v(\xi)}. \tag{11}
\]
对每个 \(b\)，\(b\cdot v\in\mathcal V_\ell\)，故
\(F_\ell+b\cdot v\) bent；其未归一化 Walsh 值模为 \(q_\ell\)。所以
\(H\) 是 Boolean。它的 Fourier 支撑正是图
\[
 \Gamma=\{(\xi,v(\xi)):\xi\in X\times Y\},\qquad
 s=|\Gamma|=q_\ell^2. \tag{12}
\]

此处线性秩与仿射维数须分开：\(\mathcal V_\ell\) 含全部 \(\xi\) 的
仿射坐标及常数函数（\(h_w\) 可逐 \(w\) 选仿射常数/\(x\)-坐标，
三角位移含各 \(z_j\)）。故图坐标函数与常数张成的函数空间恰为
\(\mathcal V_\ell\)。又常数函数本身已是 \(v\)-坐标的线性组合，因而
\[
 \dim\operatorname{span}\Gamma=D_\ell,\qquad
 \operatorname{affdim}\Gamma=D_\ell-1. \tag{13}
\]
所以这是真正的 A28 Fourier 对象，而非只在 bent 码中的计数：
\[
 \frac{s}{d^2}=\frac{q_\ell^2}{(q_\ell/8-1)^2}\longrightarrow64. \tag{14}
\]
它不构成 \(s=o(d^2)\) 反例，却提供无穷个具有真正 \(d=\Theta(\sqrt s)\)
谱维、同时任何合法超空间都没有最小翻号的实例。因此一般“近满最小翻号
覆盖”不可能作为无条件归约。它不否定另一种带条件的潜在归约，例如只对
假设维数已超线性增长的反例候选要求近满覆盖；这种归约本轮并未证明。
