# 第14轮：全局仿射方向强迫整个码的乘法稳定子平凡

## 1. 正确的命题

令 \(C\) 是 \(X=\mathbb F_2^{2k}\) 上的 full-bent 合法线性码，即 \(f_0+C\)
全为 bent。本文Walsh变换不归一化，记为W，幅值为q=2^k。设 \(M\subseteq C\) 是含常数的函数代数，并作**全模假设**
\[
 C\text{ 是 }M\text{-模}.
\tag{1}
\]
若 \(C\) 包含所有全局线性相位，则 \(M\) 只能有一个原子，故 \(M=\langle1\rangle\)。

这里 (1) 是必要量词：只知道容量定义中的某个 \(E\subseteq B\subsetneq C\) 且
\(B\) 为 \(E\)-模，并不足以推出结论。

## 2. 证明

令 \(A\) 为 \(M\) 的原子，则 \(1_A\in M\subseteq C\)。定义局部残余变换
\[
 Q=W(1_{A^c}(-1)^{f_0}),\qquad
 P_g=W(1_A(-1)^{f_0+g})\quad(g\in C_A).
\]
由 \(f_0+g\) 与 \(f_0+1_A+g\) 都 bent，
\[
 Q\pm P_g\in\{\pm q\};
\]
故逐频率
\[
 QP_g=0,\qquad Q^2+P_g^2=q^2.
\tag{2}
\]
因此 \(D_A=\{Q=0\}\) 与 \(g\in C_A\) 无关，所有 \(P_g\) 共同支撑于 \(D_A\)，
且在其上绝对值为 \(q\)。Plancherel 给
\[
 |D_A|=|A|.
\tag{3}
\]

令 \(\ell_\alpha(x)=\alpha\cdot x\) 为任意全局线性相位。由 (1)、
\(1_A\in M\)、\(\ell_\alpha\in C\)，有
\[
 1_A\ell_\alpha\in C_A.
\tag{4}
\]
在 \(A\) 上该函数正是 \(\ell_\alpha\)，故
\[
 P_{1_A\ell_\alpha}(\xi)=P_0(\xi+\alpha).
\]
局部共同支撑性给
\[
 D_A+\alpha=D_A\qquad(\forall\alpha\in X).
\tag{5}
\]
非空集合 \(D_A\) 在所有平移下不变只能为 \(X\)，故由 (3)
\(|A|=|X|\)。因而只有一个原子，\(M=\langle1\rangle\)。

## 3. 部分仿射方向

若只知道一个 \(u\)-维线性相位空间 \(U\subseteq C\)，同一证明给
\[
 D_A+U=D_A.
\]
所以 \(D_A\) 是 \(U\)-陪集的并，且
\[
 2^u\mid |D_A|=|A|.
\tag{6}
\]
若还使用 bent-difference 的最小单元下界 \(|A|\ge q\)，则原子数 \(b\) 满足
\[
 b\le\min\{q,\ |X|/2^u\}.
\tag{7}
\]

## 4. 为什么它不直接界定 \(\kappa(C)\)

容量允许 \(E\subseteq B\subsetneq C\)，只要求 \(B\) 为 \(E\)-模。此时
\(\ell_\alpha\in C\) 不代表 \(\ell_\alpha\in B\)，更不代表
\(1_A\ell_\alpha\in C\) 或 \(B\)。故 (4) 不能使用。

具体地，令
\[
 f_0(x,y)=x\cdot y,\qquad
 C=\{h(y)+a\cdot x:h:\mathbb F_2^k\to\mathbb F_2,\ a\in\mathbb F_2^k\}.
\tag{8}
\]
它包含全部全局仿射函数，且 \(f_0+C\) 全为 Maiorana--McFarland bent 函数。
取 \(E=B=\{h(y)\}\)。这是 \(C\) 的真子码中的函数代数，按 \(y\)-纤维分割为
\(q\) 个原子，给
\[
 \kappa(C)\ge2q.
\]
但 \(C\) 不是该 \(E\) 的全模：\(1_{\{y=y_0\}}x_i\notin C\)。相应局部谱
\(D_A\) 是一个频率纤维，而不受 \(x\)-频率平移不变。

因此“全模代数刚性”与“容量的部分模”必须明确区分：前者在 full affine closure 下
退化为常数，后者仍可有大容量。

本定理的全模假设不可省略，不是一般容量上界；本轮MM例子恰说明仿射闭合与大部分模块可以同时存在。
