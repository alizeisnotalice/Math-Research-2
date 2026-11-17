# 第13轮：任意模块代数下的选择器容量精确倍率

## 1. 容量及其分割表达

设 \(C\) 是有限非空集合 \(X\) 上的二元函数线性空间，且 \(1_X\in C\)。乘法均为逐点乘法。定义
\[
\kappa(C)=\max_{1\in E\subseteq B\subseteq C,\ E\text{ 为代数},\ EB\subseteq B}
(\dim B+\dim E).
\tag{1}
\]
空间有限，最大值存在。这里不要求 \(C\) 自身是代数。显然 \(\dim C+1\le\kappa(C)\)，且包含关系 \(C_1\subseteq C_2\) 给 \(\kappa(C_1)\le\kappa(C_2)\)。

有限二元函数代数 \(E\) 等于某个分割 \(\mathcal P\) 的全部分块常值函数。若 \(C_A=\{c\in C:\operatorname{supp}c\subseteq A\}\)，则固定 \(E\) 时最大可行模块为 \(\bigoplus_{A\in\mathcal P}C_A\)：任意 \(b\in B\) 可按块指标分解，反向这些局部子码的直和确为模块。因此
\[
\boxed{\kappa(C)=\max_{\mathcal P:\,1_A\in C\ (A\in\mathcal P)}
\sum_{A\in\mathcal P}(\dim C_A+1).}
\tag{2}
\]
每个分割块非空；常数与各块的单位均按其所在集合解释。

## 2. 对稳定代数可无损分块

记 \(M(C)=\{a:aC\subseteq C\}\)。因 \(1\in C\)，它是包含于 \(C\) 的含单位代数。设 \(A\subseteq M(C)\) 为含单位子代数。给任意(1)的可行对 \((E,B)\)，把 \(AE,AB\) 定义为全部乘积的线性张成。则
\[
AE\subseteq AB\subseteq C,\qquad (AE)(AB)\subseteq AB.
\]
两空间均包含原空间，且 \(AE\) 是包含 \(A\) 的代数。因此求最大值时可要求 \(A\subseteq E\)。若 \(P_i\) 为 \(A\) 的原子，\(C_i=1_{P_i}C\) 视为 \(P_i\) 上的函数码，则可行对逐原子直和，得到
\[
\kappa(C)=\sum_i\kappa(C_i).
\tag{3}
\]
特别地，对不交定义域上的函数码，容量在直和下相加。函数码在定义域上作双射换标不改变容量；重复每个点相同次数也不改变容量，因为乘法与包含关系的函数代数同构保留所有可行对。

## 3. 全参数三角选择器的精确倍率

设 \(Q=2^\ell\)，\(W=Z=\mathbb F_2^\ell\)。定义
\[
S_Q(C,A)=\left\{h_w(x)+\sum_{j=0}^{\ell-1}z_j a_{j,w_{<j}}(x):
h_w\in C,\ a_{j,p}\in A\right\}.
\tag{4}
\]
所有列出的参数独立。[第11轮构造](011c-module-selector-budget-zh.md)证明：若 \(f_0+C\) 全 bent，则基准相位 \(f_0(x)+z\cdot w\) 加(4)的任意成员仍 bent，且
\[
\dim S_Q(C,A)=Q\dim C+(Q-1)\dim A.
\tag{5}
\]
下面的容量结论是纯函数代数恒等式，本身不需要 bent 假设：
\[
\boxed{\kappa(S_Q(C,A))=Q\kappa(C).}
\tag{6}
\]

**证明。** 旧代数 \(A\) 仍稳定(4)，按(3)可在每个旧原子独立工作。这把问题化为 \(A=\langle1\rangle\)；它没有要求候选新代数原本尊重旧原子，因为第2节可先扩大代数与模块，目标值不减。

按第一对变量 \(z=z_0,w=w_0\) 分解，令 \(K=S_{Q/2}(C,\langle1\rangle)\)，则
\[
L(K)=S_Q(C,\langle1\rangle)
=\{h_w(y)+az:h_0,h_1\in K,\ a\in\mathbb F_2\}.
\]
其无 \(z\) 部分为 \(H=K\oplus K\)，故 \(\kappa(H)=2\kappa(K)\)。对任意 \(E\subseteq B\subseteq L(K)\) 的可行对，分三种情况。

若 \(E\) 含 \(e=h+z\)，对每个 \(e'\in E\cap H\)，乘积 \(ee'\) 的 \(z\) 系数是 \(e'\)，必须为常数。故 \(E\cap H=\langle1\rangle\)，并且 \(E=\langle1,e\rangle\)。任取 \(b=h_b+a_bz\in B\)，则
\[
eb=hh_b+z(h_b+a_bh+a_b).
\]
括号必须为常数，所以 \(b\in\langle1,e\rangle\)。因此 \(B=E\)，得分为 \(4\le2\kappa(K)\)。

若 \(E\subseteq H\) 但 \(B\) 含 \(h+z\)，则对每个 \(e\in E\)，乘积 \(e(h+z)\) 的 \(z\) 系数为 \(e\)，故 \(E=\langle1\rangle\)。此时得分至多
\[
\dim L(K)+1=2\dim K+2\le2\kappa(K).
\]

若 \(B\subseteq H\)，得分不超过 \(\kappa(H)=2\kappa(K)\)。三种情况穷尽全部可行对，反向下界由 \(H\subseteq L(K)\) 给出，故 \(\kappa(L(K))=2\kappa(K)\)。从 \(Q=1\) 归纳并按旧原子求和，即得(6)。\(\square\)

这包含依赖全部新 \(z,w\) 变量的代数，以及任意混合旧原子的代数；第12轮仅允许旧原子内的 \(w\)-局部分割，现在不再需要这一限制。

## 4. 整个构造闭包的最优渐近密度

从一个 bent 合法码 \(C_0\) 出发，\(q_0=\sqrt{|X|}\)。每步允许任选 \(1\in E_i\subseteq B_i\subseteq C_i\)，其中 \(E_i\) 为代数、\(B_i\) 是其模，再取
\[
C_{i+1}=S_{Q_i}(B_i,E_i),\qquad q_{i+1}=Q_iq_i.
\]
允许任意有限次、任意可变尺度及任意子码删减。由单调性和(6)，
\[
\frac{\kappa(C_{i+1})}{q_{i+1}}
=\frac{\kappa(B_i)}{q_i}\le\frac{\kappa(C_i)}{q_i}
\le\frac{\kappa(C_0)}{q_0}.
\tag{7}
\]
特别地 \(\dim C_i/q_i<\kappa(C_0)/q_0\)。定义域双射换标也不改变此结论。因此从固定种子出发，这整个构造机制不能产生无界维数密度。

这个界还是该机制的精确渐近上确界。取实现 \(\kappa(C_0)=\dim B+e\) 的一对 \(B,E\)，其中 \(e=\dim E\)，单次选择器即有
\[
\frac{\dim S_Q(B,E)}{q_0Q}
=\frac{\kappa(C_0)}{q_0}-\frac{e}{q_0Q}
\longrightarrow\frac{\kappa(C_0)}{q_0}.
\tag{8}
\]
无需依靠多轮来提高这个极限。若某个种子满足 \(\kappa(C_0)>2q_0\)，取二次幂 \(Q\) 使 \(Q(\kappa(C_0)-2q_0)>e-1\)，便得到实际全 bent 码维数超过 \(2q_0Q-1\)。这仍只是常数密度增益，不是A28反例。

## 5. 范围与尚未解决的问题

(6)—(8)控制任意删减、改选模块代数及全参数三角选择器的组合，不限制新代数的变量依赖。它不涵盖向输出码再添加原来不在其中的合法方向、使用非三角置换族或其他同步拼接。尤其没有证明任意合法码均可由该机制从统一有界容量种子生成。

一般 \(\kappa(C)=O(q)\) 仍未证明；仅由 \(\dim C+1\le\kappa(C)\)，这样的界已经包含一般全 bent 合法码的目标秩界。由固定种子闭包有界不能推断随规模变化的种子容量统一有界，也不能据此宣布一般非平坦A28完成。第11轮提供构造合法性，本轮独立审查承担容量证明，有限计算另见第13轮记录。不声称首次性。
