# 第14轮：线性扩张的容量界与单个仿射方向的大跳

## 1. 任意线性扩张的通用容量界

设 \(C\subseteq C'\subseteq\mathbb F_2^X\) 都含 \(1_X\)，并写
\[
 d=\dim(C'/C).
\]
则
\[
 \boxed{\ \kappa(C')\le(d+1)\kappa(C).\ }
\tag{1}
\]

取实现 \(\kappa(C')\) 的 admissible 分割 \(\{A_j\}\)。每个单元指示函数
\(1_{A_j}\in C'\) 在商空间 \(C'/C\simeq\mathbb F_2^d\) 中有标签 \(\lambda_j\)，且
\(\sum_j\lambda_j=0\)，因为 \(\sum_j1_{A_j}=1_X\in C\)。

把这些标签分组为极小零和组。零标签为单点组；非零极小零和组的真子集线性无关，
故每组大小至多 \(d+1\)。每组的单元并 \(U\) 满足 \(1_U\in C\)，从而给出 \(C\)
的一个粗分割。对一组 \(J\)，支撑互不交给
\[
 \bigoplus_{A\in J}C'_A\subseteq C'_U,\qquad
 \dim C'_U\le\dim C_U+d.
\]
故该组的细分容量至多
\[
 \sum_{A\in J}(\dim C'_A+1)
\le\dim C_U+d+|J|
\le(\dim C_U+1)+2d.
\]
若粗分割有 \(p\) 组，其容量 \(K\ge2p\)，因为每个 \(1_U\in C_U\)。求和便得
\[
 \kappa(C')\le K+2dp\le(d+1)K\le(d+1)\kappa(C).
\]
特别地，一维扩张至多翻倍，但不可能有一般的加性或 codimension-Lipschitz 界。

## 2. 一维仿射 completion 的真实 bent 大跳例

令 \(Y=\mathbb F_2^k\), \(q=2^k\ge4\)，取 \(X=\mathbb F_2^k\times Y\)，变量写为
\((x,y)\)，并令
\[
 f_0(x,y)=x\cdot y,\qquad
 T=Y\setminus\{0,e_1\}.
\]
令
\[
 C=\{h:Y\to\mathbb F_2:\sum_{y\in T}h(y)=0\}.
\tag{2}
\]
把 \(h\) 视为 \(X\) 上只依赖 \(y\) 的函数。对任何 \(h\)，
\(x\cdot y+h(y)\) 是 Maiorana--McFarland bent 函数，故 \(f_0+C\) 完全合法；
\(1\in C\)，因为 \(|T|=q-2\) 为偶数。

在两个自由点 \(\{0,e_1\}\) 上 \(C\) 是全函数代数，容量为 4；在 \(T\) 上是长度
\(m=q-2\) 的偶重码。后者的每个 admissible 单元大小为偶数，且
\(\dim C_A=|A|-1\)，所以每格贡献 \(|A|\)，总容量为 \(m\)。因此
\[
 \kappa(C)=4+(q-2)=q+2.
\tag{3}
\]

取全局仿射方向 \(\ell(y)=y_1\)。有
\[
 \sum_{y\in T}\ell(y)=q/2-1\equiv1\pmod2,
\]
故 \(\ell\notin C\)。而 \(C\) 在全函数空间 \(\mathbb F_2^Y\) 中余维 1，所以
\[
 C+\langle\ell\rangle=\mathbb F_2^Y,\qquad
 \kappa(C+\langle\ell\rangle)=2q.
\tag{4}
\]
因此一个真正保 bent 的**单个全局仿射方向**可使
\[
 \kappa:\ q+2\longmapsto2q,
\]
绝对跳跃为 \(q-2\)。这严格否定任何只依赖扩张余维的加性或 Lipschitz 控制。
它仍服从 (1)，且 \(\kappa/q\) 只从 \(1+2/q\) 升至 2；不是 A28 反例。

## 3. 重复 completion 的范围

固定 \(X\) 上，仿射闭包
\[
 C\longmapsto C+RM(1,\dim X)
\]
幂等，不能无限重复。上述例子的跳跃来自 \(C\) 本身没有仿射闭包，而不是从已
affine-closed 初始码反复取得增长。

在 selector/tensor 构造中若不断引入新变量，(1) 只给每轮至多乘以“新增商维数加 1”。
是否存在真实 bent 例子使这种因子多轮接近取等、并使 \(\kappa/q\) 超过 2，本文没有
建立。故该例只否定对任意合法码使用“一个仿射方向只增加常数容量”的加性估计；它不排除相对于一个已仿射闭合初始码的更强联合预算。一般选择器与补全的组合是否可统一控制，仍需单独证明。

## 4. 无损控制的兼容子类与一般加回问题

若 \(E\subseteq B\subseteq C\)，且不仅 \(B\) 是 \(E\)-模，整个 \(C\) 也是 \(E\)-模，则
\[
S_Q(B,E)+C_{\rm diag}\subseteq S_Q(C,E).
\]
第13轮精确倍率与容量单调性立即给
\[
\kappa(S_Q(B,E)+C_{\rm diag})\le Q\kappa(C).
\]
这里 \(C_{\rm diag}\) 指与新变量无关的全部旧 \(C\) 方向；只加回全局仿射方向是其子空间。特别地，\(E=\langle1\rangle\) 总满足此结论。

一般 \(E\) 只稳定删减后的 \(B\)，不一定稳定 \(C\)，不能使用上述包含。对 \(Q=2\)，更强的纯函数代数问题可以写成
\[
T=C+wB+zE,\qquad
\kappa(T)\stackrel{?}{\le}2\kappa(C).
\tag{5}
\]
如果(5)对所有含单位函数码成立，递归便控制任意 \(Q=2^\ell\) 的旧方向加回；目前本轮只提供明确范围的有限核验，不把它当作定理。

一个自然的独立核估计已经被否定：令 \(\rho:T\to C\) 为 \(w=z=0\) 处评价，不能对每个新代数/模块对 \(F\subseteq D\subseteq T\) 声称
\[
\dim(F\cap\ker\rho)+\dim(D\cap\ker\rho)\le\dim B+\dim E.
\]
反例为 \(C=B=R\) 是维数 \(m>1\) 的函数代数、\(E=\langle1\rangle\)，取
\(F=D=R+wR\)。两个核均为 \(wR\)，左边 \(2m>m+1\)。但该例的总容量仍恰为 \(2\kappa(R)\)，因此它只否定分开估计投影与核的证明方法，没有否定(5)。

本轮原题状态：一般A28开放。容量大跳例不是超平方根无穷族；通用扩张因子 \(d+1\) 不能直接替代迭代中的统一预算。
