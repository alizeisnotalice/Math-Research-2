# 第15轮：选择器后加回旧方向的容量预算

本文沿用[第13轮](013-selector-capacity-zh.md)的容量定义。所有函数空间均为有限非空集合上的二元线性空间，乘法逐点进行。变量之间独立，维数取自 \(\mathbb F_2\)。首先证明纯函数代数结论，再说明 bent 合法性及迭代范围。

## 1. 主定理

设
\[
1\in E\subseteq B\subseteq C,\qquad EE\subseteq E,\quad EB\subseteq B.
\tag{1}
\]
不要求 \(E\) 稳定整个 \(C\)。对独立布尔变量 \(w,z\)，令
\[
T=C+wB+zE.
\]
则
\[
\boxed{\kappa(T)\le\kappa(C)+\kappa(B)\le2\kappa(C).}
\tag{2}
\]
特别地，第14轮的一般加回候选现在得到证明。第14轮指出的独立投影核估计仍为假；本次把两种新变量的核一起计入 \(\kappa(B)\)，没有重新使用该错误估计。

## 2. 允许零代数的辅助容量

对任意函数线性空间 \(N\)，不要求 \(1\in N\)，定义
\[
\sigma(N)=\max\{\dim F+\dim D:
 F\subseteq D\subseteq N,\quad FF\subseteq F,\quad FD\subseteq D\}.
\tag{3}
\]
允许 \(F=0\)。这里也不要求 \(F\) 的内部单位作用在 \(D\) 上为恒等；全部要求只有式(3)。

若 \(1\in N\)，则 \(\sigma(N)=\kappa(N)\)：任意(3)的对可扩大为
\(F+\langle1\rangle\subseteq D+\langle1\rangle\)，乘法闭与模块条件仍成立，维数不减。反向包含显然。

若含单位代数 \(A\) 稳定 \(N\)，任意候选对可扩大成 \(AF,AD\)，其中乘积表示线性张成。因为 \(1\in A\)，扩张包含原对；因为 \(AA\subseteq A\)、\(FF\subseteq F\)、\(FD\subseteq D\)，扩张仍可行。故可在最大化时要求两空间都是 \(A\)-模，继而按 \(A\) 原子直和分解。记 \(N_A=1_A N\)，则
\[
\sigma(N)=\sum_A\sigma(N_A).
\tag{4}
\]
反向由各块候选对的直和给出。当 \(N\) 含单位时，式(4)也恢复第13轮的 \(\kappa\) 直和性质。

## 3. 两个核的精确联合预算

关键恒等式为
\[
\boxed{\sigma(wB+zE)=\kappa(B).}
\tag{5}
\]

对角 \(E\) 稳定 \(N=wB+zE\)，因为 \(EB\subseteq B\)。在一个 \(E\) 原子上，识别局部单位为1，令 \(K=1_A B\)，问题变为
\[
N_A=wK+z\mathbb F_2,\qquad 1\in K.
\tag{6}
\]
注意这里没有常数项：所有成员都在 \(w=z=0\) 处为零。对(6)中的任意可行对 \(F\subseteq D\)，分三种情况。

**第一种：\(F\) 含 \(f=wb+z\)。** 若 \(wa\in F\)，则
\[
f(wa)=wab+wza\in N_A
\]
强迫 \(a=0\)，因为 \(N_A\) 不含 \(wz\) 项。因此 \(F\cap wK=0\)，由 \(z\) 系数至多一维知 \(F=\langle f\rangle\)。再取任意 \(d=wb'+\alpha z\in D\)，模块条件给
\[
fd=wbb'+\alpha z+wz(b'+\alpha b)\in N_A.
\]
所以 \(b'=\alpha b\)，即 \(d=\alpha f\)。于是 \(D=F\)，得分为2，不超过 \(\kappa(K)\)。

**第二种：\(F\subseteq wK\)，但 \(D\) 含 \(wb+z\)。** 对任意 \(wa\in F\)，同一 \(wz\) 系数强迫 \(a=0\)。故 \(F=0\)，得分
\[
\dim D\le\dim K+1\le\kappa(K).
\]

**第三种：\(D\subseteq wK\)。** 去掉公共因子 \(w\)，因 \(w^2=w\) 得到 \(K\) 中的非含单位代数/模对。第2节的补单位操作给其得分不超过 \(\kappa(K)\)。

这三种情况穷尽全部候选。由(4)和 \(E\) 稳定整个 \(B\)，
\[
\sigma(N)\le\sum_A\kappa(1_A B)=\kappa(B).
\]
反向取 \(B\) 中实现 \(\kappa(B)\) 的含单位对 \(F_0\subseteq D_0\)，则
\(wF_0\subseteq wD_0\subseteq N\) 是(3)的可行对，得分恰为 \(\kappa(B)\)。式(5)成立。

## 4. 评价的像与核

取 \(T\) 中任意实现容量的含单位对 \(F\subseteq D\)。评价
\[
\rho:T\longrightarrow C,\qquad \rho(t)=t|_{w=z=0}
\]
是环境函数代数的同态限制。像对 \(\rho(F)\subseteq\rho(D)\subseteq C\) 保持含单位、代数与模块条件，因此得分至多 \(\kappa(C)\)。

核对
\[
F\cap\ker\rho\subseteq D\cap\ker\rho\subseteq wB+zE
\]
是式(3)的可行对，不必含单位。其得分由(5)至多 \(\kappa(B)\)。两次秩—核维数公式相加，给
\[
\dim F+\dim D\le\kappa(C)+\kappa(B).
\]
再用单调性 \(\kappa(B)\le\kappa(C)\)，即得主定理。

该上界不是普遍等式。本轮有限严格例给出 \((\kappa(C),\kappa(B),\kappa(T))=(6,4,9)\) 以及 \((8,6,12)\)；详见计算证书。\(B=C\) 时则由第13轮取等，故系数2不能缩小。

## 5. 任意选择器尺度

令 \(Q=2^\ell\)。沿用全参数三角选择器
\[
S_Q(B,E)=\left\{h_w(x)+\sum_{j=0}^{\ell-1}z_j a_{j,w_{<j}}(x):
 h_w\in B,\ a_{j,p}\in E\right\},
\]
并令
\[
T_Q=S_Q(B,E)+C_{\mathrm{diag}}.
\]
则有更精细的预算
\[
\boxed{Q\kappa(B)\le\kappa(T_Q)
\le\kappa(C)+(Q-1)\kappa(B)\le Q\kappa(C).}
\tag{7}
\]
下界来自 \(S_Q(B,E)\subseteq T_Q\) 与第13轮精确倍率。

上界对 \(Q\) 归纳。\(Q=1\) 时 \(T_1=C\)。按首对 \(w_0,z_0\) 分解，写
\[
K=S_{Q/2}(B,E),\qquad T_{Q/2}=K+C_{\mathrm{diag}}.
\]
两条 \(w_0\) 分支的差在 \(K\) 中，故
\[
T_Q=T_{Q/2}+w_0K+z_0E_{\mathrm{diag}}.
\]
这里 \(E_{\mathrm{diag}}\subseteq K\subseteq T_{Q/2}\)，且对角 \(E\) 稳定整个 \(K\)。式(2)及第13轮给
\[
\kappa(T_Q)\le\kappa(T_{Q/2})+\kappa(K)
\le\kappa(C)+(Q/2-1)\kappa(B)+(Q/2)\kappa(B).
\]
于是(7)成立。

## 6. bent 合法性、仿射恢复与闭包范围

若 \(f_0+C\) 是 \(X=\mathbb F_2^{2k}\) 上的共同全bent族，则对每个 \(c\in C\)，
\(f_0+c+B\) 也是共同全bent族。对它使用同一选择器合法性定理，即知基准相位
\(f_0(x)+z\cdot w\) 加 \(T_Q\) 中任一成员仍 bent。这一步没有要求 \(E\) 稳定整个 \(C\)。

令 \(q=\sqrt{|X|}\)。每轮允许从当前 \(C_i\) 中任意删减为含单位的模块 \(B_i\)、任选满足(1)的代数 \(E_i\)、作任意尺度选择器，并加回当前 \(C_i\) 的全部或任意部分旧方向，再作任意仍含单位的删减。由(7)，
\[
\frac{\kappa(C_{i+1})}{q_{i+1}}
\le\frac{\kappa(C_i)}{q_i}.
\tag{8}
\]
故从固定种子出发，这个比第13轮更大的构造闭包仍不能产生无界维数密度。

若当前 \(C_i\) 包含全部全局仿射相位，旧仿射方向在 \((C_i)_{\rm diag}\) 内；全部新增 \(w_j,z_j\) 线性方向已经在选择器中。因此“删减—选择器—恢复所有全局仿射方向”也受(8)控制。若初始码尚未仿射闭合，可先补齐一次；第14轮的大跳可能发生，但固定种子的这一次有限扩张不能在上述机制内跨轮累积为无界密度。

整个闭包的最佳渐近密度仍恰为初始 \(\kappa(C_0)/q_0\)：上界为(8)，第13轮的单次最优模块选择器已给下界。这里 \(C_0\) 指开始迭代时的码；如先补齐初始仿射方向，应取补齐后的码。对固定 \(C,B,E\)，(7)还给 \(\kappa(T_Q)/Q\to\kappa(B)\)，旧方向加回相对 \(Q\kappa(B)\) 的增益至多固定量 \(\kappa(C)-\kappa(B)\)。

以上不覆盖从当前码之外新增其他合法方向、非三角同步构造、随规模变化且容量密度无界的种子，亦不控制一般非平坦谱。**一般A28仍开放**；本轮关闭的是旧方向加回与反复仿射恢复这两种明确构造路线。
