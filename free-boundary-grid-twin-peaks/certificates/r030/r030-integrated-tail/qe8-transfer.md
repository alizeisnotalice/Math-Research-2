# QE8 的条件传递：top-two 律与两峰条件律

## 范围与前提

这是一条**条件定理**，不证明 QE8，也不证明任一模型下的远峰消失。

令 $\mathcal B$ 是 R028 的 uniform top-two 切片，$\mathbb P_B$ 是其均匀律；
稳定提升／第二峰插入所给出的权重变量记为 $L\in\{1,\ldots,M\}$。对在该
提升中保持的峰几何可测非负函数 $f$，精确变换为

$$
 \mathbb E_2 f=\frac{\mathbb E_B[L f]}{\mathbb E_B L};
 \qquad
 \mathbb P_2(E)=\frac{\mathbb E_B[L\mathbf1_E]}{\mathbb E_B L}.
 \tag{1}
$$

特别地，$E=\{\rho_n>\varepsilon n\}$ 是此类事件：它只取决于两峰及其原图距离。
令

$$
 T_t=\#\{b\in\mathcal B:L(b)\ge t\},\qquad B=|\mathcal B|,
 \tag{2}
$$

并把 $T_t=0$ 延拓到 $t>M$。这里所说的 QE8 是全局（对 $v$ 求和以后）条件

$$
 \tag{QE8}T_t\le {8B\over t(t+1)}\quad(2\le t\le M).
$$

它不是逐孔不等式；中心孔 $n=7,t=4$ 的精确 CRT 数据已经反驳了相应的逐孔加强。

## 条件定理

**命题。** 若 (QE8) 对给定的方格成立，则

$$
 \mathbb E_B L\le5.
 \tag{3}
$$

对于每个保根的峰几何事件 $E$，令 $p=\mathbb P_B(E)$。那么

$$
 {p\over5}\le\mathbb P_2(E)\le8\sqrt p.\tag{4}
$$

令 $D_n=\|x_{\rm peak}-y_{\rm peak}\|_2^2/n^2\in[0,2]$ 为项目既定的
平方欧氏峰距归一化，并定义

$$
 Q_B=\mathbb E_BD_n,\qquad Q_2=\mathbb E_2D_n,
$$

则还具有

$$
 {Q_B\over5}\le Q_2\le8\sqrt{2Q_B}.\tag{5}
$$

因此，若 QE8 对一列 $n$ 一致成立，则该列上

$$
 \mathbb P_B(\rho_n>\varepsilon n)\to0
 \quad\Longleftrightarrow\quad
 \mathbb P_2(\rho_n>\varepsilon n)\to0
 \qquad(\text{每个固定 }\varepsilon>0),
 \tag{6}
$$

且 $Q_B\to0\Longleftrightarrow Q_2\to0$。后一陈述与
$\mathbb P_2(\rho_n>\varepsilon n)\le2Q_2/\varepsilon^2$ 结合，给出通常的矩收尾。

**证明。** 尾和公式和 (QE8) 给出

$$
 \mathbb E_B L={1\over B}\sum_{t\ge1}T_t
 \le1+8\sum_{t=2}^\infty{1\over t(t+1)}=5.
 \tag{7}
$$

对任意整数 $K\ge1$，再次用尾和公式，

$$
\begin{aligned}
 \mathbb E_B[L\mathbf1_{\{L>K\}}]
 &=K{T_{K+1}\over B}+\sum_{t\ge K+1}{T_t\over B}\\
 &\le {8K\over(K+1)(K+2)}+{8\over K+1}
 ={16\over K+2}.
\end{aligned}\tag{8}
$$

因为 $L\ge1$，(1)、(3)、(8) 蕴含

$$
 {p\over5}\le\mathbb P_2(E)
 \le Kp+{16\over K+2}.
 \tag{9}
$$

当 $p>0$，取 $K=\lfloor4/\sqrt p\rfloor$。由于 $0<p\le1$，$K\ge4$，且

$$
 Kp\le4\sqrt p,\qquad {16\over K+2}<4\sqrt p,
$$

故得到 (4)。当 $p=0$，(1) 直接给出 $\mathbb P_2(E)=0$。

对 $D_n$ 重复同一分割，并使用 $0\le D_n\le2$，有

$$
 {Q_B\over5}\le Q_2
 \le KQ_B+{32\over K+2}.\tag{10}
$$

当 $Q_B>0$，正确的平衡取值是

$$
 K=\left\lfloor{8\over\sqrt{2Q_B}}\right\rfloor,
$$

而不是事件版本的取值。此时 $KQ_B\le4\sqrt{2Q_B}$ 及
$32/(K+2)<4\sqrt{2Q_B}$，于是 (5) 成立。$Q_B=0$ 时由非负性和 (1) 直接得到
$Q_2=0$。

## 限制

(4)--(6) 只说明：一旦独立地证明 top-two 律的远峰概率/矩趋零，就可以在 QE8
下传递到原来的 $K_n=2$ 条件律，反向亦然。QE8 的归一化控制本身既不证明
$\mathbb P_B(\rho_n>\varepsilon n)\to0$，也不证明 $Q_B\to0$；因此它本身不能闭合
原问题。

有限 $n$ 的 QE8 检验或逐孔反例都不改变上述条件定理的量词。
