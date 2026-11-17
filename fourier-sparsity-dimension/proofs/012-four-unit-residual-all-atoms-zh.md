# \(t=4\) 残余的全部原子数：\(\dim G\le6\) 定理

## 定理

设 \(X=\mathbb F_2^{2k}\)、\(q=2^k\)（\(k\ge2\)），\(R\subseteq X\)、\(D\subseteq\widehat X\)，且 \(|R|=|D|=4q\)。令 \(1_R\in G\le\mathbb F_2^R\)，
并设所有
\[
 1_R(-1)^{f+g}\quad(g\in G)
\]
共同 Walsh 支撑于 \(D\)，其上幅值均为 \(q\)。再假定每个
\(g\notin\{0,1_R\}\) 的重量属于 \([3q/2,5q/2]\)，且没有重量 \(q\) 的字。
则
\[
 \dim G\le6.
\]

证明使用 Gopalan--O'Donnell--Servedio--Shpilka--Wimmer, Theorem 8.1 的
Boolean Fourier-sparsity granularity 引理：若 Boolean 函数的 Fourier sparsity
为 \(s>1\)，则所有 Fourier 系数属于
\(2^{-(\lfloor\log_2s\rfloor-1)}\mathbb Z\)；sparsity 1 时 granularity 为 0。
故 sparsity 至多 15 时系数属于 \(\frac14\mathbb Z\)。
原始定理见[作者原文第21页](https://www.cs.columbia.edu/~rocco/Public/4-27-2011-gossw.pdf)。所有Walsh变换不归一化；参数行的Fourier系数按均匀概率归一化。

## 1. 原子、概率与 \(r\ge13\)

取补空间 \(H\) 使 \(G=H\oplus\langle1_R\rangle\)，令
\(m=\dim H\)。按 \(H\) 的评价模式把 \(R\) 分为 \(r\) 个原子，大小为 \(w_i\)。
对参数 \(h\in H\) 作字符平均，可从共同部分 Walsh 变换中抽出每个原子的变换 \(P_i\)，故其频率支撑包含于 \(D\)。原子函数的物理支撑大小为 \(w_i\)，而 \(|X|=q^2\)，不确定性给 \(w_i|D|\ge q^2\)，即 \(w_i\ge q/4\)，故 \(r\le16\)。令
\[
 \mu_i=\frac{w_i}{4q}.
\]
对 \(0\ne h\in H\)，重量条件给
\[
 |\widehat\mu(h)|=\left|1-\frac{\operatorname{wt}(h)}{2q}\right|\le\frac14.
\]
若反设 \(\dim G\ge7\)，则 \(m\ge6\)，Parseval 给
\[
 \sum_i\mu_i^2
\le2^{-m}\left(1+\frac{2^m-1}{16}\right)
\le\frac1{16}+\frac{15}{1024}=\frac{79}{1024}.
\tag{1}
\]
另一方面 \(\sum_i\mu_i^2\ge1/r\)，故 \(r\ge13\)。

对 \(\xi\in D\)，原子部分 Walsh 和 \(P_i(\xi)\) 给出一个 Boolean Fourier 行，
\[
 B_\xi(h)=\sum_i\frac{P_i(\xi)}q(-1)^{h_i}.
\tag{2}
\]
当 \(r\le15\)，granularity 输入说明每个系数为 \(1/4\) 的整数倍（常数行更强）。
Parseval 因而强迫每行至少一个系数绝对值 \(\ge1/2\)。令
\[
 \mathcal B=\{i:w_i\ge q/2\}.
\]
逐行有
\[
 \sum_{i\in\mathcal B}(P_i(\xi)/q)^2\ge\frac14.
\]
对 \(\xi\in D\) 求和并用每个原子的 Plancherel
\[
 \sum_{\xi\in D}(P_i(\xi)/q)^2=w_i
\]
得到
\[
 \sum_{i\in\mathcal B}w_i\ge q.
\tag{3}
\]

## 2. 排除 \(r=13,15\)

若 \(r=15\)，最小原子大小账本显示至多一个大原子，且它至多 \(q/2\)，与
(3) 矛盾。

若 \(r=13\)，只要存在一个大原子，\(\mu_i\ge1/8\)。把其余质量在 12 个原子上
用 Cauchy 最小化；函数a²+(1−a)²/12在a≥1/8上递增，故得
\[
 \sum_i\mu_i^2\ge\frac1{64}+\frac{(7/8)^2}{12}
=\frac{61}{768}>\frac{79}{1024},
\]
与 (1) 矛盾；但每行的 \(1/2\) 系数又强迫存在大原子。故 \(r=13\) 不可能。

## 3. \(r=14\) 的等号结构

由原子大小账本和 (3)，此时恰有两个大原子，且
\[
 w_{b_1}=w_{b_2}=q/2,\qquad
 w_i=q/4\ (i\in S),\quad |S|=12.
\tag{4}
\]
更细地，对所有频率求和，
\[
 \sum_{\xi\in D}\sum_{i\in\{b_1,b_2\}}(P_i(\xi)/q)^2
=w_{b_1}+w_{b_2}=q.
\]
除以 \(|D|=4q\)，每行大原子能量的平均恰为 \(1/4\)。它已逐行不小于
\(1/4\)，故每行取等。结合系数属于 \(\frac14\mathbb Z\)、小原子系数的绝对值
至多 \(1/4\)，每行必为：

* 恰有一个大原子系数为 \(\pm1/2\)，另一个为 0；
* 全部 12 个小原子系数均为 \(\pm1/4\)。

固定一行，记活跃大原子为 \(b\)，其符号为 \(\varepsilon_b\)，小原子符号为
\(\varepsilon_i\)。对 \(h\in H\) 令
\[
 T(h)=\sum_{i\in S}\varepsilon_i h_i.
\]
Boolean 值差为
\[
 B_\xi(0)-B_\xi(h)=\frac12T(h)+\varepsilon_bh_b\in\{0,\pm2\}.
\]
因此
\[
 T(h)+2\varepsilon_bh_b\equiv0\pmod4.
\tag{5}
\]
把 (5) 用于 \(h,k,h+k\)，并用
\(h_b+k_b-(h+k)_b=2h_bk_b\)，得
\[
 T(hk)\equiv0\pmod2.
\tag{6}
\]
因符号模 2 都为 1，(6) 表明 \(H\) 在 12 个小原子上的投影自正交；
(5) 还表明其每个投影字偶重，所以添入全 1 字后，整个 \(G\) 的小原子投影仍自正交。

该投影的核支撑在两个大原子上，物理总大小为 \(q\)。非零核字既不是 \(1_R\)，
又有重量至多 \(q<3q/2\)，违背重量假设。因此投影单射。长度 12 自正交码维数
至多 6，故 \(\dim G\le6\)。

## 4. \(r=16\)

此时所有原子都恰为 \(q/4\)。第 11 轮的饱和 16 原子论证适用：单行的
\(\pm1/4\) Fourier 系数导出 16 坐标自正交码，第二字符矩与重量 6
Gram 矩阵矛盾排除维数至少 7。因此仍有 \(\dim G\le6\)。

结合 \(13\le r\le16\)，定理成立。第 10 轮的实际 \(q=8,\dim G=6\) 残余表明
常数 6 不能普遍降低。

## 真实 bent 残余满足这些假设

设X=F₂^{2k}，f₀ bent，G中的函数都支撑于R，1_R∈G，且f₀+G全部bent，|R|=4q。记Q=W(1_{Rᶜ}(−1)^{f₀})、P_g=W(1_R(−1)^{f₀+g})。分别加入0和1_R得Q±P_g∈{±q}，于是Q P_g=0、Q²+P_g²=q²。因此D={ξ:Q(ξ)=0}与g无关，在D上P_g=±q；Parseval给|D|=|R|=4q。

若G没有重量q字，则每个非零g属于RM(k,2k)，其重量至少3q/2。这使用bent次数界及Reed–Muller第二重量定理，准确适用范围与来源见[第9轮](009b-residual-budget-zh.md)。对g和1_R+g同时应用即得所需区间[3q/2,5q/2]。所以本定理适用于所有这样的真实bent残余，而不仅仅是饱和原子情形。

特别地，在包含极大且全部最小组件已取尽的空间中，若剩余t=4，由第9轮联合秩分解得到dimV≤(q−1)+(q−4)+6=2q+1。这里没有声称该全局界取等，或证明一般t的O(t)残余预算。第10轮真实dimG6例子证明残余上界六本身不能普遍降低，但该例的全局空间并非包含极大。

## 粒度引理的归纳说明

这是GOSSW的经典结论，不作为本轮新定理。对环境维数归纳：若支撑已满，归一化Walsh系数的分子为偶数，分母为2ⁿ。若未满，任取零系数α和非零系数β，将函数限制到(α+β)·x=0。商频率中的这对系数相加仍为原来的β系数，限制函数仍Boolean，支撑数不增且环境维数减一。对每个β分别应用归纳即可得到floor(log₂s)−1；单字符单独为整数系数。由此s≤15时所有系数属于(1/4)ℤ。
