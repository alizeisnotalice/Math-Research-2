# Flat coset 的带符号交积与四阶缺额审计

## 设定

令 \(q=2^k\), \(s=q^2\)。令 \(C\le\mathbb F_2^s\) 含全 1 字，并设
\[
 \operatorname{wt}(f+c)\in\{(s-q)/2,(s+q)/2\}\qquad(c\in C).
\tag{1}
\]
另设 \(d_{\min}(C)\ge q\)。本记录只从这些平坦 coset 条件推出编码约束；它不是
一般 A28 的维数结论。

## 1. 全部带符号交积的 2-adic 阶梯

令 \(\varepsilon_x=(-1)^{f(x)}\)，并定义
\[
 A(c)=\sum_x\varepsilon_x(-1)^{c(x)}\in\{\pm q\}.
\]
对 \(j\ge1\) 个码字 \(c_1,\ldots,c_j\)，令
\[
 I_\varepsilon(c_1,\ldots,c_j)
 =\sum_x\varepsilon_x\prod_{\nu=1}^j c_\nu(x).
\]
有限差分恒等式给
\[
 \Delta_{c_1}\cdots\Delta_{c_j}A(0)
 =(-2)^j I_\varepsilon(c_1,\ldots,c_j).
\tag{2}
\]
左边是 \(2^j\) 个 \(\pm q\) 的带符号和。对 \(j\ge1\)，模 2 每项皆为 1，
故此和总被 \(2q\) 整除。因此
\[
 I_\varepsilon(c_1,\ldots,c_j)
 \equiv0\pmod {2^{\max(k-j+1,0)}}.
\tag{3}
\]
这里重复因子可按 \(c^r=c\) 删除。

不带符号交积与带符号交积模 2 相同。因此任意至多 \(k\) 个（不同）码字的
坐标乘积都有偶重量。特别地：

* \(k\ge1\)：所有码字偶重；
* \(k\ge2\)：\(C\subseteq C^\perp\)（自正交）；
* \(k\ge3\)：\(C\) triorthogonal；
* 一般地，得到通常意义的 \(k\)-orthogonal 奇偶条件。

注意 (3) 是**带符号**交积整除；除模 2 结论外，不能自动把它改写成无符号高阶
重量整除。

## 2. 四阶 Krawtchouk 缺额

以下假定 \(k\ge2\)，令 \(G=C^\perp\)，则 \(C\subseteq G\)。对每个 \(j\)，coset
MacWilliams 恒等式给
\[
 \sum_{\substack{g\in G\\\operatorname{wt}(g)=j}}(-1)^{f\cdot g}
 =\frac1{|C|}\sum_{c\in C}K_j(\operatorname{wt}(f+c)).
\tag{4}
\]
取 \(j=4\)。因 \(K_4(s-w)=K_4(w)\)，(1) 使右边化为
\[
 K_4((s-q)/2).
\]
利用
\[
 K_4(w)=\frac{t^4-(6s-8)t^2+3s(s-2)}{24},\qquad t=s-2w,
\]
并代入 \(t=q,q^2=s\)，得到
\[
 K_4((s-q)/2)=-\frac{s(s-1)}{12}.
\tag{5}
\]
若 \(B_4=|\{g\in G:\operatorname{wt}(g)=4\}|\)，则由 (4)--(5)
\[
 B_4\ge\frac{s(s-1)}{12}.
\tag{6}
\]

## 3. 按 \(C\)-陪集的装箱

同一 \(C\)-陪集内两个不同重量 4 字之差属于 \(C\setminus\{0\}\)，故其重量至少
\(q\)。

* 若 \(q=8\)，两个 4-集的对称差为 \(8-2t\)，故必须 \(t=0\)。每个陪集至多有
  \(s/4=16\) 个重量 4 字。
* 若 \(q\ge16\)，两个重量 4 字的差最多为 8，小于 \(q\)，故每个陪集至多一个。

零陪集 \(C\) 本身没有重量 4 字。写 \(D=\dim C\)，则
\[
 |G/C|=2^{s-2D}.
\]
由 (6) 得：

\[
\begin{array}{c|c}
q=8&
2^{s-2D}-1\ge \frac{s(s-1)}{12}\frac4s=21,\quad D\le29,\\[2mm]
q\ge16&
2^{s-2D}-1\ge\frac{s(s-1)}{12},\quad
s-2D\ge\left\lceil\log_2\left(1+\frac{s(s-1)}{12}\right)\right\rceil .
\end{array}
\tag{7}
\]

第一行用 \(s=64\)。该结论正确但很弱；例如 \(q=16\) 只把 \(s/2=128\) 压到
\(D\le121\)。

## 4. 高阶链条能与不能做什么

(4) 对每个偶 \(j\) 都成立，因而
\[
 B_j(G)\ge |K_j((s-q)/2)|.
\tag{8}
\]
若 \(2j<q\)，最小距离同样使每个 \(C\)-陪集至多有一个重量 \(j\) 字，于是可把
(8) 转成商陪集数下界。

然而这种“单一重量层 + 陪集唯一性”机制本身不可能给出 \(\dim C=O(q)\)：在
\(j<q/2\) 时，即使使用平凡最强上限 \(B_j\le\binom{s}{j}\)，也只有
\[
 \log_2 B_j\le\log_2\binom{q^2}{\lfloor q/2\rfloor}=O(q\log q).
\]
它最多把 \(D\le s/2\) 改善为 \(D\le s/2-O(q\log q)\)，而不是把维数从
\(\Theta(s)\) 降至 \(O(q)\)。重量层 \(j\ge q/2\) 需要新的同陪集 packing 控制，
本记录没有提供。

## 5. 一个连续重量枚举松弛的有理可行证书

取 \(q=8,s=64,D=16\)。考虑实数重量枚举变量 \(A_i\ge0\)，满足
\(A_0=A_{64}=1\)、\(A_i=A_{64-i}\)、奇数重量及 \(1\le i\le7\) 的变量为零，且 \(\sum A_i=2^{16}\)。令
\[
 B_j=2^{-16}\sum_{i=0}^{64}A_iK_j(i),\quad
 T_j=\tfrac12(K_j(28)+K_j(36)),\quad B_j^\pm=(B_j\pm T_j)/2.
\]
加入所有 \(0\le j\le64\) 的 \(B_j^\pm\ge0\)、\(B_2=0\)，以及自正交码的带符号包含约束
\[
 B_j^+\ge A_j\ (j\equiv0\bmod4),\qquad
 B_j^-\ge A_j\ (j\equiv2\bmod4).
\]
后两式来自 \(C\subseteq C^\perp\) 与 \(f\cdot c\equiv\operatorname{wt}(c)/2\pmod2\)；该符号关系由第1节的一阶带符号重量整除得到。

以下是上述**连续松弛**的精确有理可行点，未列变量均为零：
\[
 A_0=A_{64}=1,\quad
 A_8=A_{56}=\frac{503293}{6435},\quad
 A_{26}=A_{38}=\frac{488614}{45},\quad
 A_{30}=A_{34}=\frac{3121790}{143}.
\]
独立 Fraction 程序重算了全部65个重量的变换及不等式；特别地 \(B_4=336\)、\(B_4^+=0\)、\(B_4^-=336\)。证书与逐重量回执见[第12轮计算目录](../computation/round012/RUN_RECORD.md)。

这些数不是整数码字计数，故没有构造出16维合法码，也不是A28反例。它只说明：这里明确列出的连续重量约束，单独不能排除 \(D=16\)。没有核验全部高阶交积约束、整数可实现性或所有更强矩方法；不能把这一松弛可行性扩大为对全部编码方法的否定。
