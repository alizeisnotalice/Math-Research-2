# 第11轮：任意环境下四分之一 Fourier 网格的锐七维界

## 定理与范围

设F:F₂ⁿ→{−1,1}，使用均匀概率归一化Fourier变换。若所有系数均属于(1/4)ℤ，则

\[
 \operatorname{affdim}\operatorname{supp}\widehat F\le6,
 \qquad \dim\operatorname{span}\operatorname{supp}\widehat F\le7.
\]

两个界均能取等。不限制环境维数，不要求谱平坦。这里控制的是Fourier系数的二进分母，而非实多项式次数或F₂代数次数；没有将此有限粒度结论外推到任意粒度的A28。

## 1. 平坦十六点谱

先设支撑S有16点且每个非零系数为εᵢ/4，εᵢ=(−1)^{fᵢ}。在S上取评价码

\[
 C=\langle1,\alpha\mapsto\alpha_1,\ldots,\alpha\mapsto\alpha_n\rangle
 \le\mathbb F_2^{16}.
\]

它的维数D=affdim S+1，而线性谱维至多D。每个c∈C是输入平移与整体变号对谱符号的作用。由Fourier反演，在任一输入、特别是在零处，

\[
 \sum_{i=1}^{16}(-1)^{f_i+c_i}\in\{-4,4\},
 \qquad \operatorname{wt}(f+c)\in\{6,10\}. \tag{1}
\]

定义带符号重量T(c)=Σᵢcᵢ(−1)^{fᵢ}。它等于wt(f+c)−wt(f)，故4整除T(c)。由模二约化，C所有字为偶重量；再由

\[
 T(c+d)=T(c)+T(d)-2\sum_i(-1)^{f_i}c_i d_i
\]

得c·d=0。因此C⊂C⊥，D≤8。

反设D=8，则C=C⊥。令K₄(j)=Σᵢ(−1)ⁱ binom(j,i)binom(16−j,4−i)，直接整数计算给K₄(6)=K₄(10)=−20。将四坐标字符对C平均，得到自包含恒等式

\[
 \begin{aligned}
 -20
 &=\mathbb E_{u\in f+C}K_4(\operatorname{wt}u)\\
 &=\sum_{\operatorname{wt}y=4}(-1)^{f\cdot y}
       \mathbb E_{c\in C}(-1)^{c\cdot y}\\
 &=\sum_{\substack{y\in C^\perp\\\operatorname{wt}y=4}}(-1)^{f\cdot y}
 =\sum_{\substack{y\in C\\\operatorname{wt}y=4}}(-1)^{f\cdot y}.
 \end{aligned} \tag{2}
\]

对每个右端的y，T(y)=4−2|supp(f)∩supp(y)|被4整除，故f·y=0。右端等于C的重量4字数，非负，矛盾。于是D≤7。

## 2. 非平坦系数用重复坐标处理

一般地写aα=4F̂(α)∈ℤ。Parseval给Σαaα²=16。将每个非零频率α重复|aα|次，每次带符号sgn(aα)，得到长度

\[
 L=\sum_\alpha|a_\alpha|\le\sum_\alpha a_\alpha^2=16
\]

的带符号频率多重集。由Σαaα=4F(0)，L为偶数。重复坐标不会改变仿射函数在原支撑上的线性关系，所以其评价码仍有维数affdim S+1。

现在(1)变为wt(f+c)∈{(L−4)/2,(L+4)/2}。两种可能重量之差仍为4，因此同样的带符号重量论证给评价码自正交，D≤L/2。如果L≤14，立得D≤7。如果L=16，则Σα(aα²−|aα|)=0，迫使全部非零aα为±1，回到上一节。常数函数的单点谱也包含在本论证中。

## 3. 锐性

取七个独立的±1变量r,u,v,y₀,y₁,y₂,y₃，令

\[
 F=\frac{1+r}{2}\left(\frac{1+u}{2}y_0+\frac{1-u}{2}y_1\right)
  +\frac{1-r}{2}\left(\frac{1+v}{2}y_2+\frac{1-v}{2}y_3\right).
\]

每次恰选择一个叶变量，故F为Boolean。展开得16个互异字符，系数皆±1/4。谱点之差张成三个选择变量及三个独立叶差，仿射维数为6；线性张成另增加一个叶方向，维数为7。

这使用根节点一个选择变量、两个子节点各自独立的选择变量。共享两位地址的普通寻址函数只有六维，不能作为此处的七维取等例。附带整数脚本对全部128输入核验该例，同时核验(2)中的Krawtchouk数值；任意n的上界由上述证明承担。
