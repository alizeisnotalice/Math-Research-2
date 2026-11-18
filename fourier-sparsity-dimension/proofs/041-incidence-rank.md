# 共享频率的出现顺序与秩扣减

沿用[粘合定理](041-collision-gluing.md)的 $p_i,S_i,a_i,s_i,L,B_i$。本节允许一个共享标签牵涉任意多块，每块也可同时参与多个不同的共享标签。令 $H$ 是至少两个 $B_i$ 都含有的商频率集合，$H_i=H\cap B_i$，先假设每个 $H_i$ 中的标签在 $S_i$ 中恰有一个原像。定义 $\gamma_i=(a_i+2)/s_i$。若每块满足 $a_i+2\le2\sqrt{s_i}$，且对每个 $A\subset H_i$ 有
$$\dim\operatorname{span}\{(\beta,1):\beta\in A\}\ge\gamma_i |A|,\tag{1}$$
则 $a+2\le2\sqrt{s}$。这里不要求将有碰撞的叶子分成互不相交的重复组。特别地，只要每个碰撞叶都非字符且它所参与的共享标签在 $L^*$ 中仿射无关，(1) 即成立：$s_i\ge4$ 时 $\gamma_i\le1$。

**频谱。** 对 $\beta\in H$，令 $P_\beta=\sum_{i:\beta\in H_i}p_i$，$w_{i\beta}=p_i/P_\beta$。该单一商频率带在物理空间定位于这些块，非零，因而至少有 $1/P_\beta$ 个频率。不同商频率带互不相交。每条私有局部频率的完整 $1/p_i$ 环境提升均保留。利用 $\sum_{i:\beta\in H_i} w_{i\beta}^2/p_i=1/P_\beta$ 得
$$ s\ge\sum_i T_i/p_i,\qquad T_i=s_i-|H_i|+\sum_{\beta\in H_i}w_{i\beta}^2 .\tag{2}$$
把长度 $s_i$ 的局部向量的私有分量取1，共享分量取 $w_{i\beta}$，与全1向量作 Cauchy–Schwarz，得
$$\sqrt{s_i T_i}\ge s_i-\sum_{\beta\in H_i}(1-w_{i\beta}).\tag{3}$$

**秩。** 给每块独立分配速率为 $p_i$ 的指数时钟，按首次响铃时间排序。对叶 $i$ 的共享标签 $\beta$，某个其他含此标签的叶在它之前出现的概率恰为 $1-w_{i\beta}$；这只用了参与该标签的时钟中 $i$ 最先响铃的概率。记 $A_i$ 为这些已经出现的标签。令 $V_i=\operatorname{span}\{(\beta,1):\beta\in B_i\}$。在任一顺序下，添加 $V_i$ 时，与此前空间的交含有 $A_i$ 的齐次张成空间。故用上一证明的法向秩界与 (1)，再对顺序取期望，得
$$a+2\le\sum_i\left(a_i+2-\gamma_i\sum_{\beta\in H_i}(1-w_{i\beta})\right)
\le\sum_i\gamma_i\sqrt{s_i T_i}\le2\sum_i\sqrt{T_i}\le2\sqrt s.$$
最后一步使用 (2) 和 $\sum_i p_i=1$。这些随机顺序只用于证明确定性不等式。重复次数、不同标签的交叠形状以及各块相对质量均不限。

**多个局部原像的扩展。** 令 $r_{i\beta}=|\{\alpha\in S_i:\alpha|_L=\beta\}|$，$e_i=a_i+1-\dim V_i$。保留局部输入界，并把 (1) 替换为下列两项足够条件：
$$e_i\ge\gamma_i\sum_{\beta\in H_i}(r_{i\beta}-\sqrt{r_{i\beta}}),\qquad
\dim\operatorname{span}\{(\beta,1):\beta\in A\}\ge\gamma_i\sum_{\beta\in A}\sqrt{r_{i\beta}}\quad(A\subset H_i).$$
事实上单个共享带的物理支撑仍至多 $P_\beta$。它不恒零，因为任一参与块上的局部投影有非空谱，不同块物理支撑不交。在 (2) 中把 $|H_i|$ 换为 $\sum r_{i\beta}$；在 (3) 的向量中给该标签的 $r_{i\beta}$ 个坐标都填 $w_{i\beta}/\sqrt{r_{i\beta}}$。需要扣的局部量遂为 $\gamma_i\sum(r_{i\beta}-w_{i\beta}\sqrt{r_{i\beta}})$。其确定性部分由 $e_i$ 支付，余项由随机顺序下已有共享标签的秩支付。其余证明相同。

**两个方向交叉重复的族。** 令 $m=2^k\ge2$，以两个 $k$ 位地址 $(i,j)$ 选取 $t=m^2$ 个等质量块。数据空间取独立字符方向 $e_i,f_j,u_{ij}$。第 $(i,j)$ 块输出
$$\chi_{e_i}(y)(-1)^{((e_i+f_j)\cdot y)((e_i+u_{ij})\cdot y)}.$$
其四个数据谱点为 $e_i,f_j,u_{ij},e_i+f_j+u_{ij}$。前两点的局部系数均为 $1/2$，后两点为 $1/2,-1/2$。行标签 $e_i$ 和列标签 $f_j$ 各被 $m$ 块复用，全部块经共享标签相连，而任一商频率都未出现在所有块中。每叶两个共享标签仿射无关。每个行或列共享带在地址上是一个仿射平坦的常数倍，其谱点恰有 $m$ 个；两种私有带各贡献 $t^2$ 个频率。故 $s=2t^2+2m^2$。共同频率的真实支撑由这些平坦的变换给出，大量候选系数确实相消。公式对任意 $k$ 成立。

**共享标签的活动域恰为仿射块时。** 还有另一种处理多个原像的方式。对于有多个原像的共享标签，若所有参与块的并 $U_\beta$ 为仿射平坦，在该平坦上考虑相应带信号的内在谱。限制到任一参与子块只会合并谱点，因此内在谱至少有 $\max_i r_{i\beta}$ 点；全空间谱是其完整提升，便有
$$s_\beta\ge\max_i r_{i\beta}/P_\beta\ge\sum_i r_{i\beta}w_{i\beta}^2/p_i.$$
其他只有单原像的共享标签仍可使用原来的估计。现在每一局部原像都能在 (3) 的向量中填 $w_{i\beta}$。一个足够条件是，对每个叶 $i$ 和每个 $A\subset H_i$，投到 $A$ 的全部局部原像的齐次秩至少为 $\gamma_i\sum_{\beta\in A}r_{i\beta}$。因为这些原像提升空间的维数不超过投影总核维数 $e_i$ 加上 $A$ 的齐次秩，故对 $A=A_i$ 取期望后，局部核与旧标签交秩之和足以扣除 $\gamma_i\sum r_{i\beta}(1-w_{i\beta})$，得到同一结论。
