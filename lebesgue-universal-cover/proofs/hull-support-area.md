# 两个凸体并凸包的面积及三个定宽体的中位支撑

设 \(K,L\) 为平面凸体，已包含各自给定的刚体摆放。写 \(h_K,h_L\) 为以同一原点计算的支撑函数，\(S_K=h_K+h_K''\) 为圆周上的曲率测度，\(K\vee L=\operatorname{conv}(K\cup L)\)。积分均取整个单位圆周，角度测度为弧度。对周期 Lipschitz 函数定义
\[
 \mathcal A(h)=\frac12\int_0^{2\pi}(h^2-h'^2)\,d\theta .
\]
对凸体支撑函数，该式等于面积。可由光滑严格凸体的参数化 \(p=hu+h'u'\) 与面积积分先得光滑情形，一般情形以圆周正核卷积支撑函数，支撑函数在 \(H^1\) 收敛，凸体在 Hausdorff 距离下收敛。

**二体面积恒等式。** 对任意两个摆放，
\[
 |K\vee L|=\frac{|K|+|L|}{2}
       +\frac14\int |h_L-h_K|\,d(S_K+S_L).
 \tag{1}
\]
令 \(f=h_L-h_K,\ g=f_+\)，则 \(h_{K\vee L}=h_K+g\)。展开二次型并做弱分部积分，得
\[
 2(|K\vee L|-|K|)=2\int g\,dS_K+\int(g^2-g'^2)
                  =\int g\,d(S_K+S_L).
\]
第二式用 \(g f=g^2\) 与 \(g'f'=g'^2\) 几乎处处成立。交换 \(K,L\) 后相加即得 (1)。弱分部积分对 Lipschitz 测试函数 \(g\) 成立，因原支撑函数的二阶弱导数为有限测度，可先对测试函数卷积并用 \(H^1\) 收敛。

固定相对角、手性，平移 \(L\) 只在 \(h_L-h_K\) 中添加 \(t\cdot u\)，其中 \(u=(\cos\theta,\sin\theta)\)。设 \(\mu=S_K+S_L\)，\(b=h_L-h_K\)，则平移后最小面积等于 \(|K|\) 加上
\[
 \frac12\max\left\{\int z b\,d\mu:\ 0\le z\le1,\quad \int zu\,d\mu=0\right\}.
 \tag{2}
\]
正部的对偶表达给出其中的下界。反过来，实凸函数 \(t\mapsto\int[b+t\cdot u]_+\,d\mu\) 连续且强制，故有极小点。强制性可由并凸包包含 \(K\) 的一个固定内圆盘和 \(L+t\) 的一个指定点直接得到。在极小点，\(b+t\cdot u\) 的正区取 \(z=1\)、负区取 \(z=0\)、零区可取 \(z\in[0,1]\) 构成的向量集合 \(\{\int zu\,d\mu\}\) 是闭凸集。它各方向的支撑值就是目标的单边方向导数。最小性使它包含零，故可选一个平衡 \(z\) 达到等号。对多边形，\(\mu\) 是两体所有边的外法向原子，权重为边长；平行边保留重数。这给出只在两个平移变量上求解的有限分段线性凸函数。固定相对朝向的平移凸性已有经典工作，参见[Horváth 与 Lángi](https://arxiv.org/abs/1302.6164)。

若 \(K,L\) 的宽度都为 \(w\)，则 \(h(\theta)+h(\theta+\pi)=w\) 且 \(S_K(\theta)+S_K(\theta+\pi)=w\,d\theta\)。由于 \(|h_L-h_K|\) 为 \(\pi\) 周期，(1) 化为
\[
 |K\vee L|=\frac{|K|+|L|}{2}+\frac w4\|h_L-h_K\|_{L^1(d\theta)}.
 \tag{3}
\]
两输入的独立摆放可固定其中一个体，另一个体的第一谐波平移及正反手旋转任取。所得的最小 \(L^1\) 距离是同宽凸体合同类上的距离：刚体群在支撑差的 \(L^1\) 范数上等距，三角不等式随之成立；有界的极小序列中，平移由 \(\|t\cdot u\|_1=4|t|\) 有界，正交部分紧，距离为零便确有合同。

**三体的中位函数。** 对已摆放的三个同宽体 \(K_i\)，令 \(H=\bigvee_iK_i\)，\(m(\theta)=\operatorname{med}_i h_i(\theta)\)。则
\[
 P_H-\pi w=\frac14\sum_{i<j}\|h_i-h_j\|_1,\qquad
 2|H|=\sum_{i<j}|K_i\vee K_j|-\mathcal A(m).
 \tag{4}
\]
对每个方向，三项两两绝对差之和为其极差的两倍。用 \(\min_i h_i(\theta)=w-h_H(\theta+\pi)\) 与 \(P_H=\int h_H\)，积分即得首式。排序保持三项平方和，也在几乎处处保持导数平方和。在取值相等的集合上，相应两函数的导数几乎处处相等。因此
\(\sum |K_i|=\mathcal A(h_H)+\mathcal A(m)+\mathcal A(w-h_H(\cdot+\pi))\)，末项等于 \(|H|+\pi w^2-wP_H\)。代入 (3) 和首式给出第二式。

中位数本身只需是 Lipschitz 函数。若 \(m-w/2\) 的实 Fourier 系数为 \(a_k,b_k\)，反周期性和 Parseval 恒等式给出
\[
 \mathcal A(m)=\frac{\pi w^2}{4}
 -\frac{\pi}{2}\sum_{\substack{k\ge3\\ k\ {\rm odd}}}(k^2-1)(a_k^2+b_k^2).
\]
故 (4) 同时给出只用三项两体面积的一般下界和精确的余项。二次域内的具体双输入解见[三接触完全证明](../computation/two-polygon-exact/proof.md)。
