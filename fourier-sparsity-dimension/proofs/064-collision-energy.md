# 随机核压缩的二阶反集中

沿用同目录核压缩稿的 \(H,h,\beta,n_\beta,M,d,s\)。设 \(d>0\)，且 \(F\) 的 \(d/\sqrt s\) 在全部稀疏度至多 \(s\) 的布尔函数中达到纪录。令 \(u_\beta=n_\beta/s\)，\(P=\max u_\beta\)，\(D_2=\sum_\beta u_\beta^2\)。

**随机桶的平均损失。** 固定一个 \(n\ge2\) 点的非空原频带及其中一点 \(a\)。在随机线性频率映射 \(T^*\) 下，令 \(C_a\) 为同桶其余点数。不同的两个非零二元差向量线性独立，故其落核事件两两独立。因而
\[
 \mu=\mathbb E C_a=(n-1)/M,\qquad
 \mathbb E C_a^2=\mu^2+\mu(1-1/M).
\]
各桶点数 \(j\) 原有 \(j\) 点而像只占 1 点，总删点数为 \(\sum_a C_a/(1+C_a)\)。Cauchy 给
\[
 \mathbb E\frac{C_a}{1+C_a}\ge
 \frac{(\mathbb E C_a)^2}{\mathbb E[C_a(1+C_a)]}
 =\frac{n-1}{2M+n-2}.
\]
频率系数在桶内抵消只会进一步减少像谱大小。\(n=1\) 的平均损失为零。结合确定性桶容量，本次限制满足
\[
s-\mathbb E s'\ge \sum_{\beta:n_\beta>0}
 \max\left\{(n_\beta-M)_+,\
 \frac{n_\beta(n_\beta-1)}{2M+n_\beta-2}\right\}.
\tag{1}
\]
同时，若 \(I\) 为任一商投影谱基，\(B=\sum_{\beta\in I}(n_\beta-1)\)，则核压缩稿给 \(\mathbb E d'\ge d-h-B/(2M)\)。纪录性及 Jensen 给 \(\mathbb E d'\le (d/\sqrt s)\sqrt{\mathbb E s'}\)。由于 \(1-\sqrt{1-x}\ge x/2\)，(1) 左端除以 \(s\) 至多 \(2h/d+B/(Md)\)。

**二阶质量界。** 取不小于 \(n_{\max}\) 的最小二次幂 \(M\)。式 (1) 中的分母至多 \(3M\)，\(B\le s\)，\(M<2n_{\max}\)，所以
\[
 D_2\le \frac1s+\frac3d+\frac{12hP}{d}
 \le \frac4d+\frac{12h}{d}\sqrt{D_2}.
\]
解此二次不等式即得
\[
 \|u\|_2=\sqrt{D_2}\le\frac2{\sqrt d}+\frac{12h}{d}.
\tag{2}
\]
特别，所有维数 \(h=O(\sqrt d)\) 的方向中，均匀谱点的商标签碰撞熵 \(-\log_2 D_2\ge\log_2 d-O(1)\)。任意 \(J\) 个平行陪集的并含至多 \(s\sqrt J(2/\sqrt d+12h/d)\) 个谱点。

取两个独立均匀谱点 \(U,V\)，则对任意仿射 \(h\)-平坦 \(\gamma+H\) 有
\[
 \Pr(U+V\in\gamma+H)=\sum_\beta u_\beta u_{\beta+\pi(\gamma)}
 \le D_2 .
\]
这将结论同时用于有序谱点对的差频分布。高秩假想纪录所须满足的这一结构约束尚未给出一般平方根界。
