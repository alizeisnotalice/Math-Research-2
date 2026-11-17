# 双平坦 signed indicator 的可投影性：精确对偶归约

设 \(G=\mathbb F_2^{2k}\)、\(q=2^k\)，并令
\[
h=1_A(-1)^f,\quad |A|=2q,\quad W_h=q\varepsilon1_S,
\quad |S|=2q,\quad\varepsilon:S\to\{\pm1\}. \tag{1}
\]
这是合法 \(2q\) 翻号块的单块 plateaued 数据。

## 1. 双平坦数据真正控制的量

未归一化 Walsh 变换给两个恒等式：
\[
\sum_xh(x)h(x+a)=\sum_{\xi\in S}(-1)^{\xi(a)},\tag{2}
\]
\[
W_{1_A}(\delta)=W_{h^2}(\delta)
=\sum_{\xi\in S\cap(S+\delta)}
  \varepsilon_\xi\varepsilon_{\xi+\delta}.\tag{3}
\]
第二式用 \(W_{h^2}=|G|^{-1}(W_h*W_h)\) 及 \(q^2=|G|\)。二者都是**带符号**
计数，不能和无符号差对数 \(|A\cap(A+a)|\) 混同。

## 2. 可投影性的精确 Fourier 判别

对 \((k-1)\)-子空间 \(H\le G\)，令 \(K=H^\perp\)。下列等价：

1. \(H\cap(A+A)=\{0\}\)；
2. \(A\) 在 \(G/H\) 上单射；
3. \(A\) 是全部 \(H\)-陪集的一个 transversal；
4. \(1_A*1_H=1_G\)；
5. \(\dim K=k+1\) 且 \(W_{1_A}(\delta)=0\) 对全部 \(0\ne\delta\in K\)。

前三项是陪集差的直接等价。第四项在每点计数该陪集中的 \(A\)-点数。又
\(W_{1_H}=|H|1_K\)、\(W_{1_G}=|G|1_{\{0\}}\)，故 Fourier 变换给 4 与 5。
结合 (3)，原问题精确成为：能否在 signed intersection sum 的零集内找到一个
\((k+1)\)-维线性子空间（除零点外）。

## 3. 一个不能跳过的相消障碍

在 \(\mathbb F_2^4\) 令 \(A=\{x:x_3=0\}\)，
\(h=1_A(-1)^{x_1x_2}\)。则 \(|\operatorname{supp}h|=|\operatorname{supp}W_h|=8\)，
且全部非零 Walsh 系数模为4。但对 \(a=e_1\)，
\[
a\in A+A,\qquad |A\cap(A+a)|=8,\qquad
\sum_xh(x)h(x+a)=0.\tag{4}
\]
所以 (2) 的零不能证明不存在真实差值；signed autocorrelation 的完全相消可以发生。
该例本身来自真实合法翻号：取 f=x₁x₂+x₀x₃，这是两个独立二次对之和，故 bent；而 1_A=1+x₃，加入这个仿射相位仍 bent。限制到 A 后恰得上述 h。本机已用整数 Walsh 逐项核验这两个相位，故相消障碍并非仅发生于抽象双平坦数据。

该例仍可投影（核 \(\langle e_3\rangle\)），故它不是反例，只排除以 (2) 代替
无符号差集的证明路线。整数核验见[复现源码](../computation/round020/verify_signed_autocorrelation_obstacle_k2.py)与[实际结果](../computation/round020/signed_autocorrelation_obstacle_k2.json)。

## 4. 当前严格状态

(1)--(3) 并未推出 \(\{\delta:W_{1_A}(\delta)=0\}\) 含 \((k+1)\)-子空间。因此这里没有
一般可投影性证明，也没有把六变量 RM 枚举推广至高维。后续必须控制 (3) 的零集，或加入
MM/共同码等额外结构；不能只引用近不确定性取等的支撑大小。
