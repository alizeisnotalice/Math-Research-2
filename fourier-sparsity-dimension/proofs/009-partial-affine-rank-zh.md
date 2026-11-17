# 部分仿射 packing 的普适法向秩界

## 定理

令 \(X=\mathbb F_2^N\)，并令
\[
A_i=\{x\in X:\langle w,x\rangle=\lambda_i(w)\quad
 \text{对所有 }w\in W_i\}\qquad(i\in I)
\]
为非空、两两不交的仿射平坦，其中 \(W_i\le X^*\) 是法向空间，
\(\lambda_i\in W_i^*\)。若
\[
\dim W_i\le c\quad\text{对所有 }i,
\]
则
\[
 \dim\sum_{i\in I}W_i\le 2^c-1. \tag{T_c}
\]

这里允许各块有不同余维。这一点是归纳中必要的；若只陈述“余维恰为
\(c\)”的版本，归纳不能应用于下面出现的 \(U_i\)。取满的仿射
partition（\(|I|=2^c\)、所有余维为 \(c\)）时，这正给出通常的
\(2^c-1\) 界。结论对不覆盖空间的部分 packing 同样成立。

## 两块的精确判别

对任意两个指标，
\[
A_i\cap A_j\ne\varnothing
 \quad\Longleftrightarrow\quad
\lambda_i|_{W_i\cap W_j}=\lambda_j|_{W_i\cap W_j}. \tag{1}
\]
这是线性方程组的相容性判别。因此块不交当且仅当两个限制泛函不同。
特别地，不交的两个块必有 \(W_i\cap W_j\ne0\)，但本证明**不**使用
任何关于三个法向空间公共交的断言。

## 证明

对 \(c\) 作归纳。\(c=0\) 时，所有法向空间为零，故总秩为零，
\((T_0)\) 成立。

设 \((T_{c-1})\) 已成立。若总法向秩为零，没有可证之处。否则从原族中
取一个按包含极小的子族 \(J\)，使得
\[
\sum_{i\in J}W_i=\sum_{i\in I}W_i. \tag{2}
\]
以下只考虑 \(J\)；它仍是两两不交族，且其法向和就是原来的法向和。
对每个 \(i\in J\)，令
\[
 E_i:=\sum_{j\in J\setminus\{i\}}W_j,
 \qquad U_i:=W_i\cap E_i,
 \qquad p_i:=\dim W_i-\dim U_i. \tag{3}
\]
极小性给出 \(W_i\not\subseteq E_i\)，所以 \(p_i\ge1\)。尤其
\(\dim U_i\le c-1\)。

把 \(\lambda_i|_{U_i}\) 作为法向条件，定义扩大的仿射平坦
\[
 B_i:=\{x\in X:\langle u,x\rangle=\lambda_i(u)
 \quad\text{对所有 }u\in U_i\}. \tag{4}
\]
它非空：任意 \(U_i\) 上的线性泛函都由 \(X\) 中某个点的内积表示。
并且 \(B_i\) 仍两两不交。事实上，对 \(i\ne j\)，有
\[
 W_i\cap W_j\subseteq U_i\cap U_j, \tag{5}
\]
因为 \(W_j\subseteq E_i\) 且 \(W_i\subseteq E_j\)。原块不交，结合
(1)，可取 \(w\in W_i\cap W_j\) 使
\(\lambda_i(w)\ne\lambda_j(w)\)。由 (5)，同一个 \(w\) 已是两条
扩张条件的矛盾，故 \(B_i\cap B_j=\varnothing\)。

于是按块大小相加得到 Kraft 型不等式
\[
 \sum_{i\in J}2^{-\dim U_i}\le1. \tag{6}
\]
另一方面，由 \(\dim W_i\le c\)，
\[
 \sum_{i\in J}2^{p_i}
 =\sum_{i\in J}2^{\dim W_i-\dim U_i}
 \le\sum_{i\in J}2^{c-\dim U_i}
 \le2^c. \tag{7}
\]
对每个正整数 \(p\)，有 \(p\le2^{p-1}\)。因此
\[
 \sum_{i\in J}p_i
 \le\frac12\sum_{i\in J}2^{p_i}
 \le2^{c-1}. \tag{8}
\]

法向空间 \(U_i\) 的维数至多 \(c-1\)，且 (4) 已证明它们也对应一个
两两不交仿射族。因此归纳假设给出
\[
 \dim\sum_{i\in J}U_i\le2^{c-1}-1. \tag{9}
\]
对每个 \(i\) 选取 \(U_i\) 在 \(W_i\) 中的补空间 \(P_i\)，则
\(\dim P_i=p_i\)，并且
\[
 \sum_{i\in J}W_i\subseteq
 \Bigl(\sum_{i\in J}U_i\Bigr)+\sum_{i\in J}P_i.
\]
将 (8)、(9) 代入，得到
\[
 \dim\sum_{i\in I}W_i
 =\dim\sum_{i\in J}W_i
 \le(2^{c-1}-1)+2^{c-1}=2^c-1.
\]
这完成归纳。\(\square\)

## 对余维三的直接后果

当每个 \(\dim W_i=3\) 时，任意数量的两两不交仿射块都满足
\[
\dim\sum_i W_i\le7.
\]
故 round 8 中剩余的五块、六块情形不必再分开分类，也不需要 one-hole
completion。更一般地，若某个新块不被其他法向空间之和包含，则把它的
法向条件放宽到 \(U_i\) 后，所得块至少有原来两倍测度；上述归纳是把这
个局部增益以 Kraft 不等式严格累加的方式。

## 审计点

1. (1) 使用的是完整的限制泛函相等判别；它不是“公共交上一点取值要
   两两不同”的错误替代。
2. (5) 是保留不交性的唯一关键：它将原来给出不相容的某个法向量保留在
   两个放宽后的条件中。
3. 极小子族只用于保持法向**总和**，并不声称被删块可由其余仿射块覆盖。
4. 该证明不依赖上一轮随机阴性搜索。递归证书的有限检验仅用于核查实现与边界条件，不替代一般归纳。

## 取等与来源范围

深度c的满二叉决策树，为每个内部节点配置一个全局独立二元查询变量。
每条根叶路径定义余维c的仿射块，所有叶块构成分割；全树共2^c−1个内部
变量都出现在某个叶块的法向空间中，故总法向秩恰为2^c−1。
这个例子也说明本定理的普适常数不能改善。

本轮文献核验尚未找到同一部分分割结论与上述证明；不据检索阴性宣称首次。
完整分割的块数减一秩定理与这里的最大余维界是不同表述。本文后续应用
明确使用这个自包含归纳，不删去任何不交性或余维条件。
