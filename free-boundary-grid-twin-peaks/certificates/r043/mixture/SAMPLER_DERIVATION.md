# 固定曲线事件的精确完成计数与混合抽样

本文件使用而不重新发明
`proofs/completion-kernel-counterexample.md:43` 的完成数容斥。对诱导区域
\(H\) 的非空连通下降标号前缀 \(S\)，令
\(R=H\setminus S\)、\(F=R\setminus N_H(S)\)。余下词合法当且仅当
\(F\) 中没有顶点在其于 \(R\) 的闭邻域中最高。故按独立峰集容斥，

\[
h_H(S)=\sum_{I\subseteq F,\ I\text{ independent}}(-1)^{|I|}C_R(I),
\qquad C_R(\varnothing)=|R|!,
\]

并且

\[
C_R(I)=\frac{\sum_{v\in I}C_R(I\setminus\{v\})}
 {|I\cup N_R(I)|}.
\]

分母整除由闭邻域最大标签分解保证，程序逐项断言。若首根为允许且全图安全的
\(r\)，则有 \(h_H(\{r\})\) 个完成词，故以
\(h_H(\{r\})/a_H(A)\) 选根；给定非终止前缀后，以
\(h_H(S\cup\{v\})/h_H(S)\) 选相邻候选 \(v\)。候选完成数分割
\(h_H(S)\)，所以在理想独立均匀整数源下每个区域词严格等概率。

左右区域词独立均匀，且在 \(\binom{n^2}{|L_h|}\) 个交错中均匀选一个，
便得到固定 \(E_h\) 上的严格均匀词。无偏 `randbelow` 由两个均匀 64 位字拼成
128 位字并作拒绝采样；实现中的确定性 PRNG 因而只产生数值样本，并不把有限
伪随机样本称为严格概率估计。

若以后以精确 \(C_h/S_1\) 选曲线，混合律为
\(Q(L)=X(L)/S_1\)，其中 \(X\) 是命中的曲线数。因此有限恒等式为

\[
\mathbb P_2(X>0)=\frac{S_1}{A_2}\mathbb E_Q[1/X],
\qquad \frac{S_2}{S_1}=\mathbb E_QX.
\]

只有这些等式和整数输入是严格的；有限 PRNG 样本是数值证据。每个
`CompletionCounter` 的前缀缓存和独立集临时表各按 `--cache-mib` 给出保守的
条目预算，默认及硬上限为 128 MiB，并按曲线分别处理。该参数不是整个 worker
RSS 的硬上限：左右两个 counter 及其暂时重叠的表仍会相加，实际 RSS 必须用
`time -l` 收据测量。
