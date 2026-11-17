# R039：完成核传输给出的 $\mathbb E_1T$ 线性下界

## 严格结果

对有限连通简单图 $G$，记 $S(L)$ 为稳定提升源质量

$$
 S(L)=\sum_{v\ne u,\ v\not\sim u}\sum_{t=M_v}^{N-1}\frac1{t-1},
$$

其中 $u$ 是单峰词 $L$ 的根，$M_v$ 是 $v$ 的最高标号邻点。则

$$
 \boxed{\quad \mathbb E_1T\ge\mathbb E_1S.\quad} \tag{1}
$$

于是对自由边界 $P_n\square P_n$，由既有稳定提升的逐词下界，

$$
 \mathbb E_1T\ge B(N,4,N-1)
 =\bigl(1+3\log(3/4)+o(1)\bigr)N,\qquad N=n^2. \tag{2}
$$

这给出 $B_+/A_1\ge cN$ 的另一条完成核证明。它是严格的线性下界，但**没有**
给出题目当前需要的超线性发散因子，更不闭合原峰距问题。

## 完成数双计数

令 $a(S)$ 为以集合 $S$ 结束的合法一峰前缀数，$h(S)$ 为从连通集合 $S$ 开始的
合法完成词数。设 $|S|=k$，$w\in\partial S$，并写

$$
 c(S,w)=|N(w)\setminus N[S]|,
 \qquad M=N-k-1.
$$

每个 $v\in N(w)\setminus N[S]$ 正好对应 $T$ 的一个两步事件。因此

$$
 A_1\mathbb E_1T
 =\sum_{S,w} k\,a(S)
   \sum_{v\in N(w)\setminus N[S]}h(S\cup\{w,v\}). \tag{3}
$$

求和只需取 $1\le k\le N-2$ 的连通 $S$。另一方面，既有 first-touch 分解把
稳定提升源质量精确写成

$$
 A_1\mathbb E_1S
 =\sum_{S,w}a(S)h(S\cup\{w\})c(S,w)\,\omega_k, \tag{4}
$$

其中

$$
 \omega_k=\sum_{t=N-k}^{N-1}\frac1{t-1}
 =\sum_{j=0}^{k-1}\frac1{M+j}. \tag{5}
$$

式 (4) 是 `proofs/weighted-promotion.md` 的式 (5)；这里仅用它的精确计重，
不额外假设下一点在边界上均匀。

## 稳定前移的完成数传输

固定 $(S,w,v)$ 如上。任取 $S\cup\{w\}$ 的一个完成词，把 $v$ 稳定移到其剩余词
的第一位。因为 $v\sim w$，这是 $S\cup\{w,v\}$ 的一个合法完成词；移位只使
其他点的已见集合变大。反向时，固定一个 $S\cup\{w,v\}$ 完成词，至多能将 $v$
插回其余 $M-1$ 个点的 $M$ 个槽位。故严格有

$$
 h(S\cup\{w\})\le M\,h(S\cup\{w,v\}),
 \qquad h(S\cup\{w,v\})\ge\frac{h(S\cup\{w\})}{M}. \tag{6}
$$

代入 (3)，并对 $v$ 求和，得到

$$
 A_1\mathbb E_1T
 \ge\sum_{S,w}a(S)h(S\cup\{w\})c(S,w)\frac{k}{M}. \tag{7}
$$

而 (5) 中每一项都不超过 $1/M$，故

$$
 \omega_k\le\frac{k}{M}. \tag{8}
$$

由 (4)、(7)--(8) 即证 (1)。这是完成词上的稳定前移/有限纤维比较，不使用
Twin peaks Theorem 4.2，也不需要一个边界点的条件分布为均匀。

## 按 first-touch 秩的精确改进与为何仍停住

把 (4) 中秩为 $k$ 的非负总质量记为

$$
 W_k=\sum_{|S|=k,w}a(S)h(S\cup\{w\})c(S,w)\,\omega_k.
$$

则上述同一传输实际上保留了精确比较因子

$$
 A_1\mathbb E_1T\ge\sum_{k=1}^{N-2}\rho_{N,k}W_k,
 \qquad
 \rho_{N,k}:=\frac{k}{(N-k-1)\omega_k}. \tag{9}
$$

它严格优于只取 $\rho_{N,k}\ge1$。例如令 $M=N-k-1$，由

$$
 \frac1{M+j}=\frac1M-\frac{j}{M(M+j)}
 \le\frac1M-\frac{j}{M(M+k-1)}
$$

并对 $0\le j<k$ 求和，得到

$$
 \rho_{N,k}\ge
 \left(1-\frac{k-1}{2(N-2)}\right)^{-1}. \tag{10}
$$

这给出一条可用的、若能另行控制 $W_k$ 秩分布时的显式常数强化。

但它不自动产生发散因子。对于任意固定 $0<c<1$、$k/N\to c$，有

$$
 \rho_{N,k}\longrightarrow
 \frac{c/(1-c)}{\log(1/(1-c))}<\infty, \tag{11}
$$

而在 $k=o(N)$ 时 $\rho_{N,k}\to1$。现有静态度数/first-touch 下界只给总质量
$\sum_kW_k$ 的面积级控制，并不证明该质量必须进入 $k=N-o(N)$ 的末端秩区。因此
不能由 (9) 单独推出 $\log n$ 因子，更不能声称均匀一峰完成数律实际集中在任何指定
秩段。要继续，需要新的方格性输入：例如证明 $W_k$ 在逐尺度上有足够末端质量，或在
许多尺度给出严格优于 (6) 的平均逆像界。

这也解释了为何不能从 Theorem 4.2 的“大边界通常存在”直接升级：所缺的是
完成数诱导核在 active 两步对上的质量，而不是普通边界的基数。

## 有限核验（不承担一般证明）

`audit_T_ge_S_n3.py` 完整枚举 $P_3\square P_3$，得到

$$
 \sum T=79360,\qquad \sum S=\frac{262272}{5},
$$

故总和不等式严格通过。计数中有 $14640$ 个一峰词满足逐词 $S(L)>T(L)$，所以
(1) 必须按**平均/总和**理解，不能误写成逐词下界。
