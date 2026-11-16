# 主稿双向接口常数 60 的独立核验

主稿的双向接口式给每次矩形合并

\[
\left|g(\mathrm{union})-g(\mathrm{left})-g(\mathrm{right})\right|
\le \sum_{U\in\{L,R\}}(\mathrm{width}(U)+\mathrm{height}(U))\log|U|
+\log4.
\tag{1}
\]

(1) 用于本文的平衡递归合并：接口两侧宽度之比在 [1/2,2]，故下界交错代价至多 log 3，上界代价为 log 4。这里没有要求各矩形自身纵横比受控；只要两边均不超过一个尺度 \(L\)，每个面积对数不超过
\(2\log L\)。

对一个父阵列块，设实际宽、高为 \(w,h\le L\)，将它分为四块、先合并上行和下行、再合并
这两行。三个合并的六个移根指数之和恰为

\[
(w+2h_1)+(w+2h_2)+(2w+h)=4w+3h.
\]

故所有移根对数项至多 \(14L\log L\)；三项交错/反向常数各至多 \(\log4\)。当存在内部
节点时 \(L\ge2\)，所以

\[
14L\log L+3\log4\le15L\log(2L).
\tag{2}
\]

令 \(J=\lceil\log_2q\rceil-1\)。若 \(q=1\)，没有内部节点；若 \(q>1\)，
内部层正是 \(0\le\ell\le J\)，有

\[
2^J<q\le2^{J+1},\qquad L_\ell\le\frac{2qk}{2^\ell},\qquad
\#\mathrm{blocks}_\ell\le4^\ell.
\]

于是总误差至多

\[
30qk\sum_{\ell=0}^J2^\ell\log\frac{4qk}{2^\ell}.
\]

令 \(t=J-\ell\)，并用 \(2^J\le q\)、\(q/2^J\le2\)，该量不超过

\[
30q^2k\sum_{t\ge0}2^{-t}[\log(8k)+t\log2]
\le60q^2k\log(16k).
\tag{3}
\]

故主稿的

\[
|g(qk,qk)-q^2g(k,k)|\le60q^2k\log(16k)
\]

严格成立。令 \(q\to\infty\) 后，已证极限 \(f=\lim g(n,n)/n^2\) 满足

\[
|g(k,k)-fk^2|\le60k\log(16k),
\]

即主稿 (9a) 的速率也成立。

`audit_60_levels.py` 实际以整数/Fraction 检查 \(1\le q\le4096\) 的正确内部层范围、
\(q=1\) 的空内部节点、尺度上界和两条 dyadic 几何级数系数界，共检查 45057 个内部层。
这只是账本的有限交叉核验；(1)--(3) 是一般证明。

## 矩形版本的独立复核

对 a≤b，将 qab 方格铺成 qb×qa 个 a×b 叶块，叶项为 q²ab·g(a,b)。置 Q=qb。
第 ℓ 层两条物理边分别不超过 a⌈Q/2^ℓ⌉、b⌈qa/2^ℓ⌉，均不超过 b⌈Q/2^ℓ⌉。
于是同一账本给 60Q²b log(16b)=60q²b³ log(16b)。
除以总面积 q²a²b²，令 q→∞，即得主稿 (12) 的 60b/a² log(16b)。
某方向块数为 1 时只会减少合并次数；不要求中间矩形纵横比有界。
空内部层的条件是 Q=1，而非原缩放参数 q=1。本节由独立 Agent 逐式复核，未发现漏因子。
