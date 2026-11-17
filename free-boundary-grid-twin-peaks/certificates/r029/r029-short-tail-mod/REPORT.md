# r029：穿孔方格短尾 (h_1,h_2,h_3) 的模数滚动容斥

状态：有限精确实现与小规模交叉通过；不构成 (n\to\infty) 的结论。n=7 在程序中
需要显式 `--enable-n7`，本轮没有启动。

## 计算的候选恒等式

固定最高点 (v)，令 (H=G-v)、(M=|H|=n^2-1)、(N_v=N_G(v))、
(Z=V(H)\setminus N_v)。对每个 (I\subset Z) 的独立集，程序滚动计算
(C_H(I))，其中 `C_H(empty)=M!` 且

\[
C_H(I)=\frac{\sum_{u\in I}C_H(I\setminus\{u\})}{|N_G[I]|}.
\]

令 (U=V(H)\setminus N_G[I])、(P=Z\cap U)，并以 (e(A,U)) 表示从 (A)
到 (U) 的有向原图边数。实现累加下列候选分子并最后在有限域中除以
((M)_2,(M)_3,(M)_4)：

\[
\begin{aligned}
h_1&=\sum_I(-1)^{|I|}C_H(I)e(P,U)/(M)_2,\\
h_2&=\sum_I(-1)^{|I|}C_H(I)\sum_{q\in N_v\cap U}e(P\setminus N[q],U)/(M)_3,\\
h_3&=\sum_I(-1)^{|I|}C_H(I)\sum_{q\in N_v\cap U}\sum_{r\in(N_v\cup N(q))\cap U\setminus\{q\}}
e(P\setminus(N[q]\cup N[r]),U)/(M)_4.
\end{aligned}
\]

这些公式的独立证明不由本 C++ 实现承担。尤其不能声称每一个
`C_H(I)/(M)_(t+1)` 均为整数：当该 (I) 没有有效 prefix 时前缀边权为零。
本程序正确地先在模域累加完整分子、最后乘 falling-factorial 的逆元，因而不需要
这种错误的逐项整除假设。

`e(A,U)` 以四个方向的移位掩码分别 `popcount`，保留有向边重数；不枚举每一个
(p\in P)。

## 有效小规模核验

四个 r027 Pratt 认证素数的乘积是 252 bits，大于 (48!\) 的 203-bit 上界；故
n≤7 的非负整数尾箱可唯一 CRT 恢复。

| 范围 | 独立参照 | 结果 |
| --- | --- | --- |
| n=3，所有 9 个孔 | r028 exact prefix-tail DP | 所有 (h_1,h_2,h_3) 匹配 |
| n=4，所有 16 个孔 | r028 exact prefix-tail DP | 所有 (h_1,h_2,h_3) 匹配 |
| n=5，中心孔 v=12 | r029 shared-suffix/full prefix DP | 三项匹配：33,603,716,782,815,379,200；22,947,944,884,729,956,000；16,498,461,180,178,692,000 |

本目录的 C++ 四模重构结果与 r028 表的直接逐项比较在
`candidate-formula-crosscheck.json`；n3/n4 另由
`r029-short-tail-ie/short_tail_ie_check.py` 独立重跑候选
IE 到 r028 tail DP 的比较，输出 `independent-ie-vs-r028.json`。n5 的独立 full-DP
来源为 `tail-full-n5-center.json`，其 SHA-256 记录在前一 JSON。

开发中首次小批 `runs/candidate-small/` 发现了一个真实 bug：把四方向邻居 mask 取并集会
漏数两个不同 (p) 指向同一终点的有向边。该批 output 已保留但**明确无效**；修正后源码用
每个方向单独 popcount，所有认可 checkpoint 位于 `runs/candidate-small-v2/`。

## n=7 资源判断

中心孔的滚动数组/哈希 payload 模型仍为 2,484,019,712 B（2.3134 GiB），但短尾贡献每个
独立集额外最多需 1 个 (h_1)、4 个 (h_2)、28 个 (h_3) 的有向边计数。内存模型几乎
不变，wall time 却不能从 r028 的 85.03 s 可信外推。详见
`n7-cost-assessment.json`。本轮没有运行 n=7，亦不建议为一个已受解析质疑的全局尾收缩断言
盲目扩大有限计算。

## 复现

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread short_tail_mod.cpp -o short_tail_mod
./short_tail_mod 5 12 0 1 runs/candidate-small-v2/n5-v12-p0.json
python3 verify_candidate_formula.py
```

源码、二进制、核验结果与边界测试均列入 `SHA256SUMS`。所有小型子进程均正常退出。
