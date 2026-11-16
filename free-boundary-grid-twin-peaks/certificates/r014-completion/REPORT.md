# 提前界面根完成数增益：n=3、4 精确诊断

日期：2026-11-16（Asia/Shanghai）。本目录独立运行，未改动任何已有冻结包。

## 实现的量

对规范、非邻接、连通分支 (A,B)，令 (C=V\setminus(A\cup B))。对每一个外层
`B`，程序在诱导图 (G[V\setminus B]) 上重做完整 subset completion DP，精确得到

\[
H_A=h_{G[A\cup C]}(A),\qquad
F=\sum_{w\sim A,\,w\sim B}h_G(A\cup B\cup\{w\}),
\qquad \Gamma_A=F/H_A.
\]

反向的 (H_B,\Gamma_B) 独立重算。只要任一方向 (H>0) 就纳入可用并集；两方向均可用时，严格比较并选较大的

\[
D_*=a(A)a(B)H_A\left[\binom N{|B|}-\binom{|B|+|C|-1}{|B|}\right].
\]

程序按首次合并秩和小分支大小保存并集的精确 (W,D)。对每个规范 pair，以
(2^{128}) 的逐项向上取整累积 (W/D_*)，因此输出的是严格 dyadic 上界，不是浮点估计。

## n=3、4 结果

| n | 可用方向并集 W/A₂ | 并集内 Q | 并集 D/全 D | 严格区间 Σ W/D* |
|---:|---:|---:|---:|---:|
| 3 | 1 | 0.4895416110 | 1 | (50.2606165006074, 50.2606165006076) |
| 4 | 808981109/811769973 ≈ 0.9965644652 | 0.4230498795 | 0.9970575967 | (4771.204182463825, 4771.204182463827) |

n=4 中，按本轮补充的**正指数**阈值，\(\Gamma\le1\)、
\(\Gamma\le2^{|C|}\)、\(\Gamma\le4^{|C|}\) 的可用 pair 权重占全 (A_2) 分别为
0.78647843、0.99656447、0.99656447。此前的负指数诊断
\(\Gamma\le2^{-|C|}\)、\(\Gamma\le4^{-|C|}\) 权重均为零，现保留为单独标记，不能与正指数缩写混淆。各阈值采用“至少一个可用方向满足”的并集，未要求该方向同时最大化D*。

n=4 的最小 \(\Gamma\) 方向实例为 `A_mask=1,B_mask=1604`，选择规范 B 侧，
\(F=1612800,H=12418560\)，即 \(\Gamma=5/38\)；最大实例为
`A_mask=1,B_mask=1636`，选择规范 A 侧，\(F=322560,H=10\)，即 \(\Gamma=32256\)。
双方向均不可用的实例为 `A_mask=50,B_mask=19456`，其 (F=322560) 而两方向 (H=0)。
完整 n=4 中有 19,328 个规范 pair；其中 19,256 个可用、72 个双方向不可用，后者的
\(W=2677309440,D=15513262080\)，与全体 (A_2,D) 的差精确吻合。

每项用分母 (2^{128}) 向上取整；因此若有 `m` 个可用 pair，精确和严格位于
`(upper-m)/2^128` 与 `upper/2^128` 之间。上表给出了该窄区间，**下端本身远大于 1**；故在这份完整 n=4 pair 账本上，直接把单族上界求和确实产生严重熵损耗，不能给出有用的全局归一化结论。

## 独立核验与 n=5 评估

- n=3 实际枚举 (9!\) 个标号，直接恢复首次合并 pair，并对每个 pair 以独立递归计算两方向 completion 数；only-A、only-B、both 的 (W,D) 逐项等于 C++ 输出。
- n=4 的 11,506 个连通外层 `B` 产生 38,656 条有桥的有向记录和 19,328 个规范 pair，其中19,256个可用；JSON 中的各层和、类别和、全 (A₂,D) 相等及严格区间均由 `verify_completion_gain.py` 检查。
- n=5 仅做状态空间实测，未启动完整 DP：
  \(\sum_{B\text{ connected}}2^{25-|B|}=12,153,262,307\) 个 induced-DP 子集状态；即使额外要求 (V\setminus B) 连通仍有
  \(3,531,096,450\) 个。每个状态还要扫描边界/枚举可用 A，因此明显超出本轮 n=4 的合理扩大范围；没有据此把 n=5 当成已完成。

用于后续开发的容量校准（不是实跑结果）：n=4 的同一完整计算对应
2,168,528 个 `B` 与补图均连通的 induced-DP 状态，实际 wall 约 0.07 s；n=5 的
3,531,096,450 状态是其约 1628 倍。按该极乐观吞吐，单纯 induced-DP 已约 114 s；考虑
更大的 25 点边界扫描、可用 A 枚举与稀疏掩码映射，端到端先按数分钟到数十分钟估计，必须做前台检查点试跑。
当前全局五个 n=5 数组约 2.5 GiB；若仍让 12 个线程各持一个完整 `h_B` 数组，会另需约 6 GiB，
并和并行 GPU 作业争用 16 GiB 内存。后续应优先采用压缩 (H=V\setminus B) 子掩码数组，或将 CPU worker 限至 4，
再实际测峰值后决定是否扩大。

## 资源与复现

n=4 前台 12 线程实跑：wall 0.06 s、user 0.20 s、最大 RSS 24,559,616 bytes、0 swaps；全部进程已退出。

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
c++ -O3 -std=c++20 -pthread -I/opt/homebrew/include \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang \
  completion_gain.cpp -L/opt/homebrew/opt/gmp/lib -lgmpxx -lgmp \
  -o completion_gain
./completion_gain 3 completion-n3.json 12
./completion_gain 4 completion-n4.json 12
PYTHONDONTWRITEBYTECODE=1 python3 audit_completion_n3.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_completion_gain.py
c++ -O3 -std=c++20 estimate_n5.cpp -o estimate_n5
./estimate_n5 > n5-work-estimate.json
```

主Agent复核增补：重新编译n3/n4逐整数相同；验证脚本加入全A₂/D相等和严格区间下端>1检查；上表展示向外放宽的小数端点，精确dyadic端点以JSON为准。
