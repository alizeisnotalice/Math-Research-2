# r002：强制第二出生秩的 CPU 方差诊断

这是对原条件律的独立数值核验，不是原问题的概率认证或证明。

## 方案

在降序揭示中固定第二次 birth 的秩 `tau`：秩 1 任取；`1<t<tau` 只从当前 frontier 取；`t=tau` 只从尚未揭示的非-frontier 顶点取；其后继续只从 frontier 取。无可选顶点的路径作为零贡献 dead end。每个轨迹乘以所有被选择集合大小之积 `W`。因此，逐层样本均值估计相应的计数层；层和对应 (A_2(n))，而 (W f) 中 (f=\|X-Y\|_2^2/n^2) 的层和比值给出 (Q_n) 的 SIS 比值估计。

每层保存 `sumW, sumW2, sumWf, sumW2f, sumW2f2`。`variance_diagnostics.py` 从这些原始矩计算：全局 ESS、单样本最大贡献份额、按层质量份额、以及比值估计的 delta-method 描述性标准误。它要求调用者显式白名单完整域，并验证完整的 tau 集、每层 batch 数、样本数和唯一 seed；不会扫描或采用 `invalid-json` 目录。

## 数值设置与实际运行

- C++20，`mt19937_64`，rejection-uniform 索引；每个运行 12 threads。
- 种子：`splitmix(0xA022C001 + domain*0x9e3779b97f4a7c15 + n*1000003 + tau*1009 + batch)`。
- `n=3`：2 batches × 100,000 / tau，domain 0。
- `n=7`：4 个 domain；各为 4 × 100,000 / tau，即每域 19,200,000 条轨迹。
- `n=8`：8 个 domain；各为 4 × 200,000 / tau，即每域 50,400,000 条轨迹；总计 403,200,000 条。
- 外部独立峰判定器复核：n=3 12/12，n=8 八域各 12/12 保存标号均通过。
- 本机实测 `numeric_limits<long double>`：`digits10=15`, `max_exponent10=308`。本直接 W/W² 实现只允许 `n<=8`；n=12 必须转为 log-sum-exp/缩放实现。

## 结果（有限样本，描述性）

- n=3：(A_2\) 估计 143709.54822（精确 143808，比例 0.9993154）；(Q=0.48985816)，精确 (1100/2247=0.48954161)，差 (3.17\times10^{-4})，该域 delta SE (2.84\times10^{-4})。
- n=7：四域的简单均值 (Q=0.34358490)，域间 SD 0.00115921，域均值的描述性 SE 0.00057960；pooled 原始矩比值 0.34357585。精确值 0.34441831。
- n=8：八域的简单均值 (Q=0.33731306)，域间 SD 0.00517337，域均值的描述性 SE 0.00182906；pooled 原始矩比值 0.33734083。两个量刻意分别报告，不能混同。
- 单个 n=8 域的全局 ESS 为 591.95--1292.78，最大单样本全局贡献份额为 0.01265--0.02424；90% 质量均落在 tau=2,...,8。它们表明权重退化仍显著；绝不能由最小层 ESS 单独推出全局退化，亦不能据此构成严格置信区间。
- 与另一 GPU 探索值 (Q_8\approx0.3379) 的差，和本 CPU 域间波动量级相近；GPU 自身误差未在此包估计，所以该一致性只是一项交叉数值检查，不是数学结论。

## 复现

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang c++ -O3 -std=c++20 -pthread variance_sis.cpp -o variance_sis
./variance_sis 8 4 200000 12 n8-domain0 0
python3 variance_diagnostics.py . n8-domain0
python3 verify_saved_labels.py n8-domain0/saved-k2-labels.json
```

已完成的进程均已退出；无常驻 SIS 进程。

归档注：主 Agent 已逐项核验原交接 SHA256SUMS。当前包省去可再生探针二进制、重复诊断输出与旧校验清单；SOURCE-RUN-HASHES.txt 只保存原交接路径/哈希。当前有效载荷由课题 manifest-r011.json 记录。各域使用不同种子，不把这一点称为有限伪随机流独立性的证明。
