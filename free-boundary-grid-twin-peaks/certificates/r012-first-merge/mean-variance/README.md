# 首次合并条件根均值—方差分解

日期：2026-11-16（Asia/Shanghai）。本目录独立于此前冻结的空间矩和根深度包。

## 严格固定点区间

对每个首次合并分支 (A,B)，保存精确整数

\[
a_H,\quad M_x(H)=\sum r_x,\quad M_y(H)=\sum r_y,
\quad M_2(H)=\sum\|r\|_2^2.
\]

它们都遵从单峰删除递推。令 (L=2^{18})，每个均值坐标的严格夹逼为

\[
\frac{\lfloor LM_x/a\rfloor}{L}\le \mu_x
\le\frac{\lceil LM_x/a\rceil}{L},
\]

且 (y) 坐标同理。坐标非负，故四个坐标平方相加即给

\(W(\|\mu_A\|^2+\|\mu_B\|^2)\) 的有理上下界，公共分母为
\(L^2\)。精确的根二阶矩总和是

\[
R_2=h\binom{|S|-1}{|A|}\,[a(B)M_2(A)+a(A)M_2(B)].
\]

于是逐状态严格地取

\[
\mathrm{within}=R_2-W(\|\mu_A\|^2+\|\mu_B\|^2),\qquad
\mathrm{between}=D-\mathrm{within}.
\]

定点向上舍入可能使一个状态的 `within` 下界略负，程序仅在下界端用已知
`within >= 0` 截断为零；因此仍是保守的严格区间。每个
`(首次合并秩 s, min(|A|,|B|))` 的完整整数账本保存在 JSON。

宽度认证：n=5 中所有总和均由

\[
64L^2\,25!<2^{128}
\]

控制；程序启动时用除法形式实际检查该不等式，所有乘积以 `__uint128_t` 完成。

## 结果

下表的量均已除以 (n^2A_2(n))；表中小数仅为区间端点的近似展示，严格端点以 JSON 整数除以 L²n²A₂ 为准。

| n | (Q_n) | 组内方差项 | 组间均值差项 | 组内占 Q |
|---:|---:|---:|---:|---:|
| 3 | 0.4895416110 | [0.0536789524, 0.0536791097] | [0.4358625013, 0.4358626586] | 约 10.96% |
| 4 | 0.4228406446 | [0.0554534352, 0.0554540968] | [0.3673865479, 0.3673872094] | 约 13.11% |
| 5 | 0.3830041101 | [0.0512052524, 0.0512062813] | [0.3317978288, 0.3317988577] | 约 13.37% |

n=5 的 110 个非零 `(s,minbranch)` 层可直接检查任意秩/平衡度子族；全局
`W=A2=7063261481362958280000` 和
`D=67631454448823894520000` 与两份独立精确证书一致。

有限诊断显示：在这三个尺寸中，**组间条件根均值差是 (Q_n) 的主要部分，组内方差约为 11%–13%**。这只描述有限精确账本，不能外推为均值项、方差项或 (Q_n) 的渐近结论。

## 核验与资源

- n=3 由独立 Python `Fraction` 子集递推复算精确均值、组内、组间总量，并验证均落在 C++ 固定点区间内，且 `within + between = D` 精确成立。
- n=3/n=4/n=5 都重新恢复既有精确 `A2,D`；分层 `W,D` 和逐层上下界可由输出 JSON 加总。
- n=5 前台 12 线程实际运行：wall 2.54 s、user 12.24 s、sys 0.66 s、最大 RSS 2,754,985,984 bytes、0 swaps；结束后无残留进程。

复现：

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
c++ -O3 -std=c++20 -pthread \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang \
  first_merge_mean_variance.cpp -o first_merge_mean_variance
./first_merge_mean_variance 3 mean-variance-n3.json 12
./first_merge_mean_variance 4 mean-variance-n4.json 12
./first_merge_mean_variance 5 mean-variance-n5.json 12
PYTHONDONTWRITEBYTECODE=1 python3 audit_mean_variance_n3.py
```
