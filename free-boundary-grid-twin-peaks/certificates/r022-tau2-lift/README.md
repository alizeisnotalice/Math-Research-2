# r022 τ=2 lift 严格汇总器

对 τ=2 的 SMC 输出，令每个粒子的 (L) 为低峰首邻点的秩减二。严格插槽恒等式给出
\(Z_{K=2}=Z_{\tau=2}E(L)\)，以及 (D,F,G) 的相应 (L)-加权式。因此脚本将
原 τ=2 normalizer `Z` 与 lift 后 full-​K2 normalizer `Z*mean(L)` 明确分开，且只报告
ratio-of-sums：
\[
 Q=\frac{\sum_r z_r d_r}{N\sum_r z_r l_r},\quad
 V=\frac{\sum_r z_r f_r}{8\sum_r z_r l_r},\quad
 C=\frac{\sum_r z_r g_r}{4\sum_r z_r l_r}.
\]

`aggregate_tau2_lift.py` 强制检查 metadata、run_end、完整 replicate 集合、GPU seed、sidecar
的 source/raw SHA 及其期望 `n/P/R/tau/domain/method/move_steps`、每组整数和、容量/溢出界、每组
`F-2G=4D`、scaled/complete 字段、死行规范性和 Double 下溢。
它以 Decimal 80 从 `logZ` 恢复贡献。其 delta 标准误依据给定的有限重复 influence
式；不以 ESS 代替该公式。ESS/maxshare 仅为权重集中度诊断。

自由 (P_n\square P_n)（(n\ge2)）的 τ=2 批不允许 `dead=true`：第一根总有
非邻点可选作第二根，随后任一未满前缀都有非空外边界。因此该状态是实现不变量
失败，而非本目标可接受的零贡献；`Z=0` 仅在有限 `logZ` 真正发生 Double 下溢时允许。

实际小验证（非生产）已读取 n=3, P=17, R=2 raw 与 sidecar；结果在
`n3-p17-r2-strict-summary.json`。这既不认证有限误差，也不证明原条件律的渐近结论。

复现：

```sh
python3 -B aggregate_tau2_lift.py RAW.jsonl --sidecar RECEIPT.json --source tau2_lift.swift
python3 -B test_tau2_lift_aggregate.py
```
