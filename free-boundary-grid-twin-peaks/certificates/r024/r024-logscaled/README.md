# r024：以 log-normalizer 缩放的第二出生距离汇总

`aggregate_logscaled_distance.py` 是冻结的 r023 距离汇总器的独立
衍生版本。它不改动 Swift 终端核或原汇总器。在所有原始记录完成
SHA、metadata、种子、群组整数界、群组恒等式和精确距离阈值审计后，
它令

\[
 z_r'=\exp(\log Z_r-\max_s\log Z_s)
\]

并以 `z_r'` 计算全部比值、delta-SE 诊断和 ESS。公共缩放消去，显式地
限制实际指数的动态范围，并输出
`log_tau2_normalizer_mean`、`log_lift_full_normalizer_mean`，而不是把
极小 normalizer 强行转回 Double。

这仍是有限 SMC 的比值诊断；delta-SE 不是置信区间、混合证明或原题的
渐近结论。

## 实际审计

- 冻结 Swift source SHA-256：`01b77bca1735f07293eae5ffeea0326810fdc5ea785c7a9e99d9030cba9b0702`。
- 原 r023 聚合器 SHA-256：`cb466dcae1628f81ff2e198e65d2f95a71b198f77cd60fe8e717c9076403c358`。
- 新聚合器与真 CLI 夹具的 SHA 见 `test-logscaled-receipt.json`。
- `test_logscaled_sut.py` 实际调用新 CLI，验证手造两 replicate 权重
  1:3、整体 `logZ-1000` 且 Double 正常 underflow 后所有比值/SE 不变，
  并拒绝 NaN、非零 Z/logZ 不一致、伪 underflow、错种子和坏整数群组。
- 完成的本机 n=16 终端批与旧 Decimal-80 聚合器逐项一致至约
  \(10^{-79}\)；对应输出为 `n16-logscaled-summary.json`。

## n=32 范围审计

旧 `a022-random-tau-production-n32-r022` 的 39 条中断 random-\(\tau\)
记录并非本终端距离核：没有 metadata/run_end、没有第二出生 L-距离直方图，
故不能提交给本聚合器，也不能给出完整 \(K=2\) 距离估计。其实际
`logZ` 范围是 \([-1070.9996605068236,-918.5100414572485]\)，全体 `Z=0`；
Double 的最小 subnormal 对数约为 \(-744.4400719213812\)。这只说明
Double `Z` 列不能单独承载 normalizer；冻结 r023 聚合器已经使用
`Decimal(exp(logZ))` 重构权重，因此其比值计算也可处理这种真下溢。
新版本的附加作用是显式公共缩放和记录对数 normalizer，而不是修复旧聚合器的
下溢错误。

在冻结终端核输出完整 `metadata + R records + run_end` 后，先用本聚合器
和对应 sidecar/source 运行；再由主线审查其 source、种子域和原始记录。
