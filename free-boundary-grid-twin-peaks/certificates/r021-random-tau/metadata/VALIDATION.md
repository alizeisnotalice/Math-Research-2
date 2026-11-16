# 实际重放验证

归档的 `source/aggregate_random_tau.py` 实际读取、重放并严格验证了两份 raw
及其相邻 manifest。两次的 `replayed-summary.json` 与原有 summary 的逐字节 diff
为空；随后为保持最小包，重复的 replay 文件已删除。

- n=3 v4 raw：64 个 `globalrep`、64 个 distinct desired seed、6 个 canonical
  dead 行；负测试实际通过（错误 q、重复 seed、缺末行、伪 dead 被拒绝；
  `logZ=-1000` 的非死下溢行可汇总）。
- root n=32 audit raw：8 个 nondead 行都 `Z=0.0` 的 Double 下溢，仍从 `logZ`
  以 Decimal 80 得到与原 summary 相同的严格汇总输出。

独立 `postprocess_weight_diagnostics.py` 也实际重放所有 dispatch/seed/domain。
n=32 的 Z 与 D 权重 ESS 均显示为 1.0，最大归一化权重为 1.0（双精度显示精确
到 1）；这说明极端权重塌缩。summary 中的微小 delta-SE 不能解释为精度、概率、
渐近现象或反例证据。
