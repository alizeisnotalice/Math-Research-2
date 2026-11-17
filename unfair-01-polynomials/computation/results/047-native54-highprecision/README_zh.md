# degree 54：task 18--23 near-unknown 的独立高精数值诊断

状态：**数值证据，非严格排除或根分配证书。** 本目录只处理父方
`candidate-diagnostics.json` 中 task 18--23 标为 `nearunknown` 的 235 个 mask；没有重新枚举。

## 方法

输入 `C` 的系数是整数。输入中的根组只是 binary64 种子；脚本以 80 decimal digits 的
`mpmath` Newton 法对每个种子求解精确整数多项式 \(C(z)=0\)，并要求更新与残差不超过
\(10^{-55}\)。它独立重写了冻结的 group-0 固定到 \(B\)、其余组按大端 bit 读取的 mask
约定。每个共轭组随后按 \(\prod(1-x/z)\) 重构，且直接重算 \(AB-C\)。

根来自数值种子，故即使残差很小也不构成“该 mask 对应精确代数因子”的证明。若系数能舍入为
0/1，脚本会额外以 Python 整数卷积检查该两个舍入多项式的乘积是否**恰为** \(C\)；这会是独立的
公平分解证据。本批没有这种命中。

## 实际运行

解释器：`/Users/zzh/.venvs/main/bin/python3`；`mpmath 1.3.0`；80 dps；Newton 最多 80 步；
判定门槛 \(10^{-55}\)。命令为：

```sh
/Users/zzh/.venvs/main/bin/python3 \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round047-native54-hi/high_precision_newton_verify.py \
  --input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/gpu_partial39/degree54_native45/inputs-degree54.json \
  --diagnostics /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round047-parent-native54/candidate-diagnostics.json \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round047-native54-hi/high_precision_task18_23.json
```

实际结果：235/235 均为 `negative`；没有 `numerically_nonnegative_nonbinary`、
`numerically_fair` 或根细化失败。六个目标的最小根间距都至少
\(6.0673888160\times10^{-2}\)，所有组联合重构 \(C\) 的误差至多约
\(2.75\times10^{-69}\)，全部候选重构乘积误差至多约 \(3.02\times10^{-69}\)。

最接近零的负系数仍远离数值边界：task 18 至 23 分别不大于
\(-0.5070,-0.6271,-0.7059,-0.8003,-0.5524,-0.7098\)。这些值说明本批的原始
binary64 `nearunknown` 标签是低精度重构残差造成的，不是非负候选证据；但这不推广到
未跑 mask，也不代替代数数认证。
