# 完成的d12 GPU数值批次

m52/n92、d12，32768并行轨迹，12000步、restart2000、margin.001、seed2026111610。全6605必要前缀，每个完整cohort分配4或5次；6cohort共196608次轨迹初始化。13136.349秒完成；原父子进程已退出，session已回收。
保存96行全部由独立Fraction验证器重算，最小精确十进制残差0.0020395991615655452，不是精确反例。GPU阴性结果不是有限排除。release记录退出前尚有108790528字节张量；退出释放由实际进程消失确认，不能把退出前值写成0。

```sh
/Users/zzh/.venvs/main/bin/python verify_prefix_fraction.py --run-dir run --worker worker/gpu_factor_search_prefix.py --prefix-module worker/fractional_prefixes.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d12-fraction-replay.json
/Users/zzh/.venvs/main/bin/python audit_prefix_coverage.py --run-dir run --worker worker/gpu_factor_search_prefix.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d12-coverage-replay.json
```

本轮已实际运行两验证器，结果分别在parent-fraction.json、parent-coverage.json。PyTorch2.15.0 MPS，NumPy2.5.4，Python3.13.18。两验证器不信任发现器成功标记。
