# 完成的d13数值批次

m56/n96，d13，seed2026111612，margin0.001，32768并行初值，12000步、每2000步重启，共196608次初值分配；全12291个合法过滤前缀在每个队列都覆盖。运行3898.606秒，exit0；父3844/子3854已不在进程表，句柄已回收。

本轮实际用独立Fraction验证器重新核验5个重启快照和末端共96行保存JSON小数；最小卷积二元残差0.002069334499537945，均非精确反例。独立覆盖器重新枚举3^12前缀，验证完整分配与最终步。它们不是有限排除。

在本目录运行 `python3 verify_prefix_fraction.py --run-dir run --worker worker/gpu_factor_search_prefix.py --prefix-module worker/fractional_prefixes.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d13-fraction-recheck.json`；覆盖使用 `audit_prefix_coverage.py` 同样run/worker并给record。

可再生state.pt未提交；全部必要保存样本与实际版本在run/config.json等文件。
