# d10 完成批次与独立样本审计

参数 m44,n84,d10，32768并行初值、12000步，每2000步重启一次，margin .001，seed2026111608。六个完整cohort共196608个轨迹初值，分别覆盖全部1671个必要前缀。实际完成3184.10秒，父27324/子27325均已退出；session22180 exit0。退出前release字段不代表进程退出后仍有常驻内存。

`run/`只保存复验所需JSON；41MiB可再生的state.pt留在CodexTemp。独立Fraction审计重算96个保存样本；这些保存小数都不是精确二元乘积。最低保存残差约0.0020387854892760515，不对未保存轨迹或一般问题给出排除结论。

本目录复验（Python标准库）：

```sh
python3 audit_d10_coverage.py --run-dir run --worker worker/gpu_factor_search_prefix.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d10-coverage-recheck.json
python3 verify_prefix_fraction.py --run-dir run --worker worker/gpu_factor_search_prefix.py --prefix-module worker/fractional_prefixes.py --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/d10-fraction-recheck.json
```

主线程已从本归档副本实际运行两个命令，均退出0；结果为 `parent-archived-coverage.json` 与 `parent-archived-fraction.json`。冻结GPU源码仅用于记录和核查，勿为复验重复启动同一批次。后继d11是独立的运行中批次。
