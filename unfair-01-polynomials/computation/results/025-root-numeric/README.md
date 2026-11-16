# 第 25 轮最终归档候选包

本包将两份已独立 MP80 核验的数值归档与一份新增的精确整数账本 wrapper 合并，供主线程择白名单同步。包仅保存最小可复现实物；GPU options 缓存不在包中，数值结论与严格证书结论严格分列。

## 最小白名单

- `numeric/archive025-H8/**`：H8 的固定完整 checkpoint、输入、独立 MP80 receipt、v2 engine、再生及布局验证脚本和母证书副本。
- `numeric/archive025-H10/**`：同上，H10 五根块版本。
- `scripts/audit_integer_coverage_singular.py`：不导入搜索器的单 checkpoint 整数覆盖/mask 审计器。
- `integer/H8-integer-coverage-receipt.json`、`integer/H10-integer-coverage-receipt.json`：本轮实际重跑的精确组合账本收据。
- `scripts/verify_modular_galois.py`、`galois/H8-*`、`galois/H10-*`：本轮实际重跑的 H8/H10 模素数 Galois 检查及其输入。
- `SHA256SUMS`、`archive-map.json`、本 README：完整性与来源映射。

不在白名单中的 options JSON 是可再生产物：使用各 numeric 子归档的 README 命令生成，并以 checkpoint 内的 `options_sha256` 拒绝不匹配缓存。

## 已实际核验

H8：整数 wrapper 验证 4 块、1,482/1,482 全覆盖及 216 保存 mask；H10：5 块、2,397,558/2,397,558 全覆盖及 528 保存 mask。wrapper 检查 `C=H(x^a)H(x^b)` 的整数恒等式、每个 E、mixed-radix 码、bitmask、次数和 checkpoint/source/input SHA。

H8 的模素数检查验证了素数 2、13、5 上的不可约/平方自由分解及相应 transitive、long-cycle、transposition 角色；H10 实际重跑验证了素数 7、73、19 上的对应角色。它们是有限域 Galois 证书审计，不能替代 H8 exact interval/allocation 证书；后者由主线程的严格包另行处理。

数值 MP80 receipt 只验证数值根输入的保存样本与完整 GPU ledger 一致，**不**是严格根、符号或一般排除证明。

## 复现命令

```sh
python3 scripts/audit_integer_coverage_singular.py \
  --checkpoint numeric/archive025-H8/checkpoint/checkpoint.json \
  --input numeric/archive025-H8/input/input.json \
  --engine-source numeric/archive025-H8/scripts/stream_coprime_ds_gpu_v2.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H8-integer-recheck.json

python3 scripts/verify_modular_galois.py galois/H8-mother-certificate.json
```
