# H10 五根块 v2：第 25 轮完整 GPU 账本

保存批次 `H10-a3b4` 的固定完整账本、最小根输入、MP80 独立验证收据及复现脚本。options 缓存（822 KiB）是可再生产物，未归档。此归档仅证明保存账本、输入和数值重构相互一致；不提供严格根证书，也不构成一般猜想的排除证明。

冻结 checkpoint `status=complete`，覆盖全部 9 个 (E) 键和全部 `2,397,558` 个规范布局代码。独立验证从高精度 `lift_root` 构造每个实二次因子，未钉端点；重构 528 个 top-64（末端 16 个）保存样本，最大卷积残差为 `2.61965185130923917706482586315e-76`，最大端点偏差为 `8.43375835458441857947885924102e-81`。

通用计数为 `F_5(3) F_5(4)/2 + 2^(5-1)-2 = 2,397,558`。运行：

```sh
python3 -m pip install --target /Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps mpmath==1.3.0 numpy==2.5.5
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps \
python3 computation/results/025-H10-coprime-gpu/scripts/verify_coprime_v2_generic.py \
  --checkpoint computation/results/025-H10-coprime-gpu/checkpoint/checkpoint.json \
  --options /path/to/H10-a3b4-options.json \
  --roots-input computation/results/025-H10-coprime-gpu/input/input.json \
  --engine-source computation/results/025-H10-coprime-gpu/scripts/stream_coprime_ds_gpu_v2.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-MP80-receipt.json
```

实际再生过的 options 命令为：

```sh
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps \
python3 computation/results/025-H10-coprime-gpu/scripts/prepare_coprime_inputs.py \
  --a 3 --b 4 --mother '1,1,1,0,0,0,0,0,0,1,1' --mother-name H10 \
  --mother-certificate /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/p8_next_input/Hd10_a3b4/H10-mother-certificate.json \
  --out-input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-rebuilt-input.json \
  --out-options /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H10-rebuilt-options.json
```

该实际命令生成 input/options SHA-256 分别为 `18660e457af345c0fe0850cc4b7484bfd137e9d5399cb0162cf7e67e2a50acb6`、`3bc5f3823b6394bb8455dc582d22b6beb856a4034207d2d127e14559d0aceebb`，与 checkpoint 一致。母证书的路径字符串进入输入 JSON，需保持上述路径以复现字节哈希；归档副本仅用于保存来源内容。
