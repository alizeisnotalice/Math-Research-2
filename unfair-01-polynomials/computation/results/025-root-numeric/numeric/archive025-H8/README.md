# H8 四根块 v2：第 25 轮完整 GPU 账本

保存批次 `H8-a2b3` 的固定完整账本、最小根输入、MP80 独立验证收据及复现脚本。options 缓存（130 KiB）是可再生产物，未归档。此归档仅证明保存账本、输入和数值重构相互一致；不提供严格根证书，也不构成一般猜想的排除证明。

冻结 checkpoint `status=complete`，覆盖全部 5 个 (E) 键和全部 `1,482` 个规范布局代码。独立验证从高精度 `lift_root` 构造每个实二次因子，未钉端点；重构 216 个保存样本，最大卷积残差为 `1.02689441725419880623734590119e-76`，最大端点偏差为 `1.05421979432305232243485740513e-80`。

通用计数为 `F_4(2) F_4(3)/2 + 2^(4-1)-2 = 1,482`。运行：

```sh
python3 -m pip install --target /Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps mpmath==1.3.0 numpy==2.5.5
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps \
python3 computation/results/025-H8-coprime-gpu/scripts/verify_coprime_v2_generic.py \
  --checkpoint computation/results/025-H8-coprime-gpu/checkpoint/checkpoint.json \
  --options /path/to/H8-a2b3-options.json \
  --roots-input computation/results/025-H8-coprime-gpu/input/input.json \
  --engine-source computation/results/025-H8-coprime-gpu/scripts/stream_coprime_ds_gpu_v2.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H8-MP80-receipt.json
```

实际再生过的 options 命令为：

```sh
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/mp80-deps \
python3 computation/results/025-H8-coprime-gpu/scripts/prepare_coprime_inputs.py \
  --a 2 --b 3 --mother '1,1,1,0,0,0,0,1,1' --mother-name H8 \
  --mother-certificate /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/p8_next_input/Hd8_a2b3/H8-mother-certificate.json \
  --out-input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H8-rebuilt-input.json \
  --out-options /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H8-rebuilt-options.json
```

该实际命令生成 input/options SHA-256 分别为 `74d1287798fb12764921d72062f6ab30682ea33567682d9be61841854de43049`、`47673274aa54c49df4f6bd5511dc3e3cf8466d46304114572f83fbd3c49b7988`，与 checkpoint 一致。母证书的路径字符串进入输入 JSON，需保持上述路径以复现字节哈希；归档副本仅用于保存来源内容。
