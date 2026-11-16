# R8 四根块 v2：第 24 轮冻结 GPU 前缀归档（rev2）

归档最终应置于课题目录 `unfair-01-polynomials/computation/results/024-R8-coprime-gpu`。
它保存四个由主线程停止并回收的 GPU 扫描前缀、对应最小根输入、独立 CPU 验证收据和可复现脚本；不保存巨大的 options 缓存。`archive024` 原版保持不动；本 rev2 的唯一实质修订是 MP80 验证器从 90 位 `lift_root` 字符串重新构造每个二次因子，而不再把 float64 系数副本当成高精度输入。

每个 checkpoint 的 `status` 都是 `running`，因此本归档的结论严格限于“保存的连续前缀和 top-64 行一致”。它不排除未扫描尾部，不构成一般反例排除或原猜想的证明。

| 批次 | `(a,b)` | 冻结前缀代码数 | MP80 重构行数 | 缺少原 CPU copy 的保存行 |
|---|---:|---:|---:|---:|
| a2b9 | (2,9) | 547,707,841 | 336 | 64 |
| a2b11 | (2,11) | 1,628,286,147 | 256 | 64 |
| a3b10 | (3,10) | 1,451,673,970 | 320 | 64 |
| a3b11 | (3,11) | 974,892,308 | 320 | 64 |

四根块的完整规范布局计数是
`F_4(a) F_4(b) / 2 + 2^(4-1) - 2`；最后一项明确是终端两种模式的 `+8-2` 修正。独立验证器以 `mpmath` 1.3.0、80 位十进制精度重构每个保存行，不导入 GPU 搜索器。

每一个高精度二次因子都由 `lift_root=[re,im]` 重建为
`[re^2+im^2,-2re,1]`。验证器额外要求该结果与原 float64 `coefficients_ascending` 副本逐项相差不超过 `1e-12`，不把端点钉为 1，并要求实际归一化端点偏差不超过审计容差 `1e-7`。这仍是数值根输入审计，而不是精确代数根证书。

## 复现独立审计

先在临时目录安装纯 Python 依赖（不改系统环境）：

```sh
python3 -m pip install --target /Volumes/data/CodexTemp/T2/unfair-01-polynomials/r8-mp80-deps mpmath==1.3.0
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/r8-mp80-deps \
  python3 computation/results/024-R8-coprime-gpu/scripts/test_r8_v2_layout.py
```

对冻结前缀 `a2b9` 的命令（其他三个只替换标签）：

```sh
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/r8-mp80-deps \
python3 computation/results/024-R8-coprime-gpu/scripts/verify_r8_coprime_v2.py \
  --checkpoint computation/results/024-R8-coprime-gpu/checkpoints/a2b9-checkpoint.json \
  --options /path/to/regenerated/a2b9-options.json \
  --roots-input computation/results/024-R8-coprime-gpu/inputs/a2b9-input.json \
  --engine-source computation/results/024-R8-coprime-gpu/scripts/stream_coprime_ds_gpu_v2.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/a2b9-MP80-receipt.json \
  --allow-partial
```

`--allow-partial` 是必需的保护门：没有它，验证器会拒绝把中断前缀记作完整审计。

## 重建 options 缓存

每个 root input 的 `parent_roots_source` 指向先前保存的 R8 数值根输入；准备脚本依赖该输入取得同一组父根。数值根输入与精确Galois证书用途不同。示例：

```sh
PYTHONPATH=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/r8-mp80-deps \
python3 computation/results/024-R8-coprime-gpu/scripts/prepare_coprime_inputs.py \
  --a 2 --b 9 --mother '1,1,1,1,1,1,0,0,1' --mother-name R8 \
  --parent-roots-input /Volumes/data/CodexProjects/T2/unfair-01-polynomials/computation/results/021-R8-L9-inputs.json \
  --mother-certificate /Volumes/data/CodexProjects/T2/unfair-01-polynomials/computation/results/021-R8-galois-certificate.json \
  --out-input /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuilt-a2b9-input.json \
  --out-options /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuilt-a2b9-options.json
```

`--parent-roots-input` 必须是实际含四个 `original_quadratic` 根块的 `021-R8-L9-inputs.json`，不是 galois certificate。复现命令已实际运行：生成的 a2b9 input SHA 为 `38b8e04c13a6e95f56d90f69841042b4ec551e9cefc9c6b99bb2d6900508aa1a`，options SHA 为 `15ece3343632c65b42183551ca912e568cbe968e366c7f80e781ddba6a51e33f`，均与 checkpoint 一致。比较生成 input/options 的 SHA-256 与 checkpoint 内 `roots_input_sha256`/`options_sha256`；不匹配时不得将缓存用于该 checkpoint。`scripts/verify_coprime_layout.py` 可先审查生成的 input 与 options 布局。所有归档文件小于 10 MiB；缓存特意不归档。

## 重跑数值搜索

在本课题目录运行；`--options` 使用上述再生且哈希匹配的缓存：

```sh
python3 computation/results/024-R8-coprime-gpu/scripts/stream_coprime_ds_gpu_v2.py --options /Volumes/data/CodexTemp/T2/unfair-01-polynomials/rebuilt-a2b9-options.json --roots-input computation/results/024-R8-coprime-gpu/inputs/a2b9-input.json --checkpoint /Volumes/data/CodexTemp/T2/unfair-01-polynomials/recomputed-a2b9.json --batch-size 131072 --checkpoint-seconds 30 --checkpoint-interval-batches 256
```

原运行是在四小时单批上界内因解析闭合提前停止；以上命令本身不设墙钟限制，调用者应在单批四小时内保存检查点、发送SIGINT并等待进程退出。停止记录仍为部分覆盖。
