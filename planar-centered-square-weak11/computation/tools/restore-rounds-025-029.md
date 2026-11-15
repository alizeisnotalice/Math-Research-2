# 第25–29轮归档复现指南

适用范围：恢复五轮保存的搜索账本和最终候选，或在独立位置重新运行固定候选的 R=10 两核。恢复文件、检查哈希和复算均不等同于重新验证各轮全部搜索面积。使用归档的原脚本及其同目录依赖；不要修改原脚本或覆盖归档证书。须先确认 `/Volumes/data` 已挂载，并将缓存与输出放在数据盘。

以工作目录为课题目录、`SCRATCH` 为**尚不存在**的新目录为例：

```sh
SCRATCH=/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/restore-example
python3 computation/tools/restore_rounds_025_029.py --archive . --output "$SCRATCH"
```

工具将各轮 `manifest.json.gz` 和 `results.json.gz` 解压成原始字节，将已归档的 `summary.json`、第26–29轮 `completed.jsonl`、各轮最终输入和第8轮基线输入复制到对应历史目录结构；`RESTORE-MANIFEST.json` 逐文件记录来源、状态及 SHA256。第25轮未归档原始 `completed.jsonl`，工具只按保存的有序 results 重建 `completed.reconstructed.jsonl`，该文件不是原完成顺序，也不声称与原文件字节相同。某些账本脚本要求名为 `completed.jsonl`：如果需要试跑，先在隔离副本中复制重建文件为该名称，并将其标记为重建；不得覆盖历史证据。

每个 manifest 的 `tasks` 保存完整输入；结果表里的 `input_sha256` 是对规范 JSON 的摘要，与最终候选文件的字节 SHA256 有不同含义。工具只复制每轮最终候选，并不生成其他候选文件。若需其他候选，依据 manifest 任务的 `input` 明确生成，并按脚本的规范序列化方式核验摘要，标为重建。对早期搜索结果的复用应核查对应任务、来源哈希和尾界，不应把复用记录中前轮组合标签当成本轮组合标签。

前置关系：R25 用 R8 基准和 R23 质量输入；R26 用 R8、R25 的原 manifest/results 与最终输入；R27–29 均读取 R26 manifest/results，分别以 R25、R27、R28 输入为基线，R28/R29 还读取前一轮 summary/results 的复用记录。各轮需要 `periodic_halo.py`、`independent_periodic_halo.py`、`periodic_halo_dominated.py` 及各搜索脚本的本地 imports；R23/R25等具体引用和 SHA 以各 manifest、轮次报告和独审文件为准。

## 固定候选复算

先将证书输出目录设为一个新建、空的目录，例如 `OUT="$SCRATCH/check-r29"; mkdir "$OUT"`。从课题根目录运行：

```sh
python3 computation/periodic_halo.py computation/inputs/r029_final_0298.json --bits 40 --tail-radius 10 --output "$OUT/primary.json"
python3 computation/independent_periodic_halo.py --input computation/inputs/r029_final_0298.json --bits 40 --radius-cap 10 --output "$OUT/independent.json" --compare-primary "$OUT/primary.json"
```

将输入替换为相应的 `r025_final_0110.json`、`r026_final_0360.json`、`r027_final_0079.json`、`r028_final_0026.json` 可复算其余四轮；每轮使用各自的新输出目录。准确面积、内外比值与原归档核对，耗时、绝对路径、平台及输出文件哈希不必一致。R=10尾界及 \(q=3\) 应同时按相应证明核验。

## 搜索程序

R25 `scan_fixedmass_positions.py` 和 R26 `scan_quadrant_groups.py` 提供 `--help`，可在新目录中指定输入、前轮账本及输出。R27–29 的原脚本为保留来源证据而包含历史绝对路径及固定输出目录；**不要直接运行这些原件**。若要重跑，应复制脚本和依赖到隔离工作区，仅在副本中将 `P`、`basepath`、`targetpath`、`src` 和 `out` 指向恢复目录及新的空输出，并记录副本与修改后 SHA256；相关脚本会要求 `/Volumes/data` 挂载。搜索会使用多个进程并重算大量面积。`build_round025_ledger.py` 等历史账本脚本也含固定路径，需以同样方式在副本中调整，且R25的重建完成流不可误称原始流。

本指南和恢复工具的实际验证：Python 3.15.0 下执行恢复，成功输出26份逐项标记文件；未执行大几何重跑。输出用于可重复检查，原始轮次证书及报告仍是归档依据。
