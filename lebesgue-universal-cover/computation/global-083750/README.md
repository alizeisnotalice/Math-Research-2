# 0.83750 规范原叶外围：种子至四层的可移植包

## 结论和证据范围

固定面积分母 `2^57`，严格整数门槛 `120696470013529293`。从公开原树的 2,048 根中，以共同刚体、反射和自身对称归一化得到的 1,024 活动片内，6,710,314 个规范闭原叶中，继承的逐叶下界给出 699,060 个新靶种子。四层从这些种子依次保留 136,243、22,418、8,792、5,442 个原叶，分布于 170 个分片。故这套五体外围整数证书在全规范域排除其余 6,704,872 叶。每层耗尽预算的叶整体保留；末叶不构成低面积可行证明，也没有由本包进行六体末端关闭、缩体到最终目标的转移或一般万能覆盖下界证明。

| 层 | 输入 | 通过 | 保留 | 几何访问 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 699060 | 562817 | 136243 | 31860339 |
| 2 | 136243 | 113825 | 22418 | 46097404 |
| 3 | 22418 | 13626 | 8792 | 94226356 |
| 4 | 8792 | 3350 | 5442 | 167820059 |

本轮实际从早期公开桶、.83731 分卷和原树逐片重建 .83750 原始索引及径向算术；实际读取并取最大值合并已有 .83732、.83735、.83740、.83741、.83742、.83743、.83744 和 .83745 整叶整数回执；实际运行了本靶四层 340,004,158 次几何访问；另以 Arb 重核新靶的缩体根外界。本轮未重算全部历史几何核；历史整叶下界、径向材料及其包含证明沿用先前证明。整合时对全部 4096 片检查输入、输出、集合分割和门槛，另单独对最终余叶与第四层记录逐号比对。散列只绑定文件，不能替代面积核与其几何正确性证明。

## 输入及环境

- `COMP` 为 Math-Research-2 对应检出版本的 `lebesgue-universal-cover/computation`，本轮源检出 HEAD 为 `3d77a9e11396965fa223b9fd49ddcd67a6103cba`；四层代码 `replay_tomar_signed_translation_screen.py`、`replay_tomar_083731_scored.py`、`five_grouped_localization.cpp`、`tomar_signed_translation_screen.inc` 均没有本地改动。需要整个计算目录及其支持源码和公开资料。
- `TREE` 从 [Tomar 公开树](https://doi.org/10.5281/zenodo.23717912) 的 `cert11d_08366.tar.xz` 解包，文件 SHA-256 `6b5d28bc488edcfde9c6a49f3cb364f9086f445a58f85f166da23727036cf122`；解包后 `SHA256SUMS` SHA-256 `56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410`。
- `T32` 使用 COMP 的 `tomar_fan_table.py TREE T32 --subdiv 32` 构造；文件 SHA-256 `fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c`。
- 实际 Python 3.13.18、NumPy 2.5.5、python-flint 0.9.0、Apple clang 21.0.0，`-O3 -std=c++17 -ftrapv`，启用断言。约 12 个本机逻辑 CPU 可见；其他作业负载需实际检查。每片检查原树文件与输入哈希。

## 从种子包独立重跑四层

在任意工作路径以 `PY` 指向上述 Python，`PACKAGE` 指此目录；`INPUT` 和 `FRESH` 均为不存在的新目录，`CACHE` 是 Clang 模块缓存。工作目录路径及参数均可任意指定：

```sh
"$PY" -B "$PACKAGE/prepare_replay.py" "$INPUT"
"$PY" -B "$COMP/replay_tomar_signed_translation_screen.py" \
  "$COMP" "$TREE" "$T32" "$INPUT/seeds" "$INPUT/seed.json" "$FRESH" "$CACHE" \
  --workers 2 --until 4
"$PY" -B "$PACKAGE/verify_replay.py" "$FRESH"
```

末步逐一重新读取四卷内 8,192 个文件并与新跑结果逐字节比较，检查每层账本和最终位置。动态库二进制散列因编译环境变化可不同；几何源码、原始输入和结果记录不能变化。`prepare_replay.py` 已本轮实际运行；`verify_replay.py` 本轮以原生产输出检验全部包成员，**尚非独立几何重跑**。重新计算所得结果必须先实际完成才能声称二遍。

## 输入溯源与外域

`records/` 记录本次基准与合并报告、四层 ledger、各层审计、缩体根外 Arb 结果及最终清单。公开上游脚本为 `COMP/audit_tomar_multi_target_scores.py` 以及 COMP 内对应历史分卷、ledger 和原树；本包 `prepare_stage.py`、`merge_old_scores.py`、`merge_45.py`、`prepare_seed.py` 提供从上游已公开输入重新组包的相对参数化入口。历史 .83745 四个 `native45-screen` 卷见同仓库的 `COMP/six-input-083745/volumes/`；其他历史分卷由 COMP 提供。实际执行 `python -B PACKAGE/rebuild_seed.py COMP TREE CACHE NEW_WORK`，从原树、公开历史分卷与该目录内的上游四卷重新构造索引和种子，1024 片和两份报告与种子包逐字节相同。此重建复用历史整叶面积分数，未重算历史几何核。独立种子至第4层复算只需 `seed50.zip`，但此种子作为输入沿用历史面积证明。

圆盘及四个正定宽多边形以各自中心缩放至 `0.999999006`。`audit_root50.py COMP TREE` 的本次 160 位 Arb 计算得到四个根平移方框外界分别大于 0.85914435、0.94373319、0.96463856、0.97299840；高于本靶。另四个新径向截断的整数下界均不小于本靶门槛；根外排除的证明需要该脚本所引半径公式以及原径向材料包含证明。闭角规范化和原树完整性沿用已发表工作区证明；跨区边界按闭叶覆盖。

## 文件指纹

- seed50.zip `9738fff07eaed42db51bf27839e4ed76c457ade80d7a1e5352efb9f0e5c31905`
- screen1.zip `b78fa27548a7fc95c269d1901418aa6016a5a2e2619e9ee3e1d073eeabcceb9e`
- screen2.zip `d987d17b50ec3f37c04e9b6682952658096501943e60375cb3abb602ce066edc`
- screen3.zip `5f87230a15f255984260750a43753bbcfe0ac0308ce616d05c4016bb95530437`
- screen4.zip `902f6b22beac1a28a5e07236535285c7e26d570845db033dd8fb544bfbedfc0f`
- 四层总账 `3e83c009e5776e2af17be09d3d712d0cd5bc83236df439ed02cb83559fec1dd9`
- 最终清单 `06c7b4202f85d404c5e04cd53d19f4e4bd51a81da8e355e85f1f05a2913dcbc0`
