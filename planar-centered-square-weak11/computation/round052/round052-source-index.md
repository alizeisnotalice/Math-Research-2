# R52 来源与证据索引

每个清单列原工作区来源、目标路径、哈希和范围；内容保持原实际执行记录，不因整合重新生成。以下路径相对于 `T2/.local/workspaces/<工作区>/planar-centered-square-weak11/`。

| 工作区 | 清单或交接路径 | SHA-256 | 作用 |
| --- | --- | --- | --- |
| pcw-position-slab-impl | `computation/round052/position-mass-integration-manifest.json` | `3949a749be24f044ef7dc126932fa6b2079c15ac5b46d3a26aeabb42fbbac1eb` | 全 bbox 扫描、候选及空盒/质量斜率 |
| pcw-sol-ledger-audit | `computation/results/round052/input-path-integration-manifest.json` | `82ac54904624d2c37ee23b4b9a011d80bf3708bcf24a2a0b79e37c238b3caf83` | 输入、路径、尾界独审 |
| pcw-r001-audit | `computation/results/round052/position-mass-audit-manifest.json` | `98388ecaee3ee99c04244d7f5acac8bfaee5f2e9b992dd4b5f419c63b8f6d697` | 独立候选审计，未重跑全部空盒循环 |
| pcw-kernel-audit | `computation/results/round052-endpoint-union/integration-manifest.json` | `ed00aa3b0b5053451a6a2755501ff2935146a8943ec73d6c62e50619f44f1763` | 完整端点序、并长及混合导数 |
| pcw-r001-geometry | `computation/round052/handoff.json` | `e06c1e8a8ced63a71ef12e6a1e915d4f820be9d5ddcca05fea9f1332c3d892cd` | 三个质量切片的原运行交接；局部墙另见下一行 |
| pcw-r001-geometry | `computation/round052/wall_2d_manifest.json` | `d0be0914b05543fc7f3ca2a292aaa7e31cb3146cc545ae6f4377ca4d99d13b3f` | 二维局部墙；固定顺序和路径，压缩端点 slab 与实际脚本记录 |
| pcw-r001-geometry | `computation/round052/midpoint_slab_gzip_manifest.json` | `9da952aeef00bb707d0f8a334cc0ba5486796f22faaf06f73b37fad29db4b0a0` | 中点 slab 无损压缩，保留原字节 SHA 与解压命令 |

核心解析证明 `proofs/joint-position-mass-cell.md` 来自 pcw-astra-variation，SHA-256 `3a7497334d88ec47f515a67f465660c490303654169b789b56703f6481ee9dce`；其最终清单 `computation/results/round052/variation-integration-manifest.json` SHA-256 `1389e8ce484db91b1ff53df3bc99daa83074716190c4a425b8461f770634e935`。以上索引不等于独立复跑或一般 C₂ 界；子清单若包含源文件绝对路径，整合后应以目标路径定位。全域积分面积是否在局部墙处有折点仍待核验。
