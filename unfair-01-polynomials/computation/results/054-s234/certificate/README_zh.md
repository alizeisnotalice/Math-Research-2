# S234 紧凑证书包

本包记录七次实二次数域、非有理位置 \(\{2,3,4\}\) 的一个有限根筛和前缀树结果。它保存 58 个 root-screen field JSON 和 codes、68,593 个存活母式的死树索引，以及逐域独立 AB 重放回执。完整 HGQ trace 与数百万 rows5--7 前缀状态均没有复制进包。

`package-manifest.json` 对每个输入、冻结源、AB 回执和树索引分片给出 SHA-256。tree index 的每一行包含原母式、死层、最大 live 数、终状态 SHA 和完整 trace 的 canonical SHA。该 digest 是对既有完整树的压缩绑定，不是新的树重算。

包内可运行：

```sh
python3 verify_compact_index.py --package . --out compact-index-check.json
```

该检查核对所有包内 SHA、58 个 annulus survivor 集合与 68,593 条 tree-index row 的一一对应、AB 回执绑定及汇总层数 621,704、最大深度 71、最大 live 7。正式完整驱动另应现场重建 root 筛和 HGQ/AB 双树，不能将本 archive/index 检查写成新的树计算。

范围只限该固定 S234 母式盒；不证明一般均匀和猜想或其它非有理支持类。
