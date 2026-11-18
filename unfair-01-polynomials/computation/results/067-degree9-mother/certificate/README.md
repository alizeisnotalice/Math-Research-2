# 第67轮九次母式的双端入口与 S3/S6 完整排除

[中文覆盖与证明](proof-zh.md)。

重放只需 Python 3.13 标准库，断言须启用。输出目录应尚未存在且位于包外。

`/opt/homebrew/bin/python3.13 /path/to/package/replay.py --out /Volumes/data/CodexTemp/T2/round067-degree9-run --workers 8`

重放从空目录生成9份账本并逐字节与冻结版本比较。该包包含完整端点筛、严格根包络及四个单点位置的初段范围，完整排除仅覆盖位置3与反转位置6。第二、第四位置所附母式是首九行可行的局部探路例，固定这些母式已能由更深正向前缀树排除。

- `frozen/envelope.json` 有255种支持的端点筛结果与25个反转代表。
- `frozen/r3-full.rows.jsonl` 有1281个位置3母式和逐项根筛阶段。记录格式为 `[L,M,T,Delta,c,d,e,g,h]`，第三系数为 \((T+\sqrt\Delta)/2\)。
- `frozen/r3-trees.rows.jsonl` 有全部243个根筛通过项的正向树及每层规模。第二算法对每个二元预算反解对侧系数，逐层核对。
- `frozen/r3-independent.json` 记录从整阶坐标独立重算的完整匹配及空层分布。
- `frozen/singleton-probe.json` 与 `frozen/witness-verification.json` 记录小整数盒探路和具体双向前缀见证。
