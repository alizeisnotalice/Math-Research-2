# 次数64、模163平方次数7边界的独立搜索

[固定41位的证明与余式完备性](proof-zh.md)。

实现采用普通变量 \(x\) 中 \(K_q^2\) 的余式、五坐标投影桶及桶内排序。对投影相等的每对掩码，C++ 核全部14坐标。Python 另在负倒变量中重建65列、核仿射目标并独立复决所有投影碰撞。六个小例同时与既有 C++ 及逐掩码长除比较，共6474个掩码，有三个命中例及三个空例。

## 重放

使用 Apple clang 21.0.0、Python 3.13.18。下列输出路径应尚不存在。

`sh /path/to/package/build.sh`

`/opt/homebrew/bin/python3.13 -B /path/to/package/replay.py --out /Volumes/data/CodexTemp/T2/degree64-independent-replay --workers 4`

`/opt/homebrew/bin/python3.13 -B /path/to/package/check_receipts.py /Volumes/data/CodexTemp/T2/degree64-independent-replay --compare-deterministic /path/to/package/frozen-index.jsonl --out /Volumes/data/CodexTemp/T2/degree64-independent-audit.json`

边界模式覆盖完整 \(q=1,\ldots,162\)。每个参数记录模数、全部65列、仿射目标、所有投影及全余式掩码、两边遍历计数、耗时、峰值内存和退出码。断言必须启用。校验脚本独立核对162个检查点及散列。确定性索引不包含 PID 或耗时。

## 实际账本

最终异址重放建表169,869,312个半边并扫描339,738,624个半边。10,530列均经独立负倒余式检查，3,276个投影碰撞均经全余式复决，全余式命中为零。四个本机工作进程并行，各参数独立启动和退出。全参数阶段4.403656秒，连同小例的外层回执5.123秒，单子进程实测峰值13,877,248字节。初轮全量搜索用同样的41位，最终逐参数核对所有公共字段无差异。

`frozen-index.jsonl` 含162个参数的确定性散列和计数。按参数逐个保存的完整结果在外层 `full-tight-v2/by-q`，完整索引为 `full-tight-v2/all-q.jsonl`。

这个排除使用第65轮的低端次数界和模163平方次数恰7时的反向偏移。它保留 \(C_{12},C_{13}\) 的两个取值，也允许方程中的一般二元余位。模163平方次数更高、\(C\) 在有理数域带重因子及其他总次数仍需别证。
