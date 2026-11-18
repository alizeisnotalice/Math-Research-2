# 第 7 项缺口约束的逐行复用

[第 12 至 16 行的直接证明](../../../proofs/066-c7-forced-c12-c16.md) · [第 17 至 30 行及次数 64](../../../proofs/066-c7-forced-bits-to-30-and-degree64.md)。

独立审计及重放覆盖第 17 至 24 行。主进程审计第 25 至 30 行，在 later-endpoint-audit.json 记下源文件和散列。旧证书在相同的内部零锁下继续成立。

    python -B verify_c12_c16.py --out /Volumes/data/CodexTemp/T2/a30-c12-c16-new.json
