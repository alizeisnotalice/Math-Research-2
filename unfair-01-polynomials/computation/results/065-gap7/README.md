# 第 7 项缺口的因子次数 29 与 30

[完整证明与重放命令](../../../proofs/065-c7-degrees-29-30.md)。

Python 3.13.18、SymPy 1.14.0。判定使用整数、有理数和有限域精确算术。将输出写到包外可写路径。从本目录分别合并次数 29 的三个 gzip 部件和次数 30 的十一个 gzip 部件，使用 `--witness` 可逐项核验冻结乘子。次数 30 可从冻结的 `prefix29.30.json` 逐项重放，也可先完整重建早端并比较散列。

```sh
python -B verify_degrees_29_30.py --only 30 --scout prefix29.30.json --out /Volumes/data/CodexTemp/T2/a30-degree30-direct.json
```

主进程重建的前缀字节一致，重放次数 29 的全部显式恒等式，独立直接重算次数 30 的 77,977 个晚端字。计数与首次重算全同。详情见 `parent-replay.json`。两组冻结乘子已全部重放，验证 4,323 条投影恒等式和 29,888 条单位恒等式，详见 `archive/verification.json`。
