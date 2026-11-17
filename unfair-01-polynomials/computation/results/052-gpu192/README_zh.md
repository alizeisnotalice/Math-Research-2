# m=96、n=144 的 192 样本数值审计包

这是一个已完成 MPS 数值探索的可复现审计包，不是严格因子分解、有限排除或原猜想的
反例。它保存 12 个重启窗口各 16 条 JSON 候选（共 192 条）和两个实际审计回执；所有
“精确”表述只针对已保存的十进制 JSON 值、有限前缀池及源文件哈希。

参数为 `m=96,n=144,d=15,starts=4096,steps=72000,restart_every=6000,seed=2026111707`。
运行使用 MPS/float32；`release.json` 记录完成和释放请求。

## 覆盖语义

前缀池大小为 46,011。12 个窗口每个分配 4,096 个起点；循环窗口的并集覆盖整个池。
但 `12*4096-46011=3141`：最后一个窗口绕回并重复 3,141 个槽位。因此窗口**不是互不
相交的**，任何“disjoint”旧元数据都不作为本包事实。

## 包内实际复现命令

在 `package/` 目录运行（没有 `--input` 参数）：

```sh
python3 fraction-verifier-from-r040.py \
  --run-dir run --worker worker.py --prefix-module fractional_prefixes.py \
  --record reproduced-fraction192.json

python3 generic192-verifier.py \
  --run run --fraction-verifier fraction-verifier-from-r040.py \
  --worker worker.py --prefix-module fractional_prefixes.py \
  --prefix-plan prefix-plan.json --base-source frozen-base-worker.py \
  --out reproduced-generic192.json
```

第一个命令独立枚举 d=15 的有限三元前缀池，并以 `Fraction` 复算 192 个保存样本。第二个
命令重新核对源绑定、12 窗口循环覆盖、重启快照、保存的 JSON 值和精确十进制残差。
`manifest.json` 给出所有包内文件哈希；`receipts/` 是原始两份真实回执，
`reproduced-*.json` 是本包重放的输出。
