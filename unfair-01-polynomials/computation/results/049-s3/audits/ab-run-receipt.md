# S3 AB 精确前缀树：实际运行回执（不重跑）

本回执仅根据本 Agent 已完成的终端调用和冻结输出编写；本次没有重跑任何
全批树。所有输出、输入和程序仍在本目录，父方的 HGQ 程序未由本 Agent 修改。

## 解释器和调度

实际调用的命令名均为 `python3`，以 `nice -n 10` 在单核顺序循环中启动；调用
期间没有设置或改写 `PATH`。本回执生成时在同一 zsh 环境重新解析到：

```text
/Users/zzh/.local/bin/python3
Python 3.15.0 (main, Oct 14 2026, 00:56:12) [Clang 23.1.3]
```

标准批次的实际命令模板是：

```sh
nice -n 10 python3 /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round048-bfs/s3_exact_prefix_bfs.py \
  INPUT.json --out OUTPUT.json --depth 48 --state-cap 25000 --seconds-per-core 8
```

其中 `INPUT.json`、`OUTPUT.json` 依次如下。输出 SHA-256 是本回执生成时重新
读取所得；状态来自输出的 `status_counts`。

| 输入 | 行数 | 状态 | 输出 SHA-256 |
|---|---:|---|---|
| `inputs_snapshot.json` | 5,635 | dead 5,635 | `d139d930b0246259f80132c60ee92da16c0e2cf99afa7bcfc1f100258a52a361` |
| `inputs_delta_after192.json` | 2,983 | dead 2,983 | `d3e19ddc86c7294ce1ddec0c00e234407b978b768c3a1e6205d41f3a50ddc7f8` |
| `inputs_delta_after228.json` | 13,303 | dead 13,303 | `4e3d2fd8166b67d99bbcc504521a57eeef8606476998fac065294686989a3318` |
| `inputs_delta_after299.json` | 8,309 | dead 8,309 | `9cb86f0aceb5ce59b9767c85d2cc519374aac8d4d34005edbdd38b746daf3d7d` |
| `inputs_delta_after355.json` | 16,411 | dead 16,411 | `02376a6a4320924c1a1d92a80b577b8c19781cd21a967c957cf31992e28d0d73` |
| `inputs_delta_unique_after362.json` | 10,622 | dead 10,622 | `881afeaa0a501e65a48d9471593b3e0707abdc46d23bb79026e883865d60ebc6` |
| `inputs_delta_unique_after445.json` | 2,995 | dead 2,995 | `a9b352227ec36977cf0ddc266190480cc250e7eed56eaa617d3a91912a3692f0` |
| `inputs_delta_unique_after460.json` | 11,977 | dead 11,977 | `216a68fe86a85aa3964a06995a9136a58f3c15de35944f55be2f302c7f6ab21a` |
| `inputs_delta_unique_after507.json` | 39,577 | dead 39,576；depth-48 unknown 1 | `51e597c886bfe1a869185a9d29b7f37ccd3fd2da4af3635bc218783131e45d4e` |
| `inputs_delta_unique_final.json` | 201 | dead 201 | `5d9eb3561507d01af3cb86f73c9df905e25350f6ce2996e1b39803d2835793bd` |

早期四个“after”快照错误地把重复 `--exclude-snapshot` 传给单值参数，造成
21,922 条重复运行；这不改变任何单行树，但不能用各批行数相加作为覆盖数。
随后脚本修为可重复参数，且 `deduplicated_ledger_all669_ab_final.json` 对所有
重复核心逐层轨迹比较，冲突数为零。最终去重覆盖是 90,092，而不是上表行数和。

## 最后核心的续算

唯一深度 48 尚未死亡的行是
`[5,2,17,0,0,0,0]`。已实际调用：

```sh
nice -n 10 python3 /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round048-bfs/s3_exact_prefix_bfs.py \
  unknown_core_depth48.json --out unknown_core_depth512.json \
  --depth 512 --state-cap 25000 --seconds-per-core 30
```

它实际在第 53 层 `dead_all_prefixes`，没有达到任一阈值：

```text
output: unknown_core_depth512.json
SHA-256: 36ab343beaa9382b9e83f23054f1adb6e53029cc02de84545311b3df5aa4eb45
```

`make_bfs_override.py` 将此行的 53 层记录绑定回原批次索引 19,076；它不是
覆盖或修改历史的 48 层输出。

## 退出核验

本回执生成后执行了
`pgrep -af 's3_exact_prefix_bfs.py|python3.*round048-bfs'`。除该 `pgrep` 自身
外无匹配运行进程；因此本 Agent 启动的 AB 树进程均已退出并释放内存。
