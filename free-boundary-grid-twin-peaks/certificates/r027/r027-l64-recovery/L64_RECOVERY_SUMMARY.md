# n=24、L=64 续跑恢复与审计（r025）

状态：**有限样本/实现审计通过；不是原条件律的认证估计，更不是渐近证明。**

## 恢复链与完整性

原批（`n=24,P=131072,R=32,L=64,domain=9e29070000000000,strat`）因 2026-11-17 系统重启中断。原始只读文件恰有 logical replicate `0,1,2,3`，无 `run_end`。在核对种子公式后，以 `domain+4=9e29070000000004` 运行 28 个局部 replicate；其局部 ID `j=0,...,27` 对应原计划 logical ID `j+4`。小型 `P=257,R=2` 映射试跑已通过。

`reconstruct_canonical.py` 将这两段不改写地拼成唯一的 32 个 logical ID，并已审计：

- JSONL 共有 metadata + 32 unique replicate + run_end 共 34 行；所有 seed、有限 `logZ`、512 个 group/replicate 均齐全。
- 每个 group 的整数恒等式 `F-2G=4D` 精确为零，`D>=0`；距离直方图整数和等于对应 `L` group，前两距离桶均为零。
- 冻结聚合器的 `epsilon=1/2` 重放字节哈希与已存报告相同。

此处整数非负性和恒等式仅是本有限输出中的必要结构诊断，**不构成**问题所需的锥条件、混合认证、SMC 精度界或渐近结论。

## 有限加权诊断

冻结的 `aggregate_tau2_distance.py` 对全图 K=2 目标的 lifting 加权计算（不是未加权的 tau=2 链经验频率）给出：

| 量 | 数值 |
| --- | ---: |
| `Q` | 0.3568378028039002 |
| `deltaSE_Q` | 0.0470304137238380 |
| `E2_tau` | 2.894093160702812 |
| `deltaSE_E2_tau` | 0.195443694302450 |
| `P(d>6)` | 0.999999982238368 (`deltaSE` 2.085e-8) |
| `P(d>12)` | 0.980822052950935 (`deltaSE` 0.0131971) |
| `P(d>24)` | 0.00295015361301729 (`deltaSE` 0.00135003) |

这里的 `deltaSE` 是该有限重采样批的内部标准误诊断，**不是置信区间**。权重集中严重：tau=2 normalizer ESS 2.63959、最大份额 0.582886；full-lift ESS 2.97934、最大份额 0.523583。因此这些加权数值不能作为原条件律概率的精确估计或反证依据。

相比独立的旧 L=16 批（不可 pool），full-lift ESS 从 1.44990 提至 2.97934、最大份额从 0.8080 降至 0.5236；仍远不足以消除权重集中问题。

## 资源和退出

28-replicate 续跑自然返回 0。实测 wall 5383.25 s，最大 RSS 1,222,672,384 bytes，peak footprint 1,531,889,224 bytes，0 swap、0 received signals。驱动、GPU 子进程、`caffeinate` 和 finalizer 均已退出；没有重启原始中断批。

## 核心可复核哈希

- 原始 partial raw: `fbac0ef683f0677ed7f0cb18b4e884a3083ab20d7d3c59894bd17b26bdd990bc`
- 续跑 raw: `785468ba61844832befe0ffdb6327670d92bd30cbd6472771183fa7bf2726101`
- canonical 32 raw: `5516568b7bb3aec65b5432553a5615a74a1def5b7357e797ae50ea1796fa51ff`
- 映射账本: `e602e8792b27ebe5c3a1594120b2b449bf5added68260064da626da015d70e83`
- finalizer status: `8dd1ad16b80b3821f59d174eae21936367b07aa25b620b71771efa839ca042a0`
- 冻结终端源 / 二进制 / 聚合器: `01b77bca1735f07293eae5ffeea0326810fdc5ea785c7a9e99d9030cba9b0702` / `137aa827fd2fcfbacd9a40aaa7059ee9dd8f38a3616e0c849724d2143f10f2cd` / `cb466dcae1628f81ff2e198e65d2f95a71b198f77cd60fe8e717c9076403c358`.
