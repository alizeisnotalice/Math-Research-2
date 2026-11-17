# n=24 merge-RB r029 收尾核验

范围：对 r028 的单一有限 GPU 批进行完整性与实现不变量复核。此记录**不**认证混合、SMC 精度、原条件律概率、锥条件或渐近结论。

## 引用的 r028 原始材料（未改写）

- Raw: `../a022-merge-rb-n24-r028/raw/n24-p262144-r32-l16.jsonl`
  (`SHA-256 80afda8d82e1b72a74efaf8e7a9d4a0d188233f567d0a538172314fa04f28047`).
- 生产回执: `../a022-merge-rb-n24-r028/receipts/n24-p262144-r32-l16.receipt.json`
  (`SHA-256 a4269b388ff06c097f9c3b250040c6af3949cae36b9b7cb3e2f583365a5fe176`).
- Finalizer 回执: `../a022-merge-rb-n24-r028/receipts/n24-finalizer-status.json`
  (`SHA-256 b3613bd4a601e920221c4d30376e10b794ea3630e3f4ebac91126fcec22770a2`).

参数为 `n=24,P=262144,R=32,tau=2,L=16,strat,domain=9e270e0000000001`；冻结派生源 SHA-256 为 `3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6`，二进制为 `5de42d7ec2eb946cfb805b59df658cf2f25d87e6018d4561d3e7a7b3a2a9d598`。

## 完整性与退出

Raw 共 66 行：metadata、32 个唯一 ID `0..31` 的 replicate、32 个唯一 ID
`0..31` 的 RB 记录和 `run_end`；32 个 seed 亦均唯一、`logZ` 有限。生产回执
为 returncode 0。finalizer 已以 `driver_pid_observed_gone=82830` 写出 passed；复核时
driver 82830、time 82831、caffeinate 82832、GPU child 82833 和 finalizer 83002 均不在
进程表，故没有遗留本批 GPU/保持唤醒进程。

本轮另以冻结硬化聚合器（SHA-256
`f1593cc93db21d778289966c80539606b079d702fa3fd5c197394c27d4737389`）独立重放并通过；
输出 `n24-production-audit-replay.json` 的 SHA-256 是
`43091040f3906263cc28846dd43d69dee189d513aeb06a552228d779505852f3`，与 r028 finalizer
所记录的 audit SHA 一致。

## 有限诊断与资源

稳定 log-normalizer 加权的有限诊断为 `ESS=3.48587943`、最大权重份额
`0.35925076`；m-lift ESS `3.48984623`、最大份额 `0.35812660`。权重集中严重，因此
`Q_mRB=0.36671163`、`Q_L_legacy=0.37402068` 和 `A2/B` 的 RB 数值
`2.01149345` 不能当作已认证概率估计。

实测 production wall 为 2818.68 s，最大 RSS 2,110,537,728 bytes，peak footprint
3,116,239,344 bytes，0 swaps、0 received signals；已完成子进程均退出。
