# GPU 选父迁移的冻结实现与核验

此包将 host 上逐 child 的随机整数目标、CDF 二分选组迁至 Metal。每组256粒子的整数和仍由host构造CDF；组内扫描、mutation及随机键与原实现相同。它仍是多项式重采样，不是整数分层方案。

`source/` 是实际生产源码；`raw/` 是n2 smoke、n7 P16384 R8、n8 P32768 R16。seed=0xA022000000000000+n*1000000+tau*10000+rep。冻结原型只用于已记录合法参数；5个uint掩码仅支持N≤160，未宣称任意CLI输入安全。后续guard版本和更大粒子生产另归档。

`audit/gpu_state_parent_audit.swift` 实际dispatch旧和新Metal核，对n2 tau2..4、n3 tau2..9 P513，以及n7 tau2,8,47,48 P17，比较每秩counts、parentOut、全state及terminal moment。15例通过；主Agent独立编译重跑结果相同。rank_audit.py是另外的CPU参考，不能替代实际GPU审计。此包不认证下一轮的新full-command-buffer核。

主Agent的Decimal60逐重复比较及统计在相邻r012-tau-smc/root-aggregate-audit.json；n7新run缺恒零tau=N，核验显式补零。same-seed重跑不能计新独立样本。n8墙钟212.987→109.805秒、GPU时间戳66.176→68.171秒。时间重叠，不相加；覆盖比不是设备利用率。本轮此baseline没有连续AGX记录。

复现审计（输出和二进制放data临时目录）：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules audit/gpu_state_parent_audit.swift -o /Volumes/data/CodexTemp/T2/gpu_parent_audit
/Volumes/data/CodexTemp/T2/gpu_parent_audit audit
```

所有已完成的本轮程序均退出。源码及原始输出哈希/来源见source-index.json。
