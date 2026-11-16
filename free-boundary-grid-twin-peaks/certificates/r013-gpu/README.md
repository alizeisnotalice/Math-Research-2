# 全GPU粒子核与独立大粒子重复

全部数据保持原均匀双峰条件律的τ分解。`audit/` 是新countsFull/cdfFull/respropFull核的实际CPU/单command-buffer/逐rank GPU对照，不是第12轮旧核审计。32个mode/case、96组完整比较、30个正常case的终端原实现比较通过。n2/n3覆盖所有τ、P513跨组；n7覆盖τ2/8/47/48、P17；两模式均含死亡分支，另加显式强制死亡案例。完整GPU生产未保存全标签，仍不声称n8全标签独立峰判定。

`stratified/` 是主Agent实际n8 P131072、每τ R32、τ2..64的full strat批，实际种子按模2^64计算：0xA022000000000000+0x7600000000000000+n*1000000+τ*10000+rep。Metal字符串与审计原源逐字节相同；host仅加每rep autoreleasepool。原始metadata、每rep种子、run_end、1Hz整机IOKit AGX和本进程RSS均保存。

`multinomial/` 是旧GPU选父核的n8 P131072 R32独立对照。其seed规则不同：0x5100000000000000+n*1000000+τ*10000+rep（不再加0xA022前缀）。`calibration/` 是full strat n7 P65536 R16，种子规则为0xA022前缀加0x6300000000000000后取模2^64。字段名seed_domain不能跨实现自动解释为同一算式；实际源和逐rep seed为准。

主Agent `computation/audit_r013_gpu.py` 用Decimal60核验全部4800个记录种子唯一、逐τ重复齐全、死端保留、比值和两种经验delta SE。τ分层SE和按rep配对SE不同，不混称；二者均非严格CI。Metal核一致性、文件哈希和实际GPU审计的cases/failures也检查。结果见root-independent-check.json。

n8 strat完整wall325.167685秒、Metal时间戳累计312.554533秒、提交等待累计322.885323秒，区间重叠不能相加。312个约1Hz整机采样均值98.6058%，281次100%、308次≥95%，最后退出附近含0；不是本进程独占，也不是连续时间100%认证。host目标RSS采样最大21616KiB不等同于全部Metal/驱动显存；完成后target33656、monitor33646均不存在。其他课题设备负载仍可能存在。n8 multinomial另批wall343.667422秒、全机平均61.76%，两批模式、负载和同步方式不同，不据此给出受控吞吐速度倍数。

来源/归档哈希见source-index.json。本目录不含二进制、缓存、被中断的部分6100批或重复种子的额外计数。n8其他粒子数及n9批在后续轮次继续，不计入本报告。

在stratified目录复现（注意325秒只是本次实测，不是固定耗时）：

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules tau_smc_pool.swift -o /Volumes/data/CodexTemp/T2/tau_smc_r013
python3 -B run_with_monitor.py /Volumes/data/CodexTemp/T2/r013-telemetry.jsonl /Volumes/data/CodexTemp/T2/tau_smc_r013 8 131072 32 2 64 /Volumes/data/CodexTemp/T2/r013-n8.jsonl full strat 7600000000000000
```
