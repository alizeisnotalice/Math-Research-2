# r012 GPU-parent 大粒子基线（独立种子域）

本目录只含一次已完成的原条件律按 \(\tau\) SMC 数值基线；它不是关于远峰概率
的严格证明。模型与 r012 GPU-parent 冻结核相同，唯一的运行域变化是新的 seed
domain `0x5100000000000000`。Metal 字符串的 SHA-256 为
`886418b43df857752d956a9f5f675ac18e645d0d1329ced0fc246905eabc5d82`，与来源
`/Volumes/data/CodexTemp/T2/a022-original-smc-gpu-opt-r012/handoff-r012/source/tau_smc_integer_gpu_parent.swift`
逐字节相同；宿主侧只增加参数 guard、seed 记录及每重复的 autoreleasepool。

## 实际运行

```sh
swiftc -O -framework Metal source/tau_smc_integer_gpu_parent.swift -o tau_smc_gpu_parent
python3 source/run_with_monitor.py telemetry/n8-p131072-r32.jsonl \
  ./tau_smc_gpu_parent 8 131072 32 2 64 raw/n8-p131072-r32.jsonl
python3 source/summarize.py raw/n8-p131072-r32.jsonl > raw/n8-p131072-r32-summary.jsonl
```

- Apple M6；\(n=8\)、\(P=131072\)、每层 \(R=32\)、\(\tau=2,\ldots,64\)。
- 全部 (63\times32=2016) 个重复完成；没有因资源压力降低 R 或 P。
- `raw/n8-p131072-r32.jsonl` 为 2,018 行（metadata、2,016 重复、run_end），
  `run_end.wall_seconds=343.6674221`，退出码为 0。
- 每个重复均记录 seed；所有 seed 属于新域，不与旧 `0xA022…` 域混用。
- n=2、P=64、R=2 smoke 精确给出 \(Z=1/3,U=1/6\)，见 `raw/n2-smoke.jsonl`。

## 数值摘要

全 \(\tau\) 合并的未归一化重复均值为

\[
\bar Z=9.611435136069461\times10^{-15},\qquad
\bar U=3.2209178021046536\times10^{-15},
\]

故有限样本比值为 \(\bar U/\bar Z=0.3351130977\)。`summarize.py` 保存五个
原始矩 \(\bar Z,\bar U,\overline{Z^2},\overline{ZU},\overline{U^2}\)；把32个
重复当作理想独立重复时，plug-in delta-SE 为 0.0005958578。该 SE 只描述该
理想重复模型：有限 PRNG、SMC 退化、跨 \(\tau\) 合并和比值偏差均不由它认证。

\(\tau=63,64\) 各有32个全灭重复，按零 \(Z,U\) 保留；没有删去幸存者再归一化。
这并不是这些层理论不可行的证明。

## 资源与遥测

所有重复的 Metal hardware GPU 时间戳和为 289.157374 秒；同步
submit+wait 时间和为 342.023753 秒；它们是跨重复累计，不能与 343.667422 秒
全进程墙钟相加。`telemetry/` 每秒使用 IOKit `AGXAccelerator`
PerformanceStatistics 采集整机 device/renderer/tiler utilization 和目标进程 RSS。
共334个样本：device utilization 0--94%，均值 61.76%；RSS 5,616--17,520 KiB。
这些是整机、非独占 AGX 读数，不能解释为本进程利用率或“跑满 GPU”。`powermetrics`
在本机需超级用户权限，未调用。

运行结束后已检查 target 和 monitor 进程均不存在。

## 文件与限制

- `source/tau_smc_integer_gpu_parent.swift`：本次实际宿主源；Metal 核未变。
- `source/run_with_monitor.py`：前台子进程及一秒非独占遥测。
- `source/summarize.py`：五矩与描述性 delta-SE。
- `raw/`：n=2 smoke、n=8 原始 JSONL 和摘要；`telemetry/`：运行遥测原始 JSONL。

本次没有完整 GPU 标号轨迹，因而不把根距矩数值称为独立峰构造认证。数值结果也不与
任何旧批直接合并成严格置信区间。所有原始文件和源的 SHA-256 见 `HASHES.sha256`。
