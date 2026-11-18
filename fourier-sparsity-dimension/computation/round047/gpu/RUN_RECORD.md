# 第47轮双机九变量三次 ANF 的 Metal 实跑

2026-11-18，中国标准时间。两机运行了不相交的计数区间。每个计数值直接指定前 32 个三次 ANF 系数，其余系数由所列种子和计数值生成，故本批四亿个函数互异且次数为三。每个函数的 512 点整数 Walsh 变换由 Metal 执行。两条命令队列和四个在途缓冲分担真值生成、变换和绝对谱直方图归约。

| 设备 | 候选数 | Walsh 系数 | 程序墙时 | 本任务的 Metal 命令时间并集 | 整机利用率至少 90% 的采样 |
|---|---:|---:|---:|---:|---:|
| Mac mini，Apple M6 | 160,000,000 | 81,920,000,000 | 17.470 秒 | 17.301 秒 | 91/92 |
| MacBook，Apple M5 Pro | 240,000,000 | 122,880,000,000 | 40.031 秒 | 33.907 秒 | 206/206 |

两机通过同一 UTC 启动门调度。按两机记录的时钟，进程窗口交叠 17.470 秒，Metal 命令时间区间并集交叠 14.607 秒。该交叠使用两机各自的系统时钟换算。整机利用率来自 IOAccelerator，包含其他进程。两台设备在本批开始前已分别测得 99% 和 100%，故整机采样不能作为本任务独占 GPU 的证据。任务自身的 Metal 时间来自完成的命令缓冲，两队列各命令持续时间相加可能大于它们的区间并集。正式任务和随后抽检进程均已正常退出。

Mac mini 批次的每隔 32 个缓冲抽查首、末 ANF，共直接核对 156,672 个整数系数。MacBook 以同样方式核对 234,496 个系数。GPU 产生的全量直方图经过总量、Parseval 能量与非零系数计数对账。独立 Python 从两机的生产抽样编号各复算 8 个 ANF，并在其中各 2 个编号上另跑单函数 Metal 变换核对绝对谱直方图。独立复算结果及小对照的文件哈希在 `VERIFIED_SUMMARY.json`。负载测量与有限谱统计均限于这两个固定计数区间。

## 复现

将本目录中的 Swift 源和监测脚本放入机器上已有的临时目录。用该机的既有 Swift 模块缓存编译。下例以 Mac mini 的缓存为例。

```sh
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-module-cache \
  R47CubicWalshLoad.swift -o r47_cubic_walsh_gpu
mkdir -p m6-repro m5-repro
./r47_cubic_walsh_gpu m6-repro/pilot.json 1 32768 a4700001 8192 1
./r47_cubic_walsh_gpu m5-repro/pilot.json 50000001 32768 a4700002 8192 1
```

直接重跑各自的正式批次可以使用以下命令。每条分别在相应设备上运行。

```sh
./r47_cubic_walsh_gpu m6-repro/cubic_sweep_result.json 200000001 160000000 a4705001 32768 32
./r47_cubic_walsh_gpu m5-repro/cubic_sweep_result.json 600000001 240000000 a4705002 32768 32
```

需要同期利用率回执时，两机各运行一条对应的监测命令。两边写出 `ready.json` 后，分别将同一个晚几秒的 UTC Unix 秒数写入各自输出目录的 `start_epoch.txt`。正式输出目录在启动前应不存在旧闸门文件。

```sh
python3 monitor_cubic.py m6-repro 200000001 160000000 a4705001
python3 monitor_cubic.py m5-repro 600000001 240000000 a4705002
```

使用 Python 3.10 及以上核对归档结果。可选参数 `--sample-gpu-binary` 指向已编译程序时，会再核验并保存四个单 ANF 的 GPU 对照。常规重验不需要重跑正式 GPU 批。

```sh
python3 verify_gpu_sweep.py m6 m5-pro /path/to/rechecked_summary.json
```
