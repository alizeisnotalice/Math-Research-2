# 两步端点核 CPU 对照

最终 b2-B8.csv 使用 n8,P4096,R32,8线程、每层8个二步块、基础 seed0xA0221F00。2016记录，基础提议总数8,321,499,136。中间态可无效，终点有效才接受。

```sh
clang++ -O3 -std=c++17 -pthread rmsmc.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-b2
/Volumes/data/CodexTemp/T2/a022-r018-b2 8 4096 32 8 8 /Volumes/data/CodexTemp/T2/a022-r018-b2.csv 0xA0221F00
python3 ../r018-cpu/audit_cpu_csv.py --block-size 2 b2-B8.csv
python3 ../r018-cpu/verify_traces.py b2 3 b2-n3-trace.mh.csv
```

线程局部变量修复、旧批排除和统计量口径同 ../r018-cpu/README.md。终端不移动。独立 trace 校验覆盖两个 elementary proposals 后的接受决定；没有宣称 GPU 生产已经采用此核。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。

归档版额外加入参数及 TRACE 单线程入口检查，未改变核。原运行源和归档源哈希见 source-provenance.json。主审实际确认 TRACE 两线程被拒绝；单线程 seed0xA0222000 的 CSV 和4284事件轨迹与原源逐字节一致，见 guard-and-trace-check.json。
