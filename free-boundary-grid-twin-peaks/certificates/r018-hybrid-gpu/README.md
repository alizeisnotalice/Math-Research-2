# GPU 前缀模块审计原型

这是 CPU 负责 count/CDF/选父/添点、GPU 负责 prefix copy/append/move 的混合原型，不是全 GPU 生产版。
27份 strat-n3-t*-l*.json 覆盖 tau2..10 和 L0/4/16，P17；每份记录准确 seed。有效 GPU 状态字与 CPU 邻接重建逐层比较。事件列表是 CPU 同 key 重放，GPU 证据是最终 prefix/bitset，不能称 GPU 直接逐事件日志。

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules integrated.swift -o /Volumes/data/CodexTemp/T2/a022-r018-hybrid
/Volumes/data/CodexTemp/T2/a022-r018-hybrid 3 2 17 a0220000002e14e1 4 /Volumes/data/CodexTemp/T2/a022-r018-hybrid.json
python3 verify_trace.py /Volumes/data/CodexTemp/T2/a022-r018-hybrid.json
python3 compare_baseline.py
```

compare_baseline.py 实际编译前轮 rank 核，在临时卷运行，检查 n3 全秩 L0 原子状态及父索引。归档 old-n3-* 是实际原核输出。旧核终端存前终端状态，比较器明确按此对齐。27份验证通过不证明大尺寸混合。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。
