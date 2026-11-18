# 四设备复算
在对应子目录的临时副本运行。Apple Silicon 编译使用现有共享 Swift/Clang 模块缓存。CPU 求积用单线程 BLAS。
- local-cpu 用 Python 3.13.18 执行 `python area.py --cases 96`。
- mac-cpu 的 15 个并行分片分别执行 `python folds.py i --exponent 14 --seconds 46 --out shard-i.json`，以 0 至 14 替换 i，采用 Python 3.12.13。exact_folds.py 与 fold_precision.py 分别复核有理推导和高精度积分。
- local-metal 将 quartic.swift 以 Swift 6.4.2 的 -O 模式编译为 quartic，同目录放置 quartic.metal。运行 `./quartic 45 88 main 0`；以 `python verify.py quartic-main.json cpu-main.json` 复核。
- mac-metal 将 scan.swift 以 -O 模式编译为 dual_scan，同目录放置 dual.metal。运行 `./dual_scan 23 23 17 29 65 17 1024 1024 sync.json`，再运行 `python check_cpu.py sync.json sync.cpu.json`。用标准库运行 rational_contraction.py 可另行核对严格不等式。
