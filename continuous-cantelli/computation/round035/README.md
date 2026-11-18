# 双机计算复现

在对应子目录的临时副本中运行。结果中保留原始批次标识。采样计数分别对应不同模型与不同求积阶数，两次本机 GPU 扫描共享确定性采样前缀。

* local-cpu 使用 Python 3.13.18、NumPy 2.5.4、SciPy 1.18.1 和 mpmath 1.3.0。分别运行 `python flat_core.py --cases 120` 与 `python high_jets.py --cases 96 --out high-jets-results.json`。均为 12 个工作进程，令 OPENBLAS_NUM_THREADS 和 VECLIB_MAXIMUM_THREADS 取 1。center_certificate.py 只用标准库进行有理验证。r032_flat_jet_check.py 补充核验第34轮的势核。
* local-metal 使用 Apple Swift 6.4.2，要求当前目录内有 flat.metal。以 `swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" flat.swift -o flat` 编译，在共享编译缓存中运行。命令为 `./flat 65 96 main 0` 和 `./flat 65 144 coordinated 0`。用 `python verify.py flat-main.json cpu-main.json` 和对应的 coordinated 文件名进行独立求积。exact_witness.py 和 scaled_certificate.py 提供有理证书。
* remote-cpu 使用 Python 3.12.13、NumPy 2.5.4、SciPy 1.18.1。单线程 BLAS 环境中并行执行 15 份 `python overlap.py i --exponent 19 --order 64 --replicates 3 --out shard-i.json`，以 0 至 14 替换 i。certificate.py 复算有限区间证书。crosscheck.py 和 crosscheck_shallow.py 使用高精度求积。
* remote-metal 使用 Swift 6.4.2 编译 scan.swift，命名为 dual_scan，在含有 dual.metal 的目录运行 `./dual_scan 25 29 11 31 65 19 1024 1024 joint.json`。用 `python check_cpu.py joint.json joint.cpu.json` 核验选出的 12 个点。rational_rectangle.py 复算浅尾井的有理下界。

本机 CPU 计算的一阶与高阶批次有重叠。各设备时刻、整机与进程指标见对应摘要。
