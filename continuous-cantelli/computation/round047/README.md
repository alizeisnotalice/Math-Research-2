# 双机复现参数
在四个子目录的临时副本运行，编译使用 Apple Swift 6.4.2 的 -O 模式及已有共享模块缓存。GPU 校验脚本的两个参数依次为 GPU 输入 JSON 和 CPU 输出 JSON。
- 本机 CPU 采用单线程 BLAS，运行 python peano.py --cases 48，各进程用不同的确定性输入。
- 本机 GPU 将 quartic.swift 编译为 quartic，置于 quartic.metal 所在目录，运行 ./quartic 40 68 main 0。原版小扰动数值已发现消减误差。最小修复在 diagnostic 子目录，编译其中 quartic.swift，再运行 ./diagnostic 0 68 stable 984450415。
- MacBook CPU 在单线程 BLAS 环境并行启动 15 份 python centers.py i --exponent 14 --seconds 41 --out shard-i.json，以 0 至 14 替换 i。precision.py 和 cancellation.py 分别提供高精度数值与符号核对。
- MacBook GPU 编译 scan.swift 后在含 window.metal 的目录运行 ./window_scan 29 17 17 17 15 51 1536 4096 sync.json。check_cpu.py 使用适应性求积独立检查。
